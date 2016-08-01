from monty.os.path import zpath
import os
import json
from pymongo import MongoClient
import numpy as np
from decimal import Decimal

__author__ = 'Cormac Toher'


from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.core.firework import FireTaskBase, FWAction
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.analysis.eqn_of_state_thermal.eos_vasp_setup import ModifiedVolumeStructureSet
from pymatgen.analysis.eqn_of_state_thermal.eos_wf_run import eos_thermal_properties
from fireworks.core.firework import Firework, Workflow
from mpworks.firetasks.vasp_io_tasks import VaspWriterTask, VaspToDBTask
from mpworks.firetasks.custodian_task import get_custodian_task
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Composition
from pymatgen.matproj.snl import StructureNL
from mpworks.workflows import snl_to_wf
from mpworks.firetasks.snl_tasks import AddSNLTask
from mpworks.snl_utils.mpsnl import get_meta_from_structure, MPStructureNL
from pymatgen.core.structure import Structure
from mpworks.workflows.wf_settings import QA_VASP, QA_DB, QA_VASP_SMALL
from pymatgen.io.vasp.inputs import Poscar, Kpoints
from pymatgen import MPRester

def update_spec_force_convergence(spec, user_vasp_settings=None):
    fw_spec = spec
    update_set = {"ENCUT": 700, "EDIFF": 0.000001, "ALGO":"N", "NPAR":2, "NSW":0}
    if user_vasp_settings and user_vasp_settings.get("incar"):
            update_set.update(user_vasp_settings["incar"])
    fw_spec['vasp']['incar'].update(update_set)
    old_struct=Poscar.from_dict(fw_spec["vasp"]["poscar"]).structure
    if user_vasp_settings and user_vasp_settings.get("kpoints"):
        kpoints_density = user_vasp_settings["kpoints"]["kpoints_density"]
    else:
        kpoints_density = 7000
    k=Kpoints.automatic_density(old_struct, kpoints_density)
    fw_spec['vasp']['kpoints'] = k.as_dict()
    return fw_spec


class SetupFConvergenceTask(FireTaskBase, FWSerializable):
    _fw_name = "Setup Force Convergence Task"

    def run_task(self, fw_spec):
        incar = fw_spec['vasp']['incar']
        update_set = {"ENCUT": 600, "EDIFF": 0.00005}
        incar.update(update_set)
        kpoints = fw_spec['vasp']['kpoints']
        k = [int(round(2.5*k)) if int(round(2.5*k))%2 else int(round(2.5*k))+1 for k in kpoints['kpoints'][0]]
        kpoints['kpoints'] = [k]
        return FWAction()

class SetupEoSThermalTask(FireTaskBase, FWSerializable):
    _fw_name = "Setup EoS Thermal Task"

    def run_task(self, fw_spec):
        incar = Incar.from_file(zpath("INCAR"))
        incar.update({"NSW": 0})
        incar.update({"LCHARG": True})
        incar.update({"IBRION": -1})
        incar.write_file("INCAR")
        return FWAction()

class SetupModifiedVolumeStructTask(FireTaskBase, FWSerializable):
    _fw_name = "Setup Modified Volume Struct Task"

    def run_task(self, fw_spec):
        # Read structure from previous relaxation
        relaxed_struct = fw_spec['output']['crystal']
        modified_struct_set = ModifiedVolumeStructureSet(relaxed_struct)
        poisson_val = fw_spec['poisson_ratio']
        wf=[]
        for i, mod_struct in enumerate(modified_struct_set.modvol_structs):
            fws=[]
            connections={}
            f = Composition(mod_struct.formula).alphabetical_formula
            snl = StructureNL(mod_struct, 'Cormac Toher <cormac.toher@duke.edu>',projects=["Thermal"])
            tasks = [AddSNLTask()]
            snl_priority = fw_spec.get('priority', 1)
            spec = {'task_type': 'Add Modified Struct to SNL database', 'snl': snl.as_dict(),
                    '_queueadapter': QA_DB, '_priority': snl_priority}
            if 'snlgroup_id' in fw_spec and isinstance(snl, MPStructureNL):
                spec['force_mpsnl'] = snl.as_dict()
                spec['force_snlgroup_id'] = fw_spec['snlgroup_id']
                del spec['snl']
            fws.append(Firework(tasks, spec, name=get_slug(f + '--' + spec['task_type']), fw_id=-1000+i*10))
            connections[-1000+i*10] = [-999+i*10]
            spec = snl_to_wf._snl_to_spec(snl, parameters={'exact_structure':True})
            spec = update_spec_force_convergence(spec)
            spec['strainfactor'] = modified_struct_set.strainfactors[i]
            spec['original_task_id']=fw_spec["task_id"]
            spec['_priority'] = fw_spec['_priority']*2
            #Turn off dupefinder for modified structure
            del spec['_dupefinder']
            spec['task_type'] = "Calculate static modified structure"
            fws.append(Firework([VaspWriterTask(), SetupEoSThermalTask(), get_custodian_task(spec)],
                                spec, name=get_slug(f + '--' + fw_spec['task_type']), fw_id=-999+i*10))

            priority = fw_spec['_priority']*3
            spec = {'task_type': 'VASP db insertion', '_priority': priority,
                    '_allow_fizzled_parents': True, '_queueadapter': QA_DB, 
                    'eqn_of_state_thermal':"modified_structure", 'clean_task_doc':True,
                    'strainfactor':modified_struct_set.strainfactors[i], 'original_task_id':fw_spec["task_id"]}
            fws.append(Firework([VaspToDBTask()], spec, name=get_slug(f + '--' + spec['task_type']), fw_id=-998+i*10))
            connections[-999+i*10] = [-998+i*10]
            wf.append(Workflow(fws, connections))
#        fws=[]
#        connections={}
#        spec = {'task_type': 'Add EoS Thermal Data to DB Task', '_priority': priority,
#         '_queueadapter': QA_DB, 'poisson_ratio': poisson_val}
#        print("Calling firetask")
#        fws.append(Firework([AddEoSThermalDataToDBTask()], spec,
#                        name=get_slug(f + '--' + spec['task_type']),fw_id=-997+i*10))
#        connections[-998+i*10] = [-997+i*10]
#        wf.append(Workflow(fws, connections))
        return FWAction(additions=wf)        


class AddEoSThermalDataToDBTask(FireTaskBase, FWSerializable):
    _fw_name = "Add EoS Thermal Data to DB Task"

    def run_task(self, fw_spec):
	print("Running EoS Thermal Data")
        db_dir = os.environ['DB_LOC']
        db_path = os.path.join(db_dir, 'tasks_db.json')
        print("fw_spec keys = ", fw_spec.keys())        
	poisson_val = fw_spec['poisson_ratio']
        print("Poisson ratio = ", poisson_val)
        i = fw_spec['original_task_id']
#        i = fw_spec['task_id']

        with open(db_path) as f:
            db_creds = json.load(f)
        connection = MongoClient(db_creds['host'], db_creds['port'])
        tdb = connection[db_creds['database']]
        tdb.authenticate(db_creds['admin_user'], db_creds['admin_password'])
        tasks = tdb[db_creds['collection']]
        eos_thermal = tdb['eos_thermal']
#        i = 'mp-1265'
        ndocs = tasks.find({"original_task_id": i, 
                            "state":"successful"}).count()
#        ndocs = tasks.find({"task_id": i, 
#                            "state":"successful"}).count()
        existing_doc = eos_thermal.find_one({"relaxation_task_id" : i})
#        existing_doc = eos_thermal.find_one({"task_id" : i})
        if existing_doc:
            print "Updating: " + i
        else:
            print "New material: " + i
        d = {"analysis": {}, "error": [], "warning": []}
        d["ndocs"] = ndocs
        o = tasks.find_one({"task_id" : i},
                           {"pretty_formula" : 1, "spacegroup" : 1,
                            "snl" : 1, "snl_final" : 1, "run_tags" : 1})
        if not o:
            raise ValueError("Cannot find original task id")
        # Get energy vs. volume from strained structures
        d["strain_tasks"] = {}
        #ss_dict = {}
        volume_values = []
        energy_values = []
        #m = MPRester("o95HBW1KgfiSA0tN")	
        #elasticity_data = m.query(criteria={"task_id": "mp-1265"}, properties=["elasticity"])
        #poisson_val = elasticity_data[0]["elasticity"]["poisson_ratio"]
#        for k in tasks.find({"original_task_id": i}, 
        for k in tasks.find({"task_id": i}, 
                            {"strainfactor":1,
                             "calculations.output":1,
                             "state":1, "task_id":1}):
            defo = k['strainfactor']
            energy_values.append(item['output.final_energy'])
            volume_values.append(item['output.crystal.lattice.volume'])
            # d_ind = np.nonzero(defo - np.eye(3))
            # delta = Decimal((defo - np.eye(3))[d_ind][0])
            # Normal deformation
            # if d_ind[0] == d_ind[1]:
            #    dtype = "_".join(["d", str(d_ind[0][0]), 
            #                      "{:.0e}".format(delta)])
            # Shear deformation
            # else:
            #    dtype = "_".join(["s", str(d_ind[0] + d_ind[1]),
            #                      "{:.0e}".format(delta)])
            #sm = IndependentStrain(defo)
            #d["strain_tasks"][dtype] = {"state" : k["state"],
            #                                 "strainfactor" : defo,
            #                                 "strain" : sm.tolist(),
            #                                 "task_id": k["task_id"]}
            #if k["state"] == "successful":
            #    st = Stress(k["calculations"][-1]["output"] \
            #                ["ionic_steps"][-1]["stress"])
            #    ss_dict[sm] = st
        d["snl"] = o["snl"]
        if "run_tags" in o.keys():
            d["run_tags"] = o["run_tags"]
            for tag in o["run_tags"]:
                if isinstance(tag, dict):
                    if "input_id" in tag.keys():
                        d["input_mp_id"] = tag["input_id"]
        d["snl_final"] = o["snl_final"]
        d["pretty_formula"] = o["pretty_formula"]

        # Old input mp-id style
        if o["snl"]["about"].get("_mp_id"):
            d["material_id"] = o["snl"]["about"]["_mp_id"]

        # New style
        elif "input_mp_id" in d:
            d["material_id"] = d["input_mp_id"]
        else:
            d["material_id"] = None
        d["relaxation_task_id"] = i

        calc_struct = Structure.from_dict(o["snl_final"])
        # TODO:
        # JHM: This test is unnecessary at the moment, but should be redone
        """
        conventional = is_conventional(calc_struct)
        if conventional:
            d["analysis"]["is_conventional"] = True
        else:
            d["analysis"]["is_conventional"] = False
        """
        d["spacegroup"]=o.get("spacegroup", "Unknown")

        
        if ndocs >= 21:
            eos_thermal_dict = {}
            # Perform thermal equation of state fitting and analysis
#            eos_thermal_dict = eos_thermal_properties.eos_thermal_run(calc_struct, volume_values, energy_values, ieos=2)
            eos_thermal_dict = eos_thermal_properties.eos_thermal_run(calc_struct, volume_values, energy_values, ieos=2, idebye=0, poissonratio=poisson_val)
            
	    # Test to check if results have been calculated
	    print("Thermal conductivity = ", eos_thermal_dict["thermal_conductivity"])
            # Add equation of state results to dict
            d["temperature"] = eos_thermal_dict["temperature"]
            d["pressure"] = eos_thermal_dict["pressure"]
            d["Thermal_conductivity_temp_list"] = eos_thermal_dict["thermal_conductivity"]
            d["Debye_temperature_temp_list"] = eos_thermal_dict["Debye_temperature"]
            d["Gruneisen_parameter_temp_list"] = eos_thermal_dict["Gruneisen_parameter"]
            d["Heat_capacity_Cv_temp_list"] = eos_thermal_dict["Heat_capacity_Cv"]
            d["Heat_capacity_Cp_temp_list"] = eos_thermal_dict["Heat_capacity_Cp"]
            d["Volume_temp_list"] = eos_thermal_dict["Volume"]
            d["Bulk_modulus_temp_list"] = eos_thermal_dict["Bulk_modulus"]
            d["BM_coeffs"] = eos_thermal_dict["BM_coeffs"]
            jtdbest = eos_thermal_dict["Best_fit_temperature"]

            # Add values at specific temperatures (Debye temperature and Gruneisen parameter best fit, values at 300K)
            d["Debye_temperature"] =  d["Debye_temperature_temp_list"][eos_thermal_dict["Best_fit_temperature"]]
            d["Gruneisen_parameter"] = d["Gruneisen_parameter_temp_list"][eos_thermal_dict["Best_fit_temperature"]]
            d["Thermal_conductivity_300K"] = d["Thermal_conductivity_temp_list"][eos_thermal_dict["300K_point"]]
            d["Heat_capacity_Cv_300K"] = d["Heat_capacity_Cv_temp_list"][eos_thermal_dict["300K_point"]]
            d["Heat_capacity_Cp_300K"] = d["Heat_capacity_Cp_temp_list"][eos_thermal_dict["300K_point"]]
            d["Bulk_modulus_300K"] = d["Bulk_modulus_temp_list"][eos_thermal_dict["300K_point"]]



        else:
            d['state'] = "Fewer than 21 successful tasks completed"
            return FWAction()

        if o["snl"]["about"].get("_kpoint_density"):
            d["kpoint_density"]= o["snl"]["about"].get("_kpoint_density")

        if d["error"]:
            raise ValueError("Thermal equation of state analysis failed: {}".format(d["error"]))
        elif d["analysis"]["filter_pass"]:
            d["state"] = "successful"
        else:
            d["state"] = "filter_failed"
        eos_thermal.update({"relaxation_task_id": d["relaxation_task_id"]}, 
                           d, upsert=True)
        return FWAction()
