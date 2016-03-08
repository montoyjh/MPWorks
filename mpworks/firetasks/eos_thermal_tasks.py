from monty.os.path import zpath

__author__ = 'Cormac Toher'


from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.core.firework import FireTaskBase, FWAction
from pymatgen.io.vasp.inputs import Incar, Poscar
from pymatgen.analysis.eqn_of_state_thermal.eos_vasp_setup import ModifiedVolumeStructureSet
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
        return FWAction(additions=wf)        
