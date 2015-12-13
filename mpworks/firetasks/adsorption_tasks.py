from monty.os.path import zpath

__author__ = 'weichen'


from fireworks.utilities.fw_serializers import FWSerializable
from fireworks.core.firework import FireTaskBase, FWAction
from fireworks.core.firework import Firework, Workflow
from mpworks.firetasks.vasp_io_tasks import VaspWriterTask, VaspToDBTask
from mpworks.firetasks.custodian_task import get_custodian_task
from fireworks.utilities.fw_utilities import get_slug
from pymatgen import Composition
from pymatgen.matproj.snl import StructureNL
from mpworks.workflows import snl_to_wf
from mpworks.firetasks.snl_tasks import AddSNLTask
from mpworks.snl_utils.mpsnl import MPStructureNL
from pymatgen.core.structure import Structure
from mpworks.workflows.wf_settings import QA_VASP, QA_DB, QA_VASP_SMALL
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints
from pymatgen.analysis.adsorption import reorient_z, AdsorbateSiteFinder
from pymatgen.core.surface import generate_all_slabs
#from pymatgen.analysis.adsorption import generate_adsorption_structures

def update_spec_adsorbate(spec, user_vasp_settings=None):
    """
    This method is intended to update the spec of the workflow
    with parameters particular to adsorbate calculates e. g.
    dipole corrections
    """
    # TODO: write this method
    # TODO: Dipole corrections
    # TODO: selective dynamics?
    fw_spec = spec
    update_set = {"ISIF" : 0, 
                  "LDIPOL" : "True",
                  "IDIPOL" : 3
                 }
    fw_spec['vasp']['incar'].update(update_set)
    return fw_spec

class SetupSlabTasks(FireTaskBase, FWSerializable):
    """
    This class sets up surface calculations for a slab calculation
    """
    _fw_name = "Setup Slab Tasks"

    def run_task(self, fw_spec):
        # Read structure from previous structural relaxation
        relaxed_struct = fw_spec['output']['crystal']
        # Generate slabs 
        # TODO: probably should make this so we can customize 
        #           parameters for vacuum, repeat etc.
        slabs_confined = generate_all_slabs(relaxed_struct, 1, 6.0, 0.0,
                                            max_normal_search = 1)
        slabs = generate_all_slabs(relaxed_struct, 1, 6.0, 10.0,
                                   max_normal_search = 1)
        wf=[]
        for i, slab in enumerate(slabs):
            fws=[]
            connections={}
            miller_string = ''.join([str(n) for n in slab.miller_index])
            f = '_'.join([Composition(slab.formula).alphabetical_formula, 
                          miller_string])
            snl = StructureNL(slab, 'Joseph Montoya <montoyjh@lbl.gov>', 
                              projects=["Adsorption"])
            tasks = [AddSNLTask()]
            snl_priority = fw_spec.get('priority', 1)
            spec = {'task_type': 'Add Slab to SNL database',
                    'snl': snl.as_dict(),
                    '_queueadapter': QA_DB,
                    '_priority': snl_priority}
            if 'snlgroup_id' in fw_spec and isinstance(snl, MPStructureNL):
                spec['force_mpsnl'] = snl.as_dict()
                spec['force_snlgroup_id'] = fw_spec['snlgroup_id']
                del spec['snl']
            fws.append(Firework(tasks, spec,
                                name=get_slug(f + '--' + spec['task_type']),
                                fw_id=-1000+i*10))
            connections[-1000+i*10] = [-999+i*10]
            spec = snl_to_wf._snl_to_spec(snl,
                                          parameters={'exact_structure':True})
            spec = update_spec_adsorbate(spec)
            spec['miller_index'] = slab.miller_index
            spec['original_task_id'] = fw_spec["task_id"]
            spec['_priority'] = fw_spec['_priority']*2
            #Turn off dupefinder for deformed structure
            del spec['_dupefinder']
            spec['task_type'] = "Optimize slab structure"
            fws.append(Firework([VaspWriterTask(),
                                 get_custodian_task(spec)],
                                spec, 
                                name=get_slug(f + '--' + spec['task_type']), 
                                fw_id=-999+i*10))
            
            priority = fw_spec['_priority']*3
            spec = {'task_type': 'VASP db insertion', 
                    '_priority': priority,
                    '_allow_fizzled_parents': True, 
                    '_queueadapter': QA_DB, 
                    'miller_index':slab.miller_index, 
                    'clean_task_doc':True,
                    'original_task_id':fw_spec["task_id"]}
            fws.append(Firework([VaspToDBTask()], 
                                spec, 
                                name=get_slug(f + '--' + spec['task_type']), 
                                fw_id=-998+i*10))

            spec = {'task_type': 'Setup adsorption task', 
                    '_priority': priority,
                    '_allow_fizzled_parents': True, 
                    '_queueadapter': QA_DB, 
                    'miller_index':slab.miller_index, 
                    'clean_task_doc':True,
                    'original_task_id':fw_spec["task_id"]}
            fws.append(Firework([SetupAdsorptionTask()],
                                spec,
                                name=get_slug(f + '--' + spec['task_type']),
                                fw_id=-997+i*10))
            connections[-999+i*10] = [-998+i*10]
            wf.append(Workflow(fws, connections))
        return FWAction(additions=wf)

class SetupAdsorptionTasks(FireTaskBase, FWSerializable):
    """
    This class sets up surface calculations for a slab calculation
    and adsorbates
    """
    _fw_name = "Setup Adsorption Tasks"

    def run_task(self, fw_spec):
        # Read structure from previous structural relaxation
        relaxed_struct = fw_spec['output']['crystal']
        # Generate slabs 
        # TODO: probably should make this so we can customize 
        #           parameters for vacuum, repeat etc.
        slabs = generate_all_slabs(relaxed_struct, 1, 6.0, 10.0,
                                   max_normal_search = 1)
        wf=[]
        i = 0
        for slab in slabs:
            asf = AdsorbateSiteFinder(slab)
            fws=[]
            connections={}
            miller_string = ''.join([str(n) for n in slab.miller_index])
            f = '_'.join([Composition(slab.formula).alphabetical_formula, 
                          miller_string])
            snl = StructureNL(asf.slab, 'Joseph Montoya <montoyjh@lbl.gov>', 
                              projects=["Adsorption"])
            tasks = [AddSNLTask()]
            snl_priority = fw_spec.get('priority', 1)
            spec = {'task_type': 'Add Slab to SNL database',
                    'snl': snl.as_dict(),
                    '_queueadapter': QA_DB,
                    '_priority': snl_priority}
            if 'snlgroup_id' in fw_spec and isinstance(snl, MPStructureNL):
                spec['force_mpsnl'] = snl.as_dict()
                spec['force_snlgroup_id'] = fw_spec['snlgroup_id']
                del spec['snl']
            fws.append(Firework(tasks, spec,
                                name=get_slug(f + '--' + spec['task_type']),
                                fw_id=-1000+i*10))
            connections[-1000+i*10] = [-999+i*10]
            spec = snl_to_wf._snl_to_spec(snl,
                                          parameters={'exact_structure':True})
            spec = update_spec_adsorbate(spec)
            spec['miller_index'] = slab.miller_index
            spec['original_task_id'] = fw_spec["task_id"]
            spec['_priority'] = fw_spec['_priority']*2
            #Turn off dupefinder for deformed structure
            del spec['_dupefinder']
            spec['task_type'] = "Optimize slab structure"
            fws.append(Firework([VaspWriterTask(),
                                 get_custodian_task(spec)],
                                 spec, 
                                 name=get_slug(f + '--' + spec['task_type']), 
                                 fw_id=-999+i*10))
            
            priority = fw_spec['_priority']*3
            spec = {'task_type': 'VASP db insertion', 
                    '_priority': priority,
                    '_allow_fizzled_parents': True, 
                    '_queueadapter': QA_DB, 
                    'miller_index':slab.miller_index, 
                    'clean_task_doc':True,
                    'original_task_id':fw_spec["task_id"]}
            fws.append(Firework([VaspToDBTask()], 
                                spec, 
                                name=get_slug(f + '--' + spec['task_type']), 
                                fw_id=-998+i*10))
            connections[-999+i*10] = [-998+i*10]
            # TODO: Hard code adsorbates for now, later make this mutable
            adsorbate_dict = {'H': [[0., 0., 0.]],
                               'O': [[0., 0., 0.]],
                               'OH':[[0.0, 0.0, 0.0], 
                                     [0.37, 0.13, 1.7]],
                               'OOH':[[0., 0., 0.], 
                                      [-1.06, -0.4, 0.8], 
                                      [-0.7, -0.3, 1.7]]
                        }
            for ads in adsorbate_dict:
                structs = asf.generate_adsorption_structures(ads, adsorbate_dict[ads], 
                                                repeat = [2, 2, 1])
                for struct in structs:
                    snl = StructureNL(struct, 'Joseph Montoya <montoyjh@lbl.gov>', 
                                      projects=["Adsorption"])
                    tasks = [AddSNLTask()]
                    snl_priority = fw_spec.get('priority', 1)
                    spec = {'task_type': 'Add Slab to SNL database',
                            'snl': snl.as_dict(),
                            '_queueadapter': QA_DB,
                            '_priority': snl_priority}
                    if 'snlgroup_id' in fw_spec and isinstance(snl, MPStructureNL):
                        spec['force_mpsnl'] = snl.as_dict()
                        spec['force_snlgroup_id'] = fw_spec['snlgroup_id']
                        del spec['snl']
                    fws.append(Firework(tasks, spec,
                                        name=get_slug(ads + '_' + f + '--' + spec['task_type']),
                                        fw_id=-997+i*10))
                    connections[-997+i*10] = [-996+i*10]
                    spec = snl_to_wf._snl_to_spec(snl,
                                                  parameters={'exact_structure':True})
                    spec = update_spec_adsorbate(spec)
                    spec['adsorbate'] = ads
                    spec['miller_index'] = slab.miller_index
                    spec['original_task_id'] = fw_spec["task_id"]
                    spec['_priority'] = fw_spec['_priority']*2
                    #Turn off dupefinder for deformed structure
                    del spec['_dupefinder']
                    spec['task_type'] = "Optimize adsorbate structure"
                    fws.append(Firework([VaspWriterTask(),
                                         get_custodian_task(spec)],
                                         spec, 
                                         name=get_slug(ads + '_' + f + '--' + spec['task_type']), 
                                         fw_id=-996+i*10))
                    connections[-996+i*10] = [-995+i*10]

                    priority = fw_spec['_priority']*3
                    spec = {'task_type': 'VASP db insertion', 
                            '_priority': priority,
                            '_allow_fizzled_parents': True, 
                            '_queueadapter': QA_DB, 
                            'miller_index':slab.miller_index, 
                            'adsorbate':ads,
                            'clean_task_doc':True,
                            'original_task_id':fw_spec["task_id"]}
                    fws.append(Firework([VaspToDBTask()], 
                                        spec, 
                                        name=get_slug(f + '--' + spec['task_type']), 
                                        fw_id=-995+i*10))
                    i += 1
            wf.append(Workflow(fws, connections))
        return FWAction(additions=wf)

if __name__ == "__main__":
    from fireworks.core.launchpad import LaunchPad
    lpad = LaunchPad.from_file('/global/u1/m/montoyjh/jhm_elastic/config/config_Mendel/my_launchpad.yaml')
    fw_spec = lpad.get_fw_by_id(2130).spec
    lpad_2 = LaunchPad.from_file('/global/u1/m/montoyjh/jhm_surfaces/config/config_Mendel/my_launchpad.yaml')
    fw = Firework([SetupAdsorptionTasks()],
             fw_spec,
             name = 'test Cu',
            )
    lpad_2.add_wf(fw)
    # 2128
