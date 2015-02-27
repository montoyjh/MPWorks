from mpworks.firetasks.vasp_io_tasks import VaspCopyTask, VaspWriterTask, VaspToDBTask
from mpworks.firetasks.dielectrics_tasks import SetupDFPTDielectricsTask

__author__ = 'Ioannis Petousis'

from fireworks.utilities.fw_utilities import get_slug
from mpworks.firetasks.custodian_task import get_custodian_task
from mpworks.snl_utils.mpsnl import get_meta_from_structure, MPStructureNL
from mpworks.firetasks.snl_tasks import AddSNLTask
from mpworks.workflows.wf_settings import QA_DB, QA_VASP, QA_CONTROL
from pymatgen import Composition
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio_set import MPStaticDielectricDFPTVaspInputSet, MPVaspInputSet
from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Kpoints
from fireworks.core.firework import Firework, Workflow
from mpworks.workflows import snl_to_wf



def snl_to_wf_static_dielectrics(snl, parameters=None):
    fws = []
    connections = {}
    parameters = parameters if parameters else {}

    snl_priority = parameters.get('priority', 1)
    priority = snl_priority * 2  # once we start a job, keep going!

    f = Composition(snl.structure.composition.reduced_formula).alphabetical_formula
    
    # add the SNL to the SNL DB and figure out duplicate group
    tasks = [AddSNLTask()]
    spec = {'task_type': 'Add to SNL database', 'snl': snl.as_dict(), '_queueadapter': QA_DB, '_priority': snl_priority}
    if 'snlgroup_id' in parameters and isinstance(snl, MPStructureNL):
        spec['static_dielectrics_mpsnl'] = snl.as_dict()
        spec['static_dielectrics_snlgroup_id'] = parameters['snlgroup_id']
        del spec['snl']
    fws.append(Firework(tasks, spec, name=get_slug(f + '--' + spec['task_type']), fw_id=0))
    

    if 'force_convergence' in snl.projects:
        # run optmization with force as convergence criterion:
        spec = snl_to_wf._snl_to_spec(snl, parameters=parameters)
        mpvis = MPVaspInputSet()
        incar = mpvis.get_incar(snl.structure)

        incar.update({"EDIFF":parameters['EDIFF']})
        incar.update({"EDIFFG":parameters['EDIFFG']})
        incar.update({"ENCUT":parameters['ENCUT']})
        spec['vasp']['incar'] = incar.as_dict()
        kpoints_density = int(parameters['kpoint_density'])
        k=Kpoints.automatic_density(snl.structure, kpoints_density)
        spec['vasp']['kpoints'] = k.as_dict()
        # spec = update_spec_static_dielectrics_convergence(spec)
        del spec['_dupefinder']
        # spec['run_tags'].append("origin")
        spec['_priority'] = priority
        spec['_queueadapter'] = QA_VASP
        spec['task_type'] = "force convergence" # Change name here: delete Vasp? 
        tasks = [VaspWriterTask(), get_custodian_task(spec)]
        fws.append(Firework(tasks, spec, name=get_slug(f + '--' + spec['task_type']), fw_id=1))
        connections[0] = [1] # define fw_id=1 is dependent on completion of fw_id=0

        # insert into DB - Force optimization
        spec = {'task_type': 'VASP Force db insertion', '_priority': priority, '_allow_fizzled_parents': True, '_queueadapter': QA_DB}
        fws.append(Firework([VaspToDBTask()], spec, name=get_slug(f + '--' + spec['task_type']), fw_id=2))
        connections[1] = [2] # define fw_id=2 is dependent on completion of fw_id=1

        spec= {'task_type': 'Setup DFPT Dielectrics Task', '_priority': priority, '_queueadapter': QA_CONTROL}
        fws.append(Firework([SetupDFPTDielectricsTask()], spec, name=get_slug(f + '--' + spec['task_type']), fw_id=3))
        connections[2] = [3]

        wf_meta = get_meta_from_structure(snl.structure)
        wf_meta['run_version'] = 'May 2013 (1)'

        if '_materialsproject' in snl.data and 'submission_id' in snl.data['_materialsproject']:
            wf_meta['submission_id'] = snl.data['_materialsproject']['submission_id']

        return Workflow(fws, connections, name=Composition(snl.structure.composition.reduced_formula).alphabetical_formula, metadata=wf_meta)



    # run DFPT for static dielectrics run:
    if 'force_convergence' in snl.projects:
        relaxed_structure = spec['output']['crystal']

    spec = snl_to_wf._snl_to_spec(snl, parameters=parameters)
    mpvis = MPStaticDielectricDFPTVaspInputSet()
    incar = mpvis.get_incar(snl.structure)

    if 'force_convergence' in snl.projects:
        spec['vasp']['poscar'] = relaxed_structure
        #poscar = mpvis.get_poscar(Structure.from_dict(spec['output']['crystal']))

    incar.update({"EDIFF":parameters['EDIFF']})
    incar.update({"ENCUT":parameters['ENCUT']})
    spec['vasp']['incar'] = incar.as_dict()
    kpoints_density = int(parameters['kpoint_density'])
    k=Kpoints.automatic_density(snl.structure, kpoints_density)
    spec['vasp']['kpoints'] = k.as_dict()
    # spec = update_spec_static_dielectrics_convergence(spec)
    del spec['_dupefinder']
    # spec['run_tags'].append("origin")
    spec['_priority'] = priority
    spec['_queueadapter'] = QA_VASP
    spec['task_type'] = "Static Dielectrics Calculation" # Change name here: delete Vasp? 
    tasks = [VaspWriterTask(), get_custodian_task(spec)]
    fws.append(Firework(tasks, spec, name=get_slug(f + '--' + spec['task_type']), fw_id=1))
    connections[0] = [1] # define fw_id=1 is dependent on completion of fw_id=0

    # insert into DB - Static Dielectrics run
    spec = {'task_type': 'VASP db insertion', '_priority': priority, '_allow_fizzled_parents': True, '_queueadapter': QA_DB}
    fws.append(Firework([VaspToDBTask()], spec, name=get_slug(f + '--' + spec['task_type']), fw_id=2))
    connections[1] = [2] # define fw_id=2 is dependent on completion of fw_id=1

    wf_meta = get_meta_from_structure(snl.structure)
    wf_meta['run_version'] = 'May 2013 (1)'

    if '_materialsproject' in snl.data and 'submission_id' in snl.data['_materialsproject']:
        wf_meta['submission_id'] = snl.data['_materialsproject']['submission_id']

    return Workflow(fws, connections, name=Composition(snl.structure.composition.reduced_formula).alphabetical_formula, metadata=wf_meta)
