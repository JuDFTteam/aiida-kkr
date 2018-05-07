=========
Workflows
=========

This page can contain a simple tutorial for your code.


Density of states
+++++++++++++++++

Workflow: ``aiida_kkr.workflows.dos``

::

    from aiida_kkr.workflows.dos import kkr_dos_wc
    kkr_dos_wc.get_wf_defaults()
    kkr_dos_wc.get_inputs_template()
    
    Version of workflow: 0.3
    
    {'dos_params': {'emax': 1, 'tempr': 200, 'emin': -1, 'kmesh': [50, 50, 50], 'nepts': 61}, 
     'queue_name': '', 'walltime_sec': 3600, 'custom_scheduler_commands': '', 
     'resources': {'num_machines': 1}, 'use_mpi': False}
     
     _WorkChainSpecInputs({'_label': None, '_description': None, '_store_provenance': True, 
                           'dynamic': None, 'kkr': None, 'remote_data': None, 
                           'wf_parameters': <ParameterData: uuid: b65d1be1-e033-43e4-badc-9fe00335dace (unstored)>})

::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    
::

    from aiida.orm import load_node
    kkr_remote_folder = load_node(22852).out.remote_folder
    from aiida.orm import DataFactory
    ParameterData = DataFactory('parameter')
    workflow_settings = ParameterData(dict={'dos_params':{'emax': 1, 'tempr': 200, 'emin': -1,
                                                          'kmesh': [20, 20, 20], 'nepts': 81}})
    
::

    from aiida.work import run
    run(kkr_dos_wc, _label='test_doscal', _description='My test dos calculation with some 
                                                        description.', 
        kkr=kkrcode, remote_data=kkr_remote_folder, wf_parameters=workflow_settings)
    
    {'dos_data': <XyData: uuid: f4ef38ec-d372-46ce-9f2d-deeb5d1cfbe1 (pk: 22905)>,
     'dos_data_interpol': <XyData: uuid: 0617cb4f-722d-4bc0-a324-6fc39500bd52 (pk: 22906)>,
     'results_wf': <ParameterData: uuid: 12a2a350-be5b-478e-acf9-58e3f7b1054a (pk: 22904)>}
    
Generate KKR start potential
++++++++++++++++++++++++++++

Workflow: ``aiida_kkr.workflows.voro_start``

::

    from aiida_kkr.workflows.voro_start import kkr_startpot_wc
    kkr_startpot_wc.get_wf_defaults()
    kkr_startpot_wc.get_inputs_template()
    
    Version of workflow: 0.5
    
    {'delta_e_min': 1.0, 'delta_e_min_core_states': 1.0, 'fac_cls_increase': 1.3, 
     'walltime_sec': 3600, 'custom_scheduler_commands': '', 'check_dos': True, 
     'use_mpi': False, 'natom_in_cls_min': 79, 'dos_params': {'emax': 1, 'tempr': 200, 
                                                              'emin': -1, 'kmesh': [50, 50, 50], 
                                                              'nepts': 61}, 
     'queue_name': '', 'num_rerun': 4, 'threshold_dos_zero': 0.001, 'r_cls': 1.3, 
     'resources': {'num_machines': 1}}
    
     _WorkChainSpecInputs({'_label': None, '_description': None, '_store_provenance': True, 
                          'dynamic': None, 'calc_parameters': None, 'kkr': None, 'voronoi': None, 
                          'wf_parameters': <ParameterData: uuid: aaa6f37f-5293-4801-a0bb-672af08b6870 (unstored)>, 
                          'structure': None})

::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    vorocode = Code.get_from_string('voronoi@my_mac')
    
::

    from aiida_kkr.tools.kkr_params import kkrparams
    kkr_settings = kkrparams(NSPIN=1, LMAX=2)
    
::

    # create Copper bulk aiida Structure
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')
    
::

    from aiida.work import run
    ParameterData = DataFactory('parameter')
    run(kkr_startpot_wc, structure=Cu, voronoi=vorocode, kkr=kkrcode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()))

::

    {'last_doscal_dosdata': <XyData: uuid: d210b6fd-e415-44e2-be3d-54d16ac3c12d (pk: 22627)>,
     'last_doscal_dosdata_interpol': <XyData: uuid: 8420c34c-2a52-446b-8da1-38d679201afa (pk: 22628)>,
     'last_doscal_results': <ParameterData: uuid: 087efc99-f675-46d8-a577-5a7bc98144d4 (pk: 22626)>,
     'last_params_voronoi': <ParameterData: uuid: f9e65da2-39e4-493b-ae79-88664031be50 (pk: 22610)>,
     'last_voronoi_remote': <RemoteData: uuid: 8df1c73f-3715-4e57-9c72-67416d99bff4 (pk: 22612)>,
     'last_voronoi_results': <ParameterData: uuid: 2c8b440e-965d-4bff-9d64-2f6f829b1694 (pk: 22614)>,
     'results_vorostart_wc': <ParameterData: uuid: bf6e2a95-55e0-424b-a3ec-a7e821f08a34 (pk: 22629)>}



KKR scf cycle
+++++++++++++

Workflow: ``aiida_kkr.workflows.kkr_scf``

::

    from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
    print kkr_scf_wc.get_wf_defaults()
    print kkr_scf_wc.get_inputs_template()
    
    Version of workflow: 0.6
    
    {'strmix': 0.03, 'brymix': 0.05, 'init_pos': None, 'convergence_criterion': 1e-08, 
     'custom_scheduler_commands': '', 'convergence_setting_coarse': {'npol': 7, 'tempr': 1000.0, 
                                                                     'n1': 3, 'n2': 11, 'n3': 3, 
                                                                     'kmesh': [10, 10, 10]}, 
     'mixreduce': 0.5, 'mag_init': False, 'retreive_dos_data_scf_run': False, 
     'dos_params': {'emax': 0.6, 'tempr': 200, 'nepts': 81, 'kmesh': [40, 40, 40], 'emin': -1}, 
     'hfield': 0.02, 'queue_name': '', 'threshold_aggressive_mixing': 0.008, 
     'convergence_setting_fine': {'npol': 5, 'tempr': 600.0, 'n1': 7, 'n2': 29, 'n3': 7, 
                                  'kmesh': [30, 30, 30]}, 
     'use_mpi': False, 'nsteps': 50, 'resources': {'num_machines': 1}, 'delta_e_min': 1.0, 
     'walltime_sec': 3600, 'check_dos': True, 'threshold_switch_high_accuracy': 0.001, 
     'kkr_runmax': 5, 'threshold_dos_zero': 0.001}

    _WorkChainSpecInputs({'_label': None, '_description': None, '_store_provenance': True, 
                          'dynamic': None, 'calc_parameters': None, 'kkr': None, 'voronoi': None, 
                          'remote_data': None, 'wf_parameters': <ParameterData: uuid: b132dfc4-3b7c-42e7-af27-4083802aff40 (unstored)>, 
                          'structure': None})
    
::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    vorocode = Code.get_from_string('voronoi@my_mac')
    
::

    from aiida_kkr.tools.kkr_params import kkrparams
    kkr_settings = kkrparams(NSPIN=1, LMAX=2)
    
::

    # create Copper bulk aiida Structure
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')
    
::
    
    from aiida.orm import load_node
    kkr_startpot = load_node(22586)
    last_vorono_remote = kkr_startpot.get_outputs_dict().get('last_voronoi_remote')
    
::

    from aiida.work import run
    ParameterData = DataFactory('parameter')
    run(kkr_scf_wc, kkr=kkrcode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()), remote_data=last_vorono_remote)
    
::

    run(kkr_scf_wc, structure=Cu, kkr=kkrcode, voronoi=vorocode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()))
    
    {'final_dosdata_interpol': <XyData: uuid: 0c14146d-90aa-4eb8-834d-74a706e500bb (pk: 22872)>,
     'last_InputParameters': <ParameterData: uuid: 28a277ad-8998-4728-8296-75fd3b0c4eb4 (pk: 22875)>,
     'last_RemoteData': <RemoteData: uuid: d24cdfc1-938a-4308-b273-e0aa8697c975 (pk: 22876)>,
     'last_calc_out': <ParameterData: uuid: 1c8fab2d-e596-4874-9516-c1387bf7db7c (pk: 22874)>,
     'output_kkr_scf_wc_ParameterResults': <ParameterData: uuid: 0f21ac18-e556-49f8-aa26-55260d013fac (pk: 22878)>,
     'results_vorostart': <ParameterData: uuid: 93831550-8775-493a-907b-27a470b52dc8 (pk: 22877)>,
     'starting_dosdata_interpol': <XyData: uuid: 54fa57ad-f559-4837-ba1e-7db4ed67d5b0 (pk: 22873)>}


Equation of states
++++++++++++++++++

Workflow: ``aiida_kkr.workflows.eos``

.. warning:: Not implemented yet!


Check KKR parameter convergence
+++++++++++++++++++++++++++++++

``aiida_kkr.workflows.check_para_convergence``

.. warning:: Not implemented yet!
   

Find magnetic ground state
++++++++++++++++++++++++++

``aiida_kkr.workflows.check_magnetic_state``

.. warning:: Not implemented yet!