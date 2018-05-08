=========
Workflows
=========

This page can contain a short introduction to the workflows provided by ``aiida-kkr``.


Density of states
+++++++++++++++++

The density of states (DOS) workflow ``kkr_dos_wc`` automatically sets the right parameters in the 
input of a KKR calculation to perform a DOS calculation. The specifics of the DOS 
energy contour are set via the ``wf_parameters`` input node which contains default values
if no user input is given.

.. note::
    The default values of the ``wf_parameters`` input node can be extraced using 
    ``kkr_dos_wc.get_wf_defaults()``.

Inputs:
    * ``kkr`` (*aiida.orm.Code*): KKRcode using the ``kkr.kkr`` plugin
    * ``remote_data`` (*RemoteData*): The remote folder of the (converged) calculation whose output potential is used as input for the DOS run
    * ``wf_parameters`` (*ParameterData*, optional): Some settings of the workflow behavior (e.g. number of energy points in DOS contour etc.)
    * ``_label`` (*str*, optional): Label of the workflow
    * ``_description`` (*str*, optional): Longer description of the workflow
    
Returns:
    * ``dos_data`` (*XyData*): The DOS data on the DOS energy contour (i.e. at some finite temperature)
    * ``dos_data_interpol`` (*XyData*): The interpolated DOS from the line parallel to the real axis down onto the real axis
    * ``results_wf`` (*ParameterData*): The output node of the workflow containing some information on the DOS run

.. note::   
    The *x* and *y* arrays of the ``dos_data`` output nodes can easily be accessed using::
    
        x = dos_data_node.get_x()
        y = dos_data_node.get_y()
    
    where the returned list is of the form ``[label, numpy-array-of-data, unit]`` and the 
    *y*-array contains entries for total DOS, s-, p-, d-, ... contributions to the DOS, e.g.::
        
        [(u'interpolated dos tot', array([[...]]), u'states/eV'),
         (u'interpolated dos s', array([[...]]), u'states/eV'),
         (u'interpolated dos p', array([[...]]), u'states/eV'),
         (u'interpolated dos d', array([[...]]), u'states/eV'),
         (u'interpolated dos ns', array([[...]]), u'states/eV')]
    
                                        
Example Usage
-------------

We start by getting an installation of the KKRcode::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    
Next load the remote folder node of the previous calculation 
(here the :ref:`converged calculation of the Cu bulk test case <KKR_KKR_scf>`) 
from which we want to start the following DOS calculation::

    # import old KKR remote folder
    from aiida.orm import load_node
    kkr_remote_folder = load_node(22852).out.remote_folder
    
Then we set some settings of the workflow parameters (this step is optional)::

    # create workflow settings
    from aiida.orm import DataFactory
    ParameterData = DataFactory('parameter')
    workflow_settings = ParameterData(dict={'dos_params':{'emax': 1, 'tempr': 200, 'emin': -1,
                                                          'kmesh': [20, 20, 20], 'nepts': 81}})
    
Finally we run the workflow::

    from aiida_kkr.workflows.dos import kkr_dos_wc
    from aiida.work import run
    run(kkr_dos_wc, _label='test_doscal', _description='My test dos calculation.', 
        kkr=kkrcode, remote_data=kkr_remote_folder, wf_parameters=workflow_settings)
    
The following script can be used to plot the total interpolated DOS (in ``dos_data_interpol``) of the calculation above::

    x = dos_data_interpol.get_x()
    y = dos_data_interpol.get_y()[0] # tot only
    
    figure()
    plot(x[1][0], y[1][0], 'x-', lw=2, ms=6)
    xlabel(x[0]+' ({})'.format(x[-1]))
    ylabel(y[0]+' ({})'.format(y[-1]))
    axhline(0, color='grey')
    axvline(0, color='grey')

which will produce the following plot:

.. image:: ../images/DOS_Cu_example.png
    :width: 60%
    
    
Generate KKR start potential
++++++++++++++++++++++++++++

Workflow: ``kkr_startpot_wc``

Inputs:

::
        
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

Outputs:

::

    {'last_doscal_dosdata': <XyData: uuid: d210b6fd-e415-44e2-be3d-54d16ac3c12d (pk: 22627)>,
     'last_doscal_dosdata_interpol': <XyData: uuid: 8420c34c-2a52-446b-8da1-38d679201afa (pk: 22628)>,
     'last_doscal_results': <ParameterData: uuid: 087efc99-f675-46d8-a577-5a7bc98144d4 (pk: 22626)>,
     'last_params_voronoi': <ParameterData: uuid: f9e65da2-39e4-493b-ae79-88664031be50 (pk: 22610)>,
     'last_voronoi_remote': <RemoteData: uuid: 8df1c73f-3715-4e57-9c72-67416d99bff4 (pk: 22612)>,
     'last_voronoi_results': <ParameterData: uuid: 2c8b440e-965d-4bff-9d64-2f6f829b1694 (pk: 22614)>,
     'results_vorostart_wc': <ParameterData: uuid: bf6e2a95-55e0-424b-a3ec-a7e821f08a34 (pk: 22629)>}


                          
Example Usage
-------------

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

    from aiida_kkr.workflows.voro_start import kkr_startpot_wc
    from aiida.work import run
    ParameterData = DataFactory('parameter')
    run(kkr_startpot_wc, structure=Cu, voronoi=vorocode, kkr=kkrcode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()))

    

KKR scf cycle
+++++++++++++

Workflow: ``kkr_scf_wc``

Inputs:

::

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
    
Outputs:

::

    {'final_dosdata_interpol': <XyData: uuid: 0c14146d-90aa-4eb8-834d-74a706e500bb (pk: 22872)>,
     'last_InputParameters': <ParameterData: uuid: 28a277ad-8998-4728-8296-75fd3b0c4eb4 (pk: 22875)>,
     'last_RemoteData': <RemoteData: uuid: d24cdfc1-938a-4308-b273-e0aa8697c975 (pk: 22876)>,
     'last_calc_out': <ParameterData: uuid: 1c8fab2d-e596-4874-9516-c1387bf7db7c (pk: 22874)>,
     'output_kkr_scf_wc_ParameterResults': <ParameterData: uuid: 0f21ac18-e556-49f8-aa26-55260d013fac (pk: 22878)>,
     'results_vorostart': <ParameterData: uuid: 93831550-8775-493a-907b-27a470b52dc8 (pk: 22877)>,
     'starting_dosdata_interpol': <XyData: uuid: 54fa57ad-f559-4837-ba1e-7db4ed67d5b0 (pk: 22873)>}

          
Example Usage
-------------

Case 1: Start from structure and run voronoi calculation first 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    vorocode = Code.get_from_string('voronoi@my_mac')
    
::

    from aiida_kkr.tools.kkr_params import kkrparams
    kkr_settings = kkrparams(NSPIN=1, LMAX=2)
    
::
    
    from aiida.orm import load_node
    kkr_startpot = load_node(22586)
    last_vorono_remote = kkr_startpot.get_outputs_dict().get('last_voronoi_remote')
    
::

    from aiida_kkr.workflows.kkr_scf import kkr_scf_wc
    from aiida.work import run
    ParameterData = DataFactory('parameter')
    run(kkr_scf_wc, kkr=kkrcode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()), remote_data=last_vorono_remote)
    

Case 2: Start from previous KKR calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    # create Copper bulk aiida Structure
    fro numpy import array
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')
    
::

    run(kkr_scf_wc, structure=Cu, kkr=kkrcode, voronoi=vorocode, calc_parameters=ParameterData(dict=kkr_settings.get_dict()))
    

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