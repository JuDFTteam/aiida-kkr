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
    
Returns nodes:
    * ``dos_data`` (*XyData*): The DOS data on the DOS energy contour (i.e. at some finite temperature)
    * ``dos_data_interpol`` (*XyData*): The interpolated DOS from the line parallel to the real axis down onto the real axis
    * ``results_wf`` (*ParameterData*): The output node of the workflow containing some information on the DOS run

.. note::   
    The *x* and *y* arrays of the ``dos_data`` output nodes can easily be accessed using::
    
        x = dos_data_node.get_x()
        y = dos_data_node.get_y()
    
    where the returned list is of the form ``[label, numpy-array-of-data, unit]`` and the 
    *y*-array contains entries for total DOS, s-, p-, d-, ..., and non-spherical contributions to the DOS, e.g.::
        
        [(u'interpolated dos tot', array([[...]]), u'states/eV'),
         (u'interpolated dos s', array([[...]]), u'states/eV'),
         (u'interpolated dos p', array([[...]]), u'states/eV'),
         (u'interpolated dos d', array([[...]]), u'states/eV'),
         (u'interpolated dos ns', array([[...]]), u'states/eV')]
                                        
    Note that the output data are 2D arrays containing the atom resolved DOS, i.e. the DOS values for all atoms in the unit cell.
    
                                        
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
    
The following script can be used to plot the total interpolated DOS (in the 
``dos_data_interpol`` output node that can for example be access using 
``dos_data_interpol = <kkr_dos_wc-node>.out.dos_data_interpol`` where 
``<kkr_dos_wc-node>`` is the workflow node) of the calculation above::

    def plot_dos(dos_data_node):
        x = dos_data_node.get_x()
        y_all = dos_data_node.get_y()
        
        from matplotlib.pylab import figure, xlabel, ylabel, axhline, axvline, plot, legend, title
        
        figure()
        
        # loop over contributions (tot, s, p, d, ns)
        for y in y_all:
            if y==y_all[0]: # special line formatting for total DOS
                style = 'x-'
                lw = 3
            else:
                style = '--'
                lw = 2
            plot(x[1][0], y[1][0], style, lw=lw, ms=6, label=y[0].split('dos ')[1])
        
        # add axis labels etc                                
        xlabel(x[0]+' ({})'.format(x[-1]))                   
        ylabel(y[0].replace(' ns','')+' ({})'.format(y[-1]))
        axhline(0, color='grey', linestyle='dotted', zorder=-100)
        axvline(0, color='grey', linestyle='dotted', zorder=-100)
        legend(loc=2)
        title('DOS of bulk Cu')
    
    plot_dos(dos_data_interpol)
    
which will produce the following plot:

.. image:: ../images/DOS_Cu_example.png
    :width: 60%
    
    
Generate KKR start potential
++++++++++++++++++++++++++++

Workflow: ``kkr_startpot_wc``

Inputs:
    * ``structure`` (*StructureData*): 
    * ``voronoi`` (*Code*): 
    * ``kkr`` (*Code*): 
    * ``wf_parameters`` (*ParameterData*, optional): 
    * ``calc_parameters`` (*ParameterData*, optional): 
    * ``_label`` (*str*, optional): 
    * ``_description`` (*str*, optional): 
    
.. note::
    The default values of the ``wf_parameters`` input node can be extraced using 
    ``kkr_dos_wc.get_wf_defaults()`` and it should contain the following entries:
    
    General settings:
        * ``r_cls`` (*float*): 
        * ``natom_in_cls_min`` (*int*): 
        * ``fac_cls_increase`` (*float*): 
        * ``num_rerun`` (*int*): 
        
    Computer settings:
        * ``walltime_sec`` (*int*): 
        * ``custom_scheduler_commands`` (*str*): 
        * ``use_mpi`` (*bool*): 
        * ``queue_name`` (*str*): 
        * ``resources`` (*dict*): ``{'num_machines': 1}``
        
    Settings for DOS check of starting potential:
        * ``check_dos`` (*bool*): 
        * ``threshold_dos_zero`` (*float*): 
        * ``delta_e_min`` (*float*): 
        * ``delta_e_min_core_states`` (*float*): 
        * ``dos_params`` (*dict*): with the keys
            * ``emax`` (*float*): 
            * ``tempr`` (*float*): 
            * ``emin`` (*float*): 
            * ``kmesh`` ([*int*, *int*, *int*]): 
            * ``nepts`` (*int*): 

Output nodes:
    * ``last_doscal_dosdata`` (*XyData*): 
    * ``last_doscal_dosdata_interpol`` (*XyData*): 
    * ``last_doscal_results`` (*ParameterData*): 
    * ``last_params_voronoi`` (*ParameterData*): 
    * ``last_voronoi_remote`` (*RemoteData*): 
    * ``last_voronoi_results`` (*ParameterData*): 
    * ``results_vorostart_wc`` (*ParameterData*): 


                          
Example Usage
-------------

First load KKRcode and Voronoi codes::

    from aiida.orm import Code
    kkrcode = Code.get_from_string('KKRcode@my_mac')
    vorocode = Code.get_from_string('voronoi@my_mac')
    
Then choose some settings for the KKR specific parameters (LMAX cutoff etc.)::

    from aiida_kkr.tools.kkr_params import kkrparams
    kkr_settings = kkrparams(NSPIN=1, LMAX=2)
    
Now we create a structure node for the system we want to calculate::

    # create Copper bulk aiida Structure
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')
    
Finally we run the ``kkr_startpot_wc`` workflow (here using the defaults for the workflow settings)::

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

Case 1: Start from previous calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    

Case 2: Start from structure and run voronoi calculation first 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Workflow: ``aiida_kkr.workflows.check_para_convergence``

.. warning:: Not implemented yet!

Idea is to run checks after convergence for the following parameters:
    * RMAX
    * GMAX
    * cluster radius
    * energy contour
    * kmesh 
   

Find magnetic ground state
++++++++++++++++++++++++++

Workflow: ``aiida_kkr.workflows.check_magnetic_state``

.. warning:: Not implemented yet!

The idea is to run a Jij calculation to estimate if the ferromagnetic state is 
the ground state or not. Then the unit cell could be doubled to compute the 
antiferromagnetic state. In case of noncollinear magnetism the full Jij tensor 
should be analyzed.


Prepare KKRimp startup
++++++++++++++++++++++
 
.. warning:: Not implemented yet!

Possible inputs:
    * ``impurity_info``: i.e. scoef cluster info and Zimp etc.
    * ``reuse_host_GF``, optional: old kkrflex files that are reused in the same cluster for a different impurity (different Zimp)
    * ``parent_kkr_calc``: needed if host GF has to be computed, also needed to extract host systems infos
    * ``voro_code``: needed to get starting potential from ``impurity_info`` (i.e. Zimp)
    * ``kkr_code``: needed if host GF has to be computed
    
Possible output nodes:
    * ``output_parameters`` output node with dict containing infos of workflow run
    * ``potential_imp`` single file data needed for impurity calculation

Idea of the workflow:
    #. check if host GF is reused and verify consistency of cluster info
    #. run host GF generation if necessary
    #. run voronoi step for impurity starting potential (needs to create auxiliary structure first)
    #. use ``neworder_potential_wf`` workfunction to create single file data node


KKRimp DOS
++++++++++

.. warning:: Not implemented yet!

Idea:
    * create host GF for DOS contour if not given as input
    * run KKRimp step with dos settings


KKRimp self-consistency
+++++++++++++++++++++++

.. warning:: Not implemented yet!

Idea: 
    #. Use KKRimp startup workflow first if no parent KKRimp calculation is given
    #. run several steps of KKRimp until convergence
    #. check impurity DOS
    
    
    