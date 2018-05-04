===============
Calculations
===============

Some intro...

Prerequisite: aiida knowledge, KKR knowledge, installation of KKRcode, Voronoi, aiida-kkr, code, computer setup
output of ``verdi calculation plugins`` should contain::

    $ verdi calculation plugins
    * kkr.kkr
    * kkr.kkrimp
    * kkr.kkrimporter
    * kkr.voro
    
In the following first set up calculation creating shapefunctions and starting potential in voronoi calculation.
Example: bulk Cu

Initial setup::

    # connect to aiida db
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()


Voronoi starting potential generator
++++++++++++++++++++++++++++++++++++

``kkr.voro``


First we create an aiida structure. To connect to the aiida-db we need::
    
    # get aiida StructureData class:
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')

Then we cal create the aiida StructureData node (here for bulk Cu)::
    
    alat = 3.61 # lattice constant in Angstroem
    bravais = [[0.5*alat, 0.5*alat, 0], [0.5*alat, 0, 0.5*alat], [0, 0.5*alat, 0.5*alat]] # Bravais matrix in Ang. units
    # now create StructureData instance and set Bravais matrix and atom in unit cell
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')

    
Next we create an empty set of KKR parameters (LMAX cutoff etc. ) for voronoi code::

    # load kkrparms class which is a useful tool to create the set of input parameters for KKR-family of calculations
    from aiida_kkr.tools.kkr_params import kkrparams
    params = kkrparams(params_type='voronoi')
    
.. note:: we can find out which parameters are mandatory to be set using 
          ``missing_params = params.get_missing_keys(use_aiida=True)``
and set at least the mandatory parameters::

    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    
finally create an aiida ParameterData node and fill with the dictionary of parameters::

    ParameterData = DataFactory('parameter') # use DataFactory to get ParamerterData class
    ParaNode = ParameterData(dict=params.get_dict())

    
Now we get the voronoi code::

    from aiida.orm import Code # load aiida 'Code' class
    
    codename = 'voronoi@my_mac'
    code = Code.get_from_string(codename)

and create new instance of a VoronoiCalculation::

    voro_calc = code.new_calc()

and set resources that will be used (here serial job)::

    voro_calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})

then set structure and input parameter::

    voro_calc.use_structure(Cu)
    voro_calc.use_parameters(ParaNode)


Now we are ready to submit the calculation. For that we first need to store the input nodes nodes and then submit the calculation::

    voro_calc.store_all()
    voro_calc.submit()

.. note:: check calculation state (or use verdi calculation list -a -p1) using 
          ``voro_calc.get_state()``

After the calculation is done we can get the parsed output node::

    parsed_output_node = voro_calc.get_outputs_dict()

which  should look like this::

    some output 
    

KKR calculation for bulk and interfaces
+++++++++++++++++++++++++++++++++++++++

Plugin: ``kkr.kkr``

Three input nodes:
    * ``parameters`` KKR parameter fitting the requirements for a KKR calculation (ParameterData)
    * ``parent_folder`` parent calulation remote folder node (RemoteFolder)
    * ``code`` KKR code node (code)
    
.. note:: The parent calculation can be one of the following:
             #. Voronoi calculation, initial calculation starting from structure
             #. previous KKR calculation, e.g. preconverged calculation
          The necessary structure information is always extracted from the voronoi parent calculation. 
          In case of a continued calculation the voronoi parent is recuresively searched for.
          
Special features exist where a fourth input node is persent and which triggers special behavior of the KKR calculation:
    * ``impurity_info`` Node specifying the impurity cluster (ParameterData)
    
    
Start KKR calculation from voronoi parent
-----------------------------------------

Reuse settings from voronoi calculation::

    voronoi_calc_folder = voro_calc.out.remote_folder
    voro_params = voro_calc.inp.parameters

    
Now we update the KKR parameter set to meet the requirements for a KKR calculation (slightly different than voronoi calculation).
Thus, we create a new set of parameters for a KKR calculation and fill the already set values from the previous voronoin calculation::

    # new kkrparams instance for KKR calculation
    params = kkrparams(params_type='kkr', **voro_params.get_dict())
    
    # set the missing values
    params.set_multiple_values(RMAX=7., GMAX=65.)
    
    # create aiida ParameterData node from the KKR parameters
    ParaNode = ParameterData(dict=params.get_dict())

.. note:: You can find out which parameters are missing for the KKR calculation using ``params.get_missing_keys()``

Now we can get the KKR code and create a new calculation instance and set the input nodes accordingly::

    code = Code.get_from_string('KKRcode@my_mac')
    kkr_calc = code.new_calc()
    
    # set input Parameter, parent calulation (previous voronoi calculation), computer resources 
    kkr_calc.use_parameters(ParaNode)
    kkr_calc.use_parent_folder(voronoi_calc_folder)
    kkr_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})

We can then run the KKR calculation by again storing the input nodes and submit the calculation::

    kkr_calc.store_all()
    kkr_calc.submit()
        

Continue KKR calculation from KKR parent calculation
----------------------------------------------------

create new KKR calculation instance to continue KKR ontop of a previous KKR calclation::

    kkr_calc_continued = code.new_calc()

reuse old KKR parameters and update scf settings (default is NSTEPS=1, IMIX=0)::

    params.set_multiple_values(NSTEPS=50, IMIX=5)

and create aiida ParameterData node::

    ParaNode = ParameterData(dict=params.get_dict())

then set input nodes for calculation::

    kkr_calc_continued.use_code(code)
    kkr_calc_continued.use_parameters(ParaNode)
    kkr_calc_continued.use_parent_folder(voronoi_calc_folder)
    kkr_calc_continued.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})

store input nodes and submit calculation::

    kkr_calc_continued.store_all()
    kkr_calc_continued.submit()
    
The finished calculation should have this output: ...

    
Special run modes (additional input nodes)
------------------------------------------

Writeout host Green function (for KKRimp calculation)
Create impurity_info node, 

Here we take the remote folder of the converged calculation to reuse settings and write out Green function and tmat of the crystalline host system::

    kkr_converged_parent_folder = kkr_calc_continued.out.remote_folder

Now we cal extract the parameters of the kkr calculation and add the ``KKRFLEX`` run-option::

    kkrcalc_converged = kkr_converged_parent_folder.get_inputs()[0]
    kkr_params_dict = kkrcalc_converged.inp.parameters.get_dict()
    kkr_params_dict['RUNOPT'] = ['KKRFLEX']
    
The parameters dictionary is not passed to the aiida ParameterData node::

    ParaNode = ParameterData(dict=kkr_params_dict)
    
Now we create a new KKR calculation and set input nodes::

    code = kkrcalc_converged.get_code() # take the same code as in the calculation before
    GF_host_calc= code.new_calc()
    resources = kkrcalc_converged.get_resources()
    GF_host_calc.set_resources(resources)
    GF_host_calc.use_parameters(ParaNode)
    GF_host_calc.use_parent_folder(kkr_converged_parent_folder)
    # prepare impurity_info node containing the information about the impurity cluster
    imp_info = ParameterData(dict={'Rcut':1.01, 'imp_pos':[0], 'Zimp':[29.]}) # choose host-in-host calculation first
    # set impurity info node to calculation
    GF_host_calc.use_impurity_info(imp_info)
    
The the calculation can be submitted::

    # store input nodes and submit calculation
    GF_host_calc.store_all()
    GF_host_calc.submit()

The output should look like this ...



KKR impurity calculation
++++++++++++++++++++++++

``kkr.kkrimp``

Create impurity cluster
-----------------------

Writeout host Green function with KKRcode
-----------------------------------------

Create impurity starting potential
----------------------------------

Set input parameter for KKRimp
------------------------------

Submit impurity calculation and retrieve results
------------------------------------------------



    
KKR calculation importer
++++++++++++++++++++++++

``kkr.kkrimporter``


Get remote code and set remote working directory
------------------------------------------------


Set file names
--------------


Submit KKR importer calculation and retrieve results
----------------------------------------------------




Example scripts
+++++++++++++++

Here is a small collection of 

Full example Voronoi-KKR-KKRimp
-------------------------------

Compact script starting with structure setup, then voronoi calculation, followed by 
initial KKR claculation which is then continued for convergence. The converged calculation 
is then used to write out the host GF and a simple inmpurity calculation is performed.

Download: :download:`this example script <../examples/kkr_short_example.py>`

::

    #!/usr/bin/env python
    
    # connect to aiida db
    from aiida import load_dbenv, is_dbenv_loaded
    if not is_dbenv_loaded():
        load_dbenv()
    # load essential aiida classes
    from aiida.orm import Code
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    ParameterData = DataFactory('parameter')
    
    # load kkrparms class which is a useful tool to create the set of input parameters for KKR-family of calculations
    from aiida_kkr.tools.kkr_params import kkrparams
    
    # load some python modules
    from numpy import array
    
    # helper function
    def wait_for_it(calc, maxwait=200):
        from time import sleep
        N = 0
        print 'start waiting for calculation to finish'
        while not calc.has_finished() and N<(maxwait/2.):
            N += 1
            if N%5==0:
                print('.')
            sleep(2.)
        print('waiting done after {} seconds: {} {}'.format(N*2, calc.has_finished(), calc.has_finished_ok()))
    
    
    ###################################################
    # initial structure
    ###################################################
    
    # create Copper bulk aiida Structure
    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')
    
    
    ###################################################
    # Voronoi step (preparation of starting potential)
    ###################################################
    
    # create empty set of KKR parameters (LMAX cutoff etc. ) for voronoi code
    params = kkrparams(params_type='voronoi')
    
    # and set at least the mandatory parameters
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    
    # finally create an aiida ParameterData node and fill with the dictionary of parameters
    ParaNode = ParameterData(dict=params.get_dict())
    
    # choose a valid installation of the voronoi code
    codename = 'voronoi@my_mac'
    code = Code.get_from_string(codename)
    
    # create new instance of a VoronoiCalculation
    voro_calc = code.new_calc()
    
    # and set resources that will be used (here serial job)
    voro_calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})
    
    # then set structure and input parameter
    voro_calc.use_structure(Cu)
    voro_calc.use_parameters(ParaNode)
    
    # store all nodes and submit the calculation
    voro_calc.store_all()
    voro_calc.submit()
    
    wait_for_it(voro_calc)
    
    # for future reference
    voronoi_calc_folder = voro_calc.out.remote_folder
    voro_params = voro_calc.inp.parameters
    
    
    ###################################################
    # KKR step (single iteration)
    ###################################################
    
    # create new set of parameters for a KKR calculation and fill with values from previous voronoin calculation
    params = kkrparams(params_type='kkr', **voro_params.get_dict())
    
    # and set the missing values
    params.set_multiple_values(RMAX=7., GMAX=65.)
    
    # create aiida ParameterData node from the KKR parameters
    ParaNode = ParameterData(dict=params.get_dict())
    
    # get KKR code and create new calculation instance
    code = Code.get_from_string('KKRcode@my_mac')
    kkr_calc = code.new_calc()
    
    # set input Parameter, parent calulation (previous voronoi calculation), computer resources
    kkr_calc.use_parameters(ParaNode)
    kkr_calc.use_parent_folder(voronoi_calc_folder)
    kkr_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})
    
    # store nodes and submit calculation
    kkr_calc.store_all()
    kkr_calc.submit()
    
    # wait for calculation to finish
    wait_for_it(kkr_calc)
    
    
    ###################################################
    # 2nd KKR step (continued from previous KKR calc)
    ###################################################
    
    # create new KKR calculation instance to continue KKR ontop of a previous KKR calclation
    kkr_calc_continued = code.new_calc()
    
    # reuse old KKR parameters and update scf settings (default is NSTEPS=1, IMIX=0)
    params.set_multiple_values(NSTEPS=50, IMIX=5)
    # and create aiida ParameterData node
    ParaNode = ParameterData(dict=params.get_dict())
    
    # then set input nodes for calculation
    kkr_calc_continued.use_code(code)
    kkr_calc_continued.use_parameters(ParaNode)
    kkr_calc_continued.use_parent_folder(voronoi_calc_folder)
    kkr_calc_continued.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})
    
    # store input nodes and submit calculation
    kkr_calc_continued.store_all()
    kkr_calc_continued.submit()
    
    # wait for calculation to finish
    wait_for_it(kkr_calc_continued)
    
    
    ###################################################
    # writeout host GF (using converged calculation)
    ###################################################
    
    # take remote folder of converged calculation to reuse setting and write out Green function and tmat of the crystalline host system
    kkr_converged_parent_folder = kkr_calc_continued.out.remote_folder
    
    # extreact kkr calculation from parent calculation folder
    kkrcalc_converged = kkr_converged_parent_folder.get_inputs()[0]
    
    # extract parameters from parent calculation and update RUNOPT for KKRFLEX option
    kkr_params_dict = kkrcalc_converged.inp.parameters.get_dict()
    kkr_params_dict['RUNOPT'] = ['KKRFLEX']
    
    # create aiida ParameterData node with set parameters that are updated compared to converged parent kkr calculation
    ParaNode = ParameterData(dict=kkr_params_dict)
    
    # create new KKR calculation
    code = kkrcalc_converged.get_code() # take the same code as in the calculation before
    GF_host_calc= code.new_calc()
    
    # set resources, Parameter Node and parent calculation
    resources = kkrcalc_converged.get_resources()
    GF_host_calc.set_resources(resources)
    GF_host_calc.use_parameters(ParaNode)
    GF_host_calc.use_parent_folder(kkr_converged_parent_folder)
    
    # prepare impurity_info node containing the information about the impurity cluster
    imp_info = ParameterData(dict={'Rcut':1.01, 'imp_pos':[0], 'Zimp':[29.]}) # choose host-in-host calculation first
    # set impurity info node to calculation
    GF_host_calc.use_impurity_info(imp_info)
    
    # store input nodes and submit calculation
    GF_host_calc.store_all()
    GF_host_calc.submit()
    
    # wait for calculation to finish
    wait_for_it(GF_host_calc)
    
    
    ###################################################
    # KKRimp calculation (single iteration first)
    ###################################################
    
    SingleFileData = DataFactory('singlefile')
    pot_imp_path = '<path-to-impurity-potential>'
    potfile_imp = SingleFileData()
    potfile_imp.set_file(pot_imp_path)
    
    # needed to link to host GF writeout calculation
    GF_host_output_folder = GF_host_calc.out.remote_folder
    
    # create new KKRimp calculation
    from aiida_kkr.calculations.kkrimp import KkrimpCalculation
    kkrimp_calc = KkrimpCalculation()
    
    kkrimp_code = Code.get_from_string('KKRimp@my_mac')
    
    kkrimp_calc.use_code(kkrimp_code)
    kkrimp_calc.use_host_Greenfunction_folder(GF_host_output_folder)
    kkrimp_calc.use_impurity_potential(potfile_imp)
    kkrimp_calc.set_resources(resources)
    kkrimp_calc.set_computer(kkrimp_code.get_computer())
    
    #test = kkrimp_calc.submit_test()
    #folder = test[0]
    
    # store and submit
    kkrimp_calc.store_all()
    kkrimp_calc.submit()
    
    # wait for calculation to finish
    wait_for_it(kkrimp_calc)
    
    
    ###################################################
    # continued KKRimp calculation
    ###################################################
    
    kkrimp_parent_calc_folder = kkrimp_calc.out.remote_folder
    
    kkrimp_code = kkrimp_calc.get_code()
    kkrimp_calc_continued = KkrimpCalculation()
    kkrimp_calc_continued.use_code(kkrimp_code)
    kkrimp_calc_continued.use_host_Greenfunction_folder(GF_host_output_folder)
    kkrimp_calc_continued.use_parent_calc_folder(kkrimp_parent_calc_folder)
    kkrimp_calc_continued.set_resources(resources)
    kkrimp_calc_continued.set_computer(code.get_computer())
    kkrimp_calc_continued.use_parameters(ParameterData(dict=kkrparams(params_type='kkrimp', IMIX=5, SCFSTEPS=50).get_dict()))
    
    kkrimp_calc_continued.store_all()
    kkrimp_calc_continued.submit()
    
    wait_for_it(kkrimp_calc_continued)