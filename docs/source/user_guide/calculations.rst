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
    # load essential aiida classes
    from aiida.orm import Code
    from aiida.orm import DataFactory
    StructureData = DataFactory('structure')
    ParameterData = DataFactory('parameter')
    
    # load kkrparms class which is a useful tool to create the set of input parameters for KKR-family of calculations
    from aiida_kkr.tools.kkr_params import kkrparams
    
    # load some python modules
    from numpy import array
    import os
    
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



Voronoi starting potential generator
++++++++++++++++++++++++++++++++++++

``kkr.voro``


Create an aiida structure
-------------------------

Create aiida StructureData node for bulk Cu::

    alat = 3.61 # lattice constant in Angstroem
    bravais = alat*array([[0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]) # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0,0,0], symbols='Cu')

    
Create KKR parameter set
------------------------

create empty set of KKR parameters (LMAX cutoff etc. ) for voronoi code::

    params = kkrparams(params_type='voronoi')
    
we can find out which parameters are mandatory to be set::

    missing_params = params.get_missing_keys(use_aiida=True)
    print('missing mandatory parameters:')
    print(missing_params)
    
and set at least the mandatory parameters::

    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    
finally create an aiida ParameterData node and fill with the dictionary of parameters::

    ParaNode = ParameterData(dict=params.get_dict())

    
Get Voronoi code
----------------

choose a valid installation of the voronoi code::

    codename = 'voronoi@my_mac'
    code = Code.get_from_string(codename)

    
Create calculation
------------------

create new instance of a VoronoiCalculation::

    voro_calc = code.new_calc()

and set resources that will be used (here serial job)::

    voro_calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})

then set structure and input parameter::

    voro_calc.use_structure(Cu)
    voro_calc.use_parameters(ParaNode)


Submit calculation and retrieve results
---------------------------------------

store all nodes and submit the calculation::

    voro_calc.store_all()
    voro_calc.submit()
    # wait for calc to finish
    wait_for_it(voro_calc)

check calculation state (or use verdi calculation list -a -p1)::

    print voro_calc.get_state()

get parsed output node::

    parsed_output_node = voro_calc.get_outputs_dict()

    

KKR calculation for bulk and interfaces
+++++++++++++++++++++++++++++++++++++++

``kkr.kkr``

Reuse settings from voronoi calculation::

    voronoi_calc_folder = voro_calc.out.remote_folder
    voro_params = voro_calc.inp.parameters

    
Update KKR parameter set
------------------------

create new set of parameters for a KKR calculation and fill with values from previous voronoin calculation::

    params = kkrparams(params_type='kkr', **voro_params.get_dict())
    
    # find out which parameters are missing for the KKR calculation
    print params.get_missing_keys()
    # and set the missing values
    params.set_multiple_values(RMAX=7., GMAX=65.)
    
    # create aiida ParameterData node from the KKR parameters
    ParaNode = ParameterData(dict=params.get_dict())


Get code and create calculation
------------------------------

get KKR code and create new calculation instance::

    code = Code.get_from_string('KKRcode@my_mac')
    kkr_calc = code.new_calc()
    
    # set input Parameter, parent calulation (previous voronoi calculation), computer resources 
    kkr_calc.use_parameters(ParaNode)
    kkr_calc.use_parent_folder(voronoi_calc_folder)
    kkr_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})

    
Submit calculation and retrieve results
---------------------------------------

store nodes and submit calculation::

    kkr_calc.store_all()
    kkr_calc.submit()
    
    wait_for_it(kkr_calc)
    

Continue KKR calculation frmo KKR parent
----------------------------------------

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


    
KKR calculation importer
++++++++++++++++++++++++

``kkr.kkrimporter``


Get remote code and set remote working directory
------------------------------------------------


Set file names
--------------


Submit KKR importer calculation and retrieve results
----------------------------------------------------





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





