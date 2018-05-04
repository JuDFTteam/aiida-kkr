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
