#!/usr/bin/env python

# connect to aiida db
from aiida import load_profile
load_profile()
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
def wait_for_it(calc, maxwait=300):
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
### !!! adapt to your code name !!! ###
codename = 'voronoi@my_mac'
code = Code.get_from_string(codename)

# create new instance of a VoronoiCalculation
voro_calc = code.new_calc()

# and set resources that will be used (here serial job)
voro_calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})

### !!! use queue name if necessary !!! ###
# voro_calc.set_queue_name('<quene_name>')

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
# KKR step (20 iterations simple mixing)
###################################################

# create new set of parameters for a KKR calculation and fill with values from previous voronoin calculation
params = kkrparams(params_type='kkr', **voro_params.get_dict())

# and set the missing values
params.set_multiple_values(RMAX=7., GMAX=65.)

# choose 20 simple mixing iterations first to preconverge potential (here 5% simple mixing)
params.set_multiple_values(NSTEPS=20, IMIX=0, STRMIX=0.05)

# create aiida ParameterData node from the KKR parameters
ParaNode = ParameterData(dict=params.get_dict())

# get KKR code and create new calculation instance
### !!! use your code name !!! ###
code = Code.get_from_string('KKRcode@my_mac')
kkr_calc = code.new_calc()

# set input Parameter, parent calulation (previous voronoi calculation), computer resources
kkr_calc.use_parameters(ParaNode)
kkr_calc.use_parent_folder(voronoi_calc_folder)
kkr_calc.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})

### !!! use queue name if necessary !!! ###
# kkr_calc.set_queue_name('<quene_name>')

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
kkr_calc_parent_folder = kkr_calc.out.remote_folder # parent remote folder of previous calculation
kkr_calc_continued.use_parent_folder(kkr_calc_parent_folder)
kkr_calc_continued.set_resources({'num_machines': 1, 'num_mpiprocs_per_machine':1})

### !!! use queue name if necessary !!! ###
# kkr_calc_continued.set_queue_name('<quene_name>')

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

### !!! use queue name if necessary !!! ###
# GF_host_calc.set_queue_name('<quene_name>')

# prepare impurity_info node containing the information about the impurity cluster
imp_info = ParameterData(dict={'Rcut':1.01, 'ilayer_center':0, 'Zimp':[79.]})
# set impurity info node to calculation
GF_host_calc.use_impurity_info(imp_info)

# store input nodes and submit calculation
GF_host_calc.store_all()
GF_host_calc.submit()

# wait for calculation to finish
wait_for_it(GF_host_calc)


######################################################################
# KKRimp calculation (20 simple mixing iterations  for preconvergence)
######################################################################

# first create impurity start pot using auxiliary voronoi calculation

# creation of the auxiliary styructure:
# use an aiida workfunction to keep track of the provenance
from aiida.work import workfunction as wf
@wf
def change_struc_imp_aux_wf(struc, imp_info): # Note: works for single imp at center only!
    from aiida.common.constants import elements as PeriodicTableElements
    _atomic_numbers = {data['symbol']: num for num, data in PeriodicTableElements.iteritems()}

    new_struc = StructureData(cell=struc.cell)
    isite = 0
    for site in struc.sites:
        sname = site.kind_name
        kind = struc.get_kind(sname)
        pos = site.position
        zatom = _atomic_numbers[kind.get_symbols_string()]
        if isite == imp_info.get_dict().get('ilayer_center'):
            zatom = imp_info.get_dict().get('Zimp')[0]
        symbol = PeriodicTableElements.get(zatom).get('symbol')
        new_struc.append_atom(position=pos, symbols=symbol)
        isite += 1

    return new_struc

new_struc = change_struc_imp_aux_wf(voro_calc.inp.structure, imp_info)

# then Voronoi calculation for auxiliary structure
### !!! use your code name !!! ###
codename = 'voronoi@my_mac'
code = Code.get_from_string(codename)
voro_calc_aux = code.new_calc()
voro_calc_aux.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})
voro_calc_aux.use_structure(new_struc)
voro_calc_aux.use_parameters(kkrcalc_converged.inp.parameters)
voro_calc_aux.store_all()
voro_calc_aux.submit()
### !!! use queue name if necessary !!! ###
# voro_calc_aux.set_queue_name('<quene_name>')

# wait for calculation to finish
wait_for_it(voro_calc_aux)

# then create impurity startpot using auxiliary voronoi calc and converged host potential

from aiida_kkr.tools.common_workfunctions import neworder_potential_wf

potname_converged = kkrcalc_converged._POTENTIAL
potname_imp = 'potential_imp'
neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.out.retrieved.get_abs_path('scoef'), skiprows=1)[:,3]-1]
potname_impvorostart = voro_calc_aux._OUT_POTENTIAL_voronoi
replacelist_pot2 = [[0,0]]

settings_dict = {'pot1': potname_converged,  'out_pot': potname_imp, 'neworder': neworder_pot1,
                 'pot2': potname_impvorostart, 'replace_newpos': replacelist_pot2, 'label': 'startpot_KKRimp',
                 'description': 'starting potential for Au impurity in bulk Cu'}
settings = ParameterData(dict=settings_dict)

startpot_Au_imp_sfd = neworder_potential_wf(settings_node=settings,
                                            parent_calc_folder=kkrcalc_converged.out.remote_folder,
                                            parent_calc_folder2=voro_calc_aux.out.remote_folder)

# now create KKRimp calculation and run first (some simple mixing steps) calculation

# needed to link to host GF writeout calculation
GF_host_output_folder = GF_host_calc.out.remote_folder

# create new KKRimp calculation
from aiida_kkr.calculations.kkrimp import KkrimpCalculation
kkrimp_calc = KkrimpCalculation()

### !!! use your code name !!! ###
kkrimp_code = Code.get_from_string('KKRimp@my_mac')

kkrimp_calc.use_code(kkrimp_code)
kkrimp_calc.use_host_Greenfunction_folder(GF_host_output_folder)
kkrimp_calc.use_impurity_potential(startpot_Au_imp_sfd)
kkrimp_calc.set_resources(resources)
kkrimp_calc.set_computer(kkrimp_code.get_computer())

# first set 20 simple mixing steps
kkrimp_params = kkrparams(params_type='kkrimp')
kkrimp_params.set_multiple_values(SCFSTEPS=20, IMIX=0, MIXFAC=0.05)
ParamsKKRimp = ParameterData(dict=kkrimp_params.get_dict())
kkrimp_calc.use_parameters(ParamsKKRimp)

# store and submit
kkrimp_calc.store_all()
kkrimp_calc.submit()

# wait for calculation to finish
wait_for_it(kkrimp_calc)


###################################################
# continued KKRimp calculation until convergence
###################################################

kkrimp_calc_converge = kkrimp_code.new_calc()
kkrimp_calc_converge.use_parent_calc_folder(kkrimp_calc.out.remote_folder)
kkrimp_calc_converge.set_resources(resources)
kkrimp_calc_converge.use_host_Greenfunction_folder(kkrimp_calc.inp.GFhost_folder)

kkrimp_params = kkrparams(params_type='kkrimp', **kkrimp_calc.inp.parameters.get_dict())
kkrimp_params.set_multiple_values(SCFSTEPS=99, IMIX=5, MIXFAC=0.05)
ParamsKKRimp = ParameterData(dict=kkrimp_params.get_dict())
kkrimp_calc_converge.use_parameters(ParamsKKRimp)

### !!! use queue name if necessary !!! ###
# kkrimp_calc_converge.set_queue_name('<quene_name>')

# store and submit
kkrimp_calc_converge.store_all()
kkrimp_calc_converge.submit()

wait_for_it(kkrimp_calc_converge)
