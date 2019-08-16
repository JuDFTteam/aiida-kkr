#!/usr/bin/env python

# connect to aiida db
from aiida import load_profile
load_profile()
# load essential aiida classes
from aiida.orm import DataFactory, load_node
ParameterData = DataFactory('parameter')


# some settings:
#DOS contour (in Ry units), emax=EF+dE_emax:
emin, dE_emax, npt = -0.2, 0.1, 101
# kkrimp parent (converged imp pot, needs to tbe a KKRimp calculation node)
kkrimp_calc_converge = load_node(25025)

# derived quantities:
GF_host_calc = kkrimp_calc_converge.inp.GFhost_folder.inp.remote_folder
kkr_converged_parent_folder = GF_host_calc.inp.parent_calc_folder

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

################################################################################################

# first host GF with DOS contour
from aiida_kkr.tools.kkr_params import kkrparams
params = kkrparams(**GF_host_calc.inp.parameters.get_dict())
params.set_multiple_values(EMIN=emin, EMAX=GF_host_calc.res.fermi_energy+dE_emax, NPOL=0, NPT1=0, NPT2=npt, NPT3=0)
ParaNode = ParameterData(dict=params.get_dict())

code = GF_host_calc.get_code() # take the same code as in the calculation before
GF_host_doscalc= code.new_calc()
resources = GF_host_calc.get_resources()
GF_host_doscalc.set_resources(resources)
GF_host_doscalc.use_parameters(ParaNode)
GF_host_doscalc.use_parent_folder(kkr_converged_parent_folder)
GF_host_doscalc.use_impurity_info(GF_host_calc.inp.impurity_info)

# store and submit
GF_host_doscalc.store_all()
GF_host_doscalc.submit()

# wait for calculation to finish
print 'host GF calc for DOS contour'
wait_for_it(GF_host_doscalc)

# then KKRimp step using the converged potential

kkrimp_doscalc = kkrimp_calc_converge.get_code().new_calc()
kkrimp_doscalc.use_host_Greenfunction_folder(GF_host_doscalc.out.remote_folder)
kkrimp_doscalc.use_parent_calc_folder(kkrimp_calc_converge.out.remote_folder)
kkrimp_doscalc.set_resources(kkrimp_calc_converge.get_resources())

# set to DOS settings
params = kkrparams(params_type='kkrimp', **kkrimp_calc_converge.inp.parameters.get_dict())
params.set_multiple_values(RUNFLAG=['lmdos'], SCFSTEPS=1)
ParaNode = ParameterData(dict=params.get_dict())

kkrimp_doscalc.use_parameters(ParaNode)

# store and submit calculation
kkrimp_doscalc.store_all()
kkrimp_doscalc.submit()

# wait for calculation to finish

print 'KKRimp calc DOS'
wait_for_it(kkrimp_doscalc)

# Finally plot the DOS:

# get interpolated DOS from GF_host_doscalc calculation:
from masci_tools.io.common_functions import interpolate_dos
dospath_host = GF_host_doscalc.out.retrieved.get_abs_path('')
ef, dos, dos_interpol = interpolate_dos(dospath_host, return_original=True)
dos, dos_interpol = dos[0], dos_interpol[0]

# read in impurity DOS
from numpy import loadtxt
impdos0 = loadtxt(kkrimp_doscalc.out.retrieved.get_abs_path('out_lmdos.interpol.atom=01_spin1.dat'))
impdos1 = loadtxt(kkrimp_doscalc.out.retrieved.get_abs_path('out_lmdos.interpol.atom=13_spin1.dat'))
# sum over spins:
impdos0[:,1:] = impdos0[:,1:]*2
impdos1[:,1:] = impdos1[:,1:]*2

# plot bulk and impurity DOS
from matplotlib.pyplot import figure, fill_between, plot, legend, title, axhline, axvline, xlim, ylim, ylabel, xlabel, title, show
figure()
fill_between((dos_interpol[:,0]-ef)*13.6, dos_interpol[:,1]/13.6, color='lightgrey', lw=0, label='bulk Cu')
plot((impdos0[:,0]-ef)*13.6, impdos0[:,1]/13.6, label='Au imp')
plot((impdos0[:,0]-ef)*13.6, impdos1[:,1]/13.6, label='1st Cu neighbor')
plot((impdos0[:,0]-ef)*13.6, (impdos1[:,1]-dos_interpol[:,1])/dos_interpol[:,1], '--', label='relative difference in 1st Cu neighbor')
legend()
title('DOS of Au impurity embedded into bulk Cu')
axhline(0, lw=1, color='grey')
axvline(0, lw=1, color='grey')
xlim(-8, 1)
ylim(-0.5,8.5)
xlabel('E-E_F (eV)')
ylabel('DOS (states/eV)')
show()
