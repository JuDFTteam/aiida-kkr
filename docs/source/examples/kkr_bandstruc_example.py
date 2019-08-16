#!/usr/bin/env python

# connect to aiida db
from aiida import load_profile
load_profile()
# load essential aiida classes
from aiida.orm import Code, DataFactory, load_node
StructureData = DataFactory('structure')
ParameterData = DataFactory('parameter')

# helper function:
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


# some settings (parent calculations):
    
# converged KKR calculation (taken form bulk Cu KKR example)
kkr_calc_converged = load_node(24951)
# previous DOS calculation started from converged KKR calc (taken from KKRimp DOS example, i.e. GF host calculation with DOS contour)
host_dos_calc = load_node(25030)


# generate kpoints for bandstructure calculation

from aiida_kkr.calculations.voro import VoronoiCalculation
struc, voro_parent = VoronoiCalculation.find_parent_structure(kkr_calc_converged.out.remote_folder)

from aiida.tools.data.array.kpoints import get_explicit_kpoints_path
kpts = get_explicit_kpoints_path(struc).get('explicit_kpoints')


# run bandstructure calculation

# create bandstructure calculation reusing old settings (including same computer and resources in this example)
kkrcode = kkr_calc_converged.get_code()
kkrcalc = kkrcode.new_calc()
kkrcalc.use_kpoints(kpts) # pass kpoints as input
kkrcalc.use_parent_folder(kkr_calc_converged.out.remote_folder)
kkrcalc.set_resources(kkr_calc_converged.get_resources())
# change parameters to qdos settings (E range and number of points)
from aiida_kkr.tools.kkr_params import kkrparams
qdos_params = kkrparams(**kkr_calc_converged.inp.parameters.get_dict()) # reuse old settings
# reuse the same emin/emax settings as in DOS run (extracted from input parameter node)
qdos_params.set_multiple_values(EMIN=host_dos_calc.inp.parameters.get_dict().get('EMIN'), 
                                EMAX=host_dos_calc.inp.parameters.get_dict().get('EMAX'), 
                                NPT2=100)
kkrcalc.use_parameters(ParameterData(dict=qdos_params.get_dict()))

# store and submit calculation
kkrcalc.store_all()
kkrcalc.submit()

wait_for_it(kkrcalc, maxwait=600)


# plot results

# extract kpoint labels
klbl = kpts.labels
# fix overlapping labels (nicer plotting)
tmp = klbl[2]
tmp = (tmp[0], '\n'+tmp[1]+' ')
klbl[2] = tmp
tmp = klbl[3]
tmp = (tmp[0], '  '+tmp[1])
klbl[3] = tmp

#plotting of bandstructure and previously calculated DOS data

# load DOS data
from masci_tools.io.common_functions import interpolate_dos
dospath_host = host_dos_calc.out.retrieved.get_abs_path('')
ef, dos, dos_interpol = interpolate_dos(dospath_host, return_original=True)
dos, dos_interpol = dos[0], dos_interpol[0]

# load qdos file and reshape
from numpy import loadtxt, sum, log
qdos_file = kkrcalc.out.retrieved.get_abs_path('qdos.01.1.dat')
q = loadtxt(qdos_file)
nepts = len(set(q[:,0]))
data = q[:,5:].reshape(nepts, len(q)/nepts, -1)
e = (q[::len(q)/nepts, 0]-ef)*13.6

# plot bandstructure
from matplotlib.pyplot import figure, pcolormesh, show, xticks, ylabel, axhline, axvline, gca, title, plot, ylim, xlabel, suptitle
figure(figsize=((8, 4.8)))
pcolormesh(range(len(q)/nepts), e, log(sum(abs(data), axis=2)), lw=0)
xticks([i[0] for i in klbl], [i[1] for i in klbl])
ylabel('E-E_F (eV)')
axhline(0, color='lightgrey', lw=1)
title('band structure')

# plot DOS on right hand side of bandstructure plot
axBand = gca()
from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(axBand)
axDOS = divider.append_axes("right", 1.2, pad=0.1, sharey=axBand)

plot(dos_interpol[:,1]/13.6, (dos_interpol[:,0]-ef)*13.6)

ylim(e.min(), e.max())

axhline(0, color='grey', lw=1)
axvline(0, color='grey', lw=1)

axDOS.yaxis.set_tick_params(labelleft=False, labelright=True, right=True, left=False)
xlabel('states/eV')

title('DOS')
suptitle(struc.get_formula(), fontsize=16)

show()
