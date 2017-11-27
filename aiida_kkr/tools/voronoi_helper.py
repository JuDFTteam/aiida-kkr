#!/usr/bin/python

from __future__ import print_function
import sys

# redefine raw_input for python 3/2.7 compatilbility
if sys.version_info[0] >= 3:
    def raw_input(msg):
        return input(msg)


def get_corestates_from_potential(potfile='potential'):
    """Read core states from potential file"""
    from scipy import zeros
    txt = open(potfile).readlines()

    #get start of each potential part
    istarts = [iline for iline in range(len(txt)) if 'POTENTIAL' in txt[iline]]

    n_core_states = [] #number of core states per potential
    e_core_states = [] #energies of core states
    l_core_states = [] #angular momentum index, i.e. 0=s, 1=p etc...
    for ipot in range(len(istarts)):
        line = txt[istarts[ipot]+6]
        n = int(line.split()[0])
        n_core_states.append(n)
        elevels = zeros(n) #temp array for energies
        langmom = zeros(n, dtype=int) #temp array for angular momentum index
        for icore in range(n):
            line = txt[istarts[ipot]+7+icore].split()
            langmom[icore] = int(line[0])
            elevels[icore] = float(line[1].replace('D', 'E'))
        e_core_states.append(elevels)
        l_core_states.append(langmom)

    return n_core_states, e_core_states, l_core_states


def get_highest_core_state(nstates, energies, lmoments):
    """Find highest lying core state from list of core states, needed to find and check energy contour"""
    idx = energies.argmax()
    lval = lmoments[idx]
    nquant = sum(lmoments == lval) + lval
    level_descr = '%i%s'%(nquant, 'spdfgh'[lval])

    return lval, energies[idx], level_descr


def get_valence_min(outfile='out_voronoi'):
    """Construct minimum of energy contour (between valence band bottom and core states)"""
    from scipy import array
    txt = open(outfile).readlines()
    searchstr = 'All other states are above'
    valence_minimum = array([float(line.split(':')[1].split()[0]) for line in txt if searchstr in line])
    return valence_minimum


def check_voronoi_output(potfile, outfile):
    """Read output from voronoi code and create guess of energy contour"""
    from scipy import zeros
    #analyse core levels, minimum of valence band and their difference
    ncore, ecore, lcore = get_corestates_from_potential(potfile=potfile)
    e_val_min = get_valence_min(outfile=outfile)

    #print a table that summarizes the result
    e_core_max = zeros(len(ncore))
    print('pot    Highest core-level     low. val. state    diff')
    for ipot in range(len(ncore)):
        if ncore[ipot] > 0:
            lval, emax, descr = get_highest_core_state(ncore[ipot], ecore[ipot], lcore[ipot])
            e_core_max[ipot] = emax
            print('%3i     %2s %10.6f            %6.2f         %6.2f'%(ipot+1, descr, emax, e_val_min[ipot], e_val_min[ipot]-emax))
        else:
            print('%3i     << no core states >>'%(ipot+1))
            # set to some large negative number for check to not give false positive in case of empty cells
            e_core_max[ipot] = -100

    #get hint for energy integration:
    emin_guess = e_val_min.min()-0.2

    #make a check
    ediff_min = 1.5
    if (emin_guess-e_core_max.max()) < ediff_min:
        print_warning('Highest core level is too close to valence band: check results carefully!')

    return emin_guess

