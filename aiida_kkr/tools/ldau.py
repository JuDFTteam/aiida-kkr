# -*- coding: utf-8 -*-
"""
This module contains tools used in LDA+U calculations
"""


def get_ldaupot_text(ldau_settings, ef_Ry, natom, initialize=True, return_luj=False):
    """
    create the text for the ldaupot file
    """
    from masci_tools.io.common_functions import get_Ry2eV

    eV2Ry = 1. / get_Ry2eV()

    # these lists are extracted from ldau_settings node and then written to the ldaupot file
    iatoms_ldau = []
    lopt = []
    ueff = []
    jeff = []
    eref = []

    # extract values from ldau_settings
    for key, val in ldau_settings.items():
        if 'iatom' in key:
            iatoms_ldau.append(int(key.split('=')[1]))
            lopt.append(val['L'])
            # add values in Ry units
            jeff.append(val['J'] * eV2Ry)
            ueff.append(val['U'] * eV2Ry)
            eref.append(val['Eref_EF'] * eV2Ry + ef_Ry)

    if initialize:
        # this means we initialize this file
        ldaurun = 0
    else:
        # this means wldau etc are reused (need to be added to the file)
        ldaurun = 1

    # collect text which is written to ldaupot
    txt = [f'{ldaurun} ']
    txt_lopt, txt_jeff, txt_ueff, txt_eref = [], [], [], []
    ii = 0
    for iatom in range(natom):
        if iatom not in iatoms_ldau:
            txt_lopt += [f'{-1} ']
            txt_jeff += [f'{0.0} ']
            txt_ueff += [f'{0.0} ']
            txt_eref += [f'{0.0} ']
        else:
            txt_lopt += [f'{lopt[ii]} ']
            txt_jeff += [f'{jeff[ii]} ']
            txt_ueff += [f'{ueff[ii]} ']
            txt_eref += [f'{eref[ii]} ']
            ii += 1
    txt += ['\n'] + txt_lopt + ['\n'] + txt_ueff + ['\n'] + txt_jeff + ['\n'] + txt_eref
    txt += ['\nwldau\nuldau\nphi\n']

    # add initial matrices
    if initialize and 'initial_matrices' in ldau_settings.keys():
        # save for consistency check
        nldauatoms = ii

        # change first number from 0 to 1 to signal reading-in instead of calculating initial matrices
        txt[0] = '1 '

        # remove last dummy line, will be replaced with starting values now
        txt[-1] = '\n'
        txt_wldau = ['wldau\n']
        txt_uldau = ['uldau\n']
        txt_phi = ['phi\n']
        ii = 0
        for iatom in range(natom):
            if iatom in iatoms_ldau:
                txt_wldau += [f'atom {iatom+1}\n'] + ldau_settings['initial_matrices'][f'iatom={iatom}']['wldau']
                txt_uldau += [f'atom {iatom+1}\n'] + ldau_settings['initial_matrices'][f'iatom={iatom}']['uldau']
                txt_phi += [f'atom {iatom+1}\n'] + ldau_settings['initial_matrices'][f'iatom={iatom}']['phi']
                ii += 1  # count number of atoms

        # consistency check
        if nldauatoms != ii:
            raise ValueError('initial_matrices input inconsistent')

        # add additional lines to txt
        txt = txt + txt_wldau + txt_uldau + txt_phi

    if not return_luj:
        return txt
    else:
        return txt, lopt, jeff, ueff, eref, iatoms_ldau
