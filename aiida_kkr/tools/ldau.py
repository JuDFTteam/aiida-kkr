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


def get_ldaumatrices(retrieved):
    """
    Take retrieved folder of KkrimpCalculation and extract ldaupot file
    If it exists we extract the LDAU matrices 'wldau', 'uldau' and 'phi'
    """
    from aiida.common.folders import SandboxFolder
    from aiida_kkr.calculations import KkrimpCalculation

    has_ldaupot_file = False
    txt_dict_ldaumats = {}

    # create Sandbox to extract ldaupot file there
    with SandboxFolder() as tempfolder:
        # extract ldaupot file to tempfolder if it exists
        has_ldaupot_file = KkrimpCalculation.get_ldaupot_from_retrieved(retrieved, tempfolder)
        # now read ldau matrices and store in txt_dict_ldaumats
        if has_ldaupot_file:
            with tempfolder.open(KkrimpCalculation._LDAUPOT + '_old') as ldaupot_file:
                # read file and look for keywords to identify where matrices are stored in the file
                txt = ldaupot_file.readlines()
                ii = 0
                for line in txt:
                    if 'wldau' in line:
                        iwldau = ii
                    if 'uldau' in line:
                        iuldau = ii
                    if 'phi' in line:
                        iphi = ii
                    ii += 1
                # save matrices to output dict
                txt_dict_ldaumats['wldau'] = txt[iwldau + 1:iuldau]  # pylint: disable=possibly-used-before-assignment
                txt_dict_ldaumats['uldau'] = txt[iuldau + 1:iphi]  # pylint: disable=possibly-used-before-assignment
                txt_dict_ldaumats['phi'] = txt[iphi + 1:]  # pylint: disable=possibly-used-before-assignment

    return has_ldaupot_file, txt_dict_ldaumats


def get_LDAU_initmatrices_dict(txts_ldaumat1, offset=0):
    """
    Put the initial matrices text content of `get_ldaumatrices` output in a
    Dictionnary form sorted by atoms.
    """

    # Loop through characters ('wldau', 'uldau', 'phi')
    for char in ['wldau', 'uldau', 'phi']:
        atom_line = []
        iatoms = []

        # Save atom number and lines in which atoms information starts
        for iline in range(len(txts_ldaumat1[char])):
            if txts_ldaumat1[char][iline][:4] == 'atom':
                iatom = txts_ldaumat1[char][iline].split()[1]
                iatoms.append(iatom)
                atom_line.append(iline)
        atom_line.append(len(txts_ldaumat1[char]))

        mat = []
        # Extract matrices based on atom lines
        for iatom in range(len(atom_line) - 1):
            mat.append(txts_ldaumat1[char][atom_line[iatom] + 1:atom_line[iatom + 1]])

        # Assign matrices to corresponding variables based on the character
        if char == 'wldau':
            wldaumat = mat
        elif char == 'uldau':
            uldaumat = mat
        elif char == 'phi':
            phimat = mat

    # Create a dictionary to store the organized data
    LDAU_initmatrices_dict = {}

    # Fill the dictionary with atom-wise information
    for iatom in range(len(iatoms)):
        LDAU_initmatrices_dict[f'iatom={int(iatoms[iatom])-1+offset}'] = {}
        LDAU_initmatrices_dict[f'iatom={int(iatoms[iatom])-1+offset}']['wldau'] = wldaumat[iatom]  # pylint: disable=possibly-used-before-assignment
        LDAU_initmatrices_dict[f'iatom={int(iatoms[iatom])-1+offset}']['uldau'] = uldaumat[iatom]  # pylint: disable=used-before-assignment
        LDAU_initmatrices_dict[f'iatom={int(iatoms[iatom])-1+offset}']['phi'] = phimat[iatom]  # pylint: disable=used-before-assignment

    return LDAU_initmatrices_dict
