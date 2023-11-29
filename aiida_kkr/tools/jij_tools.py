"""
Tools for Jij calculations and parsing
"""

# parsing jij output

import numpy as np
from masci_tools.io.common_functions import search_string, get_aBohr2Ang, get_Ry2eV
from aiida.orm import ArrayData, StructureData, Dict, load_node, Bool
from aiida.engine import calcfunction, submit
from aiida_kkr.calculations import KkrCalculation
from aiida_kkr.tools import find_parent_structure


def get_sites(structure):
    """Get all sites also for a CPA structure"""
    sites = []  # for CPA
    for site in structure.sites:
        sitekind = structure.get_kind(site.kind_name)
        for ikind in range(len(sitekind.symbols)):
            sites.append(site)
    return sites


def get_jijs_shells(jij_calc, verbose=False):
    """read the jij.atom files from the retrieved and determine if the file has DMI or not"""

    # read the jij.atom files
    jijs_shells_all = []
    for jijfile in [i for i in jij_calc.outputs.retrieved.list_object_names() if 'Jij' in i]:
        i_index = int(jijfile.split('.atom')[1])
        if verbose:
            print('load jij.atom', i_index)
        with jij_calc.outputs.retrieved.open(jijfile) as f:
            try:
                tmp = np.loadtxt(f)
            except:
                # this means we could not read the file which happens, for example,
                # in case of CPA where "&" characters are inserted to mark the end of block
                f.seek(0)
                txt = f.readlines()
                txt = [line for line in txt if '&' not in line]
                tmp = np.loadtxt(txt)
        # sort by radius
        jij_atom = tmp[tmp[:, 0].argsort()]
        if len(jij_atom[0]) > 4:
            # add i index, needed for mapping to spirit data
            # only needed for DMI mode
            jij_atom[:, -1] = i_index
        if verbose:
            print('shape jij.atom', jij_atom.shape)
        jijs_shells_all += jij_atom.tolist()

    if verbose:
        print('jijs_shells.shape', len(jijs_shells_all))

    if len(jijs_shells_all) == 0:
        raise ValueError('Did not extract any Jij shells')

    # format of jijs_shells differs in old and new solver:
    if len(jijs_shells_all[0]) == 4:
        # - isotropic exchange constants only (old solver). The columns refer to:
        #   * [0] |Rij| in units of the lattice constant
        #   * [1] Jij in Ryd
        #   * [2] shell-number
        #   * [3] atom type of atom j
        #   * [4] atom type of atom i
        dmimode = False
    else:
        # - full exchange tensor (new solver). The columns refer to:
        #   * [0] ∣Rij∣ in units of the lattice constant
        #   * [1] Jij (isotropic part) in Ryd
        #   * [2] Dij (anti-symmetric DMI part) in Ryd
        #   * [3] Sij (diagonal traceless part) in Ryd
        #   * [4] Aij (off-diagonal symmetric part) in Ryd
        #   * [5-7] Rj− Ri (3 component vector)
        #   * [8] atom type of atom j
        #   * [9] atom type of atom i
        dmimode = True

    return jijs_shells_all, dmimode


def expand_jijs_iso(jij_calc, jijs_shells_all, shells, alat):
    """expand jijs from shells to all pairs which is easier to use in spirit
    This is the roune that deals with isotropic interactions (i.e. using the old solver)"""

    # extract basis vectors ot the structure, needed to convert from cartesian coordinates (x,y,z)
    # to relative coordinates (da, db, dc) such that (x, y, z)^T = da*a + db*b + dc*c
    # with a,b,c being the three unit vectors
    structure = find_parent_structure(jij_calc)
    cell = np.array(structure.cell)
    cell_positions = np.array([i.position for i in structure.sites])

    # transformation matrix from absolute to realtive coordinates
    U2rel = np.linalg.inv(cell.transpose())

    # expand shells data to all structure that spirit can understand (i.e. take all pairs)
    jijs_expanded, positions_expanded = [], []
    for jshell in jijs_shells_all:
        for s in shells[int(jshell[3] - 1)]:
            # i and j indices in the unit cell
            iatom, jatom = int(s[0] - 1), int(s[1] - 1)

            # absolute distance
            # Rx, Ry, Rz = alat*(s[7]-s[4]), alat*(s[8]-s[5]), alat*(s[9]-s[6])
            Rx, Ry, Rz = alat * s[7], alat * s[8], alat * s[9]

            # substract sub-lattice distance to find correct unit cell scalings
            R = np.array([Rx, Ry, Rz])
            R -= cell_positions[jatom]
            # R -= (cell_positions[jatom] - cell_positions[i])
            x, y, z = R

            # calculate da, db, dc ( i.e. the mutiplicities of the unit vectors)
            da, db, dc = np.array(np.round(np.dot(U2rel, [x, y, z])), dtype=int)

            #                      i,      j,   da, db, dc, Jij (meV)
            jijs_expanded.append([iatom, jatom, da, db, dc, jshell[1] * get_Ry2eV() * 1000])

            positions_expanded.append(np.array([Rx, Ry, Rz]) - cell_positions[iatom])

    # convert to numpy array
    jijs_expanded = np.array(jijs_expanded)
    positions_expanded = np.array(positions_expanded)

    return jijs_expanded, positions_expanded


def expand_jijs_dmi(jijs_shells_x, jijs_shells_y, jijs_shells_z, shells, cell, sites, alat, verbose=False):
    """Bring output Jijs of the new solver into the right form
    (KKR output already is expanded in shells but it needs to be brought into the right form
    by combining the x,y,z calculations and the Rvec information)
    This will give the full Jij tensor and prodiuce the columns that spirit can understand
    """

    # conversion factor from Ry to meV
    Ry2meV = 1000 * get_Ry2eV()

    # transformation matrix from absolute to realtive coordinates
    U2rel = np.linalg.inv(cell.transpose())

    jijs_expanded, positions_expanded = [], []
    for ipos, jz in enumerate(jijs_shells_z):
        # also load the information from m||x and m||y
        jx = jijs_shells_x[ipos]
        jy = jijs_shells_y[ipos]

        if verbose:
            print(ipos, jx, jy, jz)

        # get i and j indices of the atoms
        # -1 to convert from fortran output to python indexing (starting at 0)
        iatom = int(jz[-1]) - 1
        jatom = int(jz[-2]) - 1

        # Rvec = pos_i - pos_j + R_ij[j], with R_ij[j] given in the Jij.atom files (in alat units)
        dx = alat * jz[5]
        dy = alat * jz[6]
        dz = alat * jz[7]
        x = sites[iatom].position[0] - sites[jatom].position[0] + dx
        y = sites[iatom].position[1] - sites[jatom].position[1] + dy
        z = sites[iatom].position[2] - sites[jatom].position[2] + dz

        # calculate da, db, dc ( i.e. the mutiplicities of the unit vectors)
        da, db, dc = np.array(np.round(np.dot(U2rel, [x, y, z])), dtype=int)

        if verbose:
            if ipos == 0:
                print(cell)
            print(iatom, jatom, dx, dy, dz, x, y, z, da, db, dc)

        #collect full Jij tensor
        # 1: J, 2: D, 3: S, 4: A
        Jxx = (jz[1] + jz[3] + jy[1] + jy[3]) / 2. * Ry2meV
        Jxy = jz[4] * Ry2meV
        Jxz = jy[4] * Ry2meV
        Jyx = jz[4] * Ry2meV
        Jyy = (jz[1] - jz[3] + jx[1] + jx[3]) / 2. * Ry2meV
        Jyz = jx[4] * Ry2meV
        Jzx = jy[4] * Ry2meV
        Jzy = jx[4] * Ry2meV
        Jzz = (jx[1] - jx[3] + jy[1] - jy[3]) / 2. * Ry2meV
        # calculate DMI vector
        Dx = jx[2] * Ry2meV
        # for y component no minus sign: ERROR in documentation!
        #Dy = -jy[2] * Ry2meV
        Dy = jy[2] * Ry2meV
        Dz = jz[2] * Ry2meV

        # collect data in big arrays
        #                       i,     j,   da, db, dc,   Jij (meV), Dij vector(meV), full Jij tensor (meV)
        jijs_expanded.append([
            iatom, jatom, da, db, dc, (Jxx + Jyy + Jzz) / 3., Dx, Dy, Dz, Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz
        ])
        positions_expanded.append([dx, dy, dz])

    # convert collected data to numpy array
    jijs_expanded = np.array(jijs_expanded)
    positions_expanded = np.array(positions_expanded)

    return jijs_expanded, positions_expanded


def get_jij_structure(structure, jijs_expanded, jij_calc):
    """make auxiliary structure that has only the ij sites from the Jij calculation"""
    struc_jij_sites = StructureData(cell=structure.cell)
    struc_jij_sites.pbc = structure.pbc
    all_sites_jij = set(list(jijs_expanded[:, 0]) + list(jijs_expanded[:, 1]))
    isite, icount, mappings = -1, 0, []  # for mapping to the sites of the reduced structure
    mappings_back = {}
    for site in structure.sites:
        sitekind = structure.get_kind(site.kind_name)
        for ikind, symbol in enumerate(sitekind.symbols):
            isite += 1
            # take only structues for which Jij couplings are extracted
            if isite in all_sites_jij:
                mappings.append([icount, isite])
                mappings_back[isite] = icount
                icount += 1
                if ikind == 0:
                    struc_jij_sites.append_atom(
                        position=site.position, symbols=sitekind.symbols, weights=sitekind.weights
                    )

    # add spin moments extra (used by spirit plugin)
    calc_conv = jij_calc.inputs.parent_folder.get_incoming(node_class=KkrCalculation).first().node
    mu_s_all = calc_conv.outputs.output_parameters['magnetism_group']['spin_moment_per_atom']
    mu_s = []
    for ij in mappings:
        mu_s.append(mu_s_all[ij[1]])
    struc_jij_sites.set_extra('spin_moments', mu_s)

    return struc_jij_sites, mappings_back, mu_s


def get_shells(jij_calc, verbose=False):
    # read the shells information that is needed to map to the complete list of pairs
    with jij_calc.outputs.retrieved.open('shells.dat') as f:
        txt = f.readlines()
    nshell = int(txt.pop(0).split()[0])
    naez = len(find_parent_structure(jij_calc).sites)
    ioffset = 0
    shells = []
    for ishell in range(nshell):
        nat = int(txt[ioffset].split()[1])
        shell = []
        for iline in range(nat):
            shell.append(txt[ioffset + 1 + iline].split())
        shells.append(np.array(shell, dtype=float))
        ioffset += 1 + nat
    if verbose:
        print('found shells:', shells)
    return shells


@calcfunction
def parse_jij_calc_data(
    jij_calc_retrieved, jij_calc_x_retrieved=None, jij_calc_y_retrieved=None, verbose=lambda: Bool(False)
):
    """
    Parse the output of a Jij calculation from the retreived folder of a KkrCalculation

    :params jij_calc_retrieved: retrieved folder output of a KkrCalculation which ran with the Jij inputs (m||z is assumed)
    :params jij_calc_x_retrieved: like jij_calc_retrieved but for m||x (only needed for new solver)
    :params jij_calc_y_retrieved: like jij_calc_retrieved but for m||y (only needed for new solver)
    :params verbose: True/False can be used to print debugging output

    :returns:
      {'jij_data': jij_data, 'structure_jij_sites': struc_jij_sites}
      where
      * jij_data is and ArrayData that contains the expanded Jij's (see 'array_descriptions' extra for more details)
      * struc_jij_sites is the reduced structure that contains onlt the atoms which have Jij couplings
        (comes from the input of the KkrCalculation). The mappings to the original structure (and their i,j indices) is given as an extra
    """
    verbose = verbose.value

    # extract kkrCalculation from retreived child
    jij_calc = jij_calc_retrieved.get_incoming(node_class=KkrCalculation).first().node
    if verbose:
        print('jij(z) calculation:', jij_calc.uuid)

    # extract basis vectors ot the structure, needed to convert from cartesian coordinates (x,y,z)
    # to relative coordinates (da, db, dc) such that (x, y, z)^T = da*a + db*b + dc*c
    # with a,b,c being the three unit vectors
    structure = find_parent_structure(jij_calc)
    cell = np.array(structure.cell)
    natyp = len(get_sites(structure))
    if verbose:
        print('found structure with natyp:', natyp)

    # in KKR everything is scaled with the lattice constant
    alat = jij_calc.outputs.output_parameters['alat_internal'] * get_aBohr2Ang()

    # read jij.atom files
    jijs_shells_z, dmimode = get_jijs_shells(jij_calc, verbose)
    if verbose:
        print('shape jij_shells_z:', len(jijs_shells_z), dmimode)

    # read jij.atom files if a calculation in x and y are given in the input
    dmimode_x, dmimode_y = False, False
    if dmimode and jij_calc_x_retrieved is not None:
        jij_calc_x = jij_calc_x_retrieved.get_incoming(node_class=KkrCalculation).first().node
        if verbose:
            print('jij(x) calculation:', jij_calc_x.uuid)
        jijs_shells_x, dmimode_x = get_jijs_shells(jij_calc_x)
        if verbose:
            print('shape jij_shells_x:', len(jijs_shells_x), dmimode_x)
        if not dmimode_x:
            raise ValueError('jij_calc_x is not a DMI calculation (i.e. used old solver)')
    if dmimode and jij_calc_y_retrieved is not None:
        jij_calc_y = jij_calc_y_retrieved.get_incoming(node_class=KkrCalculation).first().node
        if verbose:
            print('jij(y) calculation:', jij_calc_y.uuid)
        jijs_shells_y, dmimode_y = get_jijs_shells(jij_calc_y)
        if verbose:
            print('shape jij_shells_y:', len(jijs_shells_y), dmimode_y)
        if not dmimode_y:
            raise ValueError('jij_calc_y is not a DMI calculation (i.e. used old solver)')
    # consistency check
    if dmimode and not (dmimode_x and dmimode_y):
        raise ValueError('Found dmimode but not all calculations (x,y and z) were given correctly')

    # read the shells information
    shells = get_shells(jij_calc, verbose)
    # print('shells', shells)

    # expand shells data to a structure that spirit can understand (i.e. take all pairs)
    if not dmimode:
        if verbose:
            print('expand isotropic Jijs')
        jijs_expanded, positions_expanded = expand_jijs_iso(jij_calc, jijs_shells_z, shells, alat)
    else:
        if verbose:
            print('expand anisotropic Jijs, Dij etc')
        jijs_expanded, positions_expanded = expand_jijs_dmi(
            jijs_shells_x, jijs_shells_y, jijs_shells_z, shells, cell, get_sites(structure), alat, verbose=verbose
        )

    # sort arrays
    isort = np.lexsort(jijs_expanded[:, :5][:, ::-1].transpose())
    jijs_expanded = jijs_expanded[isort]
    positions_expanded = positions_expanded[isort]

    # create an auxiliary structure that contains only the sites which are used in the Jij step
    # (i.e. we drop all sites where we don't have couplings)
    struc_jij_sites, mappings_back, mu_s = get_jij_structure(structure, jijs_expanded, jij_calc)

    if verbose:
        print(f'reduced structure has {len(struc_jij_sites.sites)} sites')

    # change i,j indices to match smaller structure
    # and correct sign if mu_i and mu_j are oriented antiparallel
    for ii, iatom in enumerate(jijs_expanded[:, 0]):
        jatom = jijs_expanded[ii, 1]
        jijs_expanded[ii, 0] = mappings_back[int(iatom)]
        jijs_expanded[ii, 1] = mappings_back[int(jatom)]
        # correct sign
        smu_i = np.sign(mu_s[int(jijs_expanded[ii, 0])])
        smu_j = np.sign(mu_s[int(jijs_expanded[ii, 1])])
        sign_factor = smu_i * smu_j
        jijs_expanded[ii, 5] *= sign_factor  # change sign of Jij
        jijs_expanded[ii, 6:9] *= sign_factor  # change sign of Dij

    # now collect the outputs in AiiDA Array objects
    jij_data = ArrayData()
    # jij_data.set_array('Jij_shells', jijs_shells_z)
    # if dmimode:
    #     # also add x and y shells output
    #     jij_data.set_array('Jij_shells_x', jijs_shells_x)
    #     jij_data.set_array('Jij_shells_y', jijs_shells_y)
    jij_data.set_array('Jij_expanded', jijs_expanded)
    jij_data.set_array('positions_expanded', positions_expanded)
    # add description to extras
    jij_data.extras['array_descriptions'] = {
        #         'Jij_shells': """Jij output in the shells that KKR found.
        #     The format differs for the type of calculation:
        #     - isotropic exchange constants only (old solver). The columns refer to:
        #       * [0] |Rij| in units of the lattice constant
        #       * [1] Jij in Ryd
        #       * [2] shell-number
        #       * [3] atom type of atom j
        #     - full exchange tensor (new solver). The columns refer to:
        #       * [0] ∣Rij∣ in units of the lattice constant
        #       * [1] Jij (isotropic part) in Ryd
        #       * [2] Dij (anti-symmetric DMI part) in Ryd
        #       * [3] Sij (diagonal traceless part) in Ryd
        #       * [4] Aij (off-diagonal symmetric part) in Ryd
        #       * [5-7] Rj - Ri (3 component vector)
        #       * [8] atom type of atom j

        #     if the full exchange tensor is calculated the Jij_shells_x and Jij_shells_y also exist""",
        'Jij_expanded':
        'i, j, da, db, dc, Jij (meV) [, Dij vector (x, y, z in meV), full Jij tensor (xx, xy, xz, yx, yy, yz, zx, zy, zz in meV)]',
        'positions_expanded': 'x, y, z (Ang.) positions of all pairs in Jij_expanded',
    }

    # add extras to generated structure for quick access
    struc_jij_sites.extras['mappings_ij'] = mappings_back
    struc_jij_sites.extras['uuid_struc_orig'] = structure.uuid
    struc_jij_sites.extras['uuid_jij_data'] = jij_data.uuid

    # return dict (link_label: node)
    return {'jij_data': jij_data, 'structure_jij_sites': struc_jij_sites}


def parse_jij_calc(jij_calc, jij_calc_x=None, jij_calc_y=None, verbose=False):
    """
    Parse outcome of Jij calculation and return jij_data and structure_jij_sites.
    Calculate not only Jij but also Dij vector if 3 calculations for m ||z, m||x and m||y are given.

    :param jij_calc: KkrCalculation with Jij inputs for m || z
    :param jij_calc_x: optional, KkrCalculation with Jij inputs for m || x
    :param jij_calc_y: optional, KkrCalculation with Jij inputs for m || y
    :param verbose: boolean that activate verbose writeout during parsing

    :return jij_data: ArrayData that contains the expanded Jij's (see 'array_descriptions' extra of the node for more details)
    :return struc_jij_sites: the reduced structure that contains onlt the atoms which have Jij couplings
        (comes from the input of the KkrCalculation). The mappings to the original structure (and their i,j indices) is given as an extra.
    """

    if not jij_calc.is_finished_ok:
        raise ValueError('Jij (z) calculation not finished ok (yet)!')
    if jij_calc_x is not None and not jij_calc_x.is_finished_ok:
        raise ValueError('Jij (x) calculation not finished ok (yet)!')
    if jij_calc_y is not None and not jij_calc_y.is_finished_ok:
        raise ValueError('Jij (y) calculation not finished ok (yet)!')

    if 'parsed_jij_data' in jij_calc.extras:
        # load from extras
        if verbose:
            print('found parsed jij data, load it now')
        jij_data_uuid = jij_calc.extras['parsed_jij_data']['jij_data_uuid']
        struc_jij_sites_uuid = jij_calc.extras['parsed_jij_data']['struc_jij_sites_uuid']
        jij_data = load_node(jij_data_uuid)
        structure_jij_sites = load_node(struc_jij_sites_uuid)
    else:
        # no old data found, do parsing here
        # the following calcfunction keeps the data provenance by taking the retrieved folders as inputs
        if jij_calc_x is None:
            jij_calc_x_retrieved = None
        else:
            jij_calc_x_retrieved = jij_calc_x.outputs.retrieved
        if jij_calc_y is None:
            jij_calc_y_retrieved = None
        else:
            jij_calc_y_retrieved = jij_calc_y.outputs.retrieved
        jij_parsed = parse_jij_calc_data(
            jij_calc.outputs.retrieved, jij_calc_x_retrieved, jij_calc_y_retrieved, verbose=Bool(verbose)
        )
        jij_data, structure_jij_sites = jij_parsed['jij_data'], jij_parsed['structure_jij_sites']

        # then save as extra to be able to load it in the next run
        if verbose:
            print('finished parsing jij data, save it as extra now')
        jij_calc.set_extra(
            'parsed_jij_data', {
                'jij_data_uuid': jij_data.uuid,
                'struc_jij_sites_uuid': structure_jij_sites.uuid
            }
        )

    # finally return result
    return jij_data, structure_jij_sites


def load_or_submit_cpa_jijcalcs(
    scf_cpa_wf,
    JIJSITEI=None,
    JIJSITEJ=None,
    tempr=400.,
    rclustz=0.9,
    kmesh=[100, 100, 100],
    NATOMIMPD=500,
    NSHELD=2000,
    JIJRAD=2.0,
    options=None,
    uuid_x=None,
    uuid_y=None,
    uuid_z=None
):
    """Load or submit Jij calculation for CPA"""

    scf_remote = scf_cpa_wf.outputs.last_RemoteData
    last_calc = scf_remote.get_incoming(node_class=KkrCalculation).first().node

    builder = last_calc.get_builder_restart()
    builder.parent_folder = scf_remote

    # set Jij parameters
    para_Jij = {k: v for k, v in scf_cpa_wf.outputs.last_InputParameters.get_dict().items() if v}
    para_Jij['TEMPR'] = tempr  # slightly reduce temperature
    para_Jij['RCLUSTZ'] = rclustz  # increase cluster radius
    para_Jij['BZDIVIDE'] = kmesh  # increase k-points
    para_Jij['NSTEPS'] = 1  # one-shot
    para_Jij['NATOMIMPD'] = NATOMIMPD  # array dimension
    para_Jij['NSHELD'] = NSHELD  # array dimension
    para_Jij['KPOIBZ'] = np.product(kmesh)  # array dimension
    # add 'XCPL' runopt to list of runopts (activates Jij calculation)
    runopts = para_Jij.get('RUNOPT', [])
    runopts.append('XCPL    ')
    para_Jij['RUNOPT'] = runopts
    # set Jij parameters
    # i and j index for Jij calculation in internal units
    # uses site index (i.e. needs to be <=10)
    if JIJSITEI is not None:
        para_Jij['JIJSITEI'] = JIJSITEI
        if JIJSITEJ is None:
            JIJSITEJ = JIJSITEI
    if JIJSITEJ is not None:
        para_Jij['JIJSITEJ'] = JIJSITEJ
    para_Jij['JIJRAD'] = JIJRAD  # radius in lattice constants up to which the Jijs are calculated

    builder.parameters = Dict(para_Jij)

    # starting angles for 3 directions, needed to extract full Jij tensor

    Nsites = len(get_sites(find_parent_structure(last_calc)))

    init_angles_x = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [90.0 for i in range(Nsites)],
        'phi': [0.0 for i in range(Nsites)],
    })

    init_angles_y = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [90.0 for i in range(Nsites)],
        'phi': [90.0 for i in range(Nsites)],
    })

    init_angles_z = Dict({
        'fix_dir': [True for i in range(Nsites)],
        'theta': [0.0 for i in range(Nsites)],
        'phi': [0.0 for i in range(Nsites)],
    })

    # submit m||z calculation
    builder.initial_noco_angles = init_angles_z
    builder.metadata.label = 'Jij_z'
    if options is not None:
        builder.metadata.options = options
    if uuid_z is not None:
        calc_Jij_z = load_node(uuid_z)
    else:
        calc_Jij_z = submit(builder)
        print('submitted z', calc_Jij_z)

    builder.initial_noco_angles = init_angles_y
    builder.metadata.label = 'Jij_y'
    if options is not None:
        builder.metadata.options = options
    if uuid_y is not None:
        calc_Jij_y = load_node(uuid_y)
    else:
        calc_Jij_y = submit(builder)
        print('submitted y', calc_Jij_y)

    builder.initial_noco_angles = init_angles_x
    builder.metadata.label = 'Jij_x'
    if options is not None:
        builder.metadata.options = options
    if uuid_x is not None:
        calc_Jij_x = load_node(uuid_x)
    else:
        calc_Jij_x = submit(builder)
        print('submitted x', calc_Jij_x)

    return calc_Jij_x, calc_Jij_y, calc_Jij_z
