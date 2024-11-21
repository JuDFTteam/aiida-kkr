# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools doing STM-like scans around impurity clusters
"""

import numpy as np
from aiida import orm, engine
from aiida_kkr.tools import find_parent_structure
from aiida_kkr.tools.combine_imps import get_scoef_single_imp
from aiida_kkr.tools.imp_cluster_tools import pos_exists_already, combine_clusters
from masci_tools.io.common_functions import get_alat_from_bravais

__copyright__ = (u'Copyright (c), 2023, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.2'
__contributors__ = (u'Philipp Rüßmann', u'Raffaele Aliberti')

##############################################################################
# combine impurty clusters


def convert_to_imp_cls(host_structure, imp_info):
    """
    convert imp info to rcls form
    """
    if 'imp_cls' in imp_info.get_dict():
        clust1 = imp_info['imp_cls']
        imp_info_cls = imp_info
    else:
        # convert Zimp, Rcut info to imp_cls info
        clust1 = get_scoef_single_imp(host_structure, imp_info)
        imp_info_cls = orm.Dict({'imp_cls': clust1, 'Zimp': imp_info['Zimp'], 'Rimp_rel': [[0., 0., 0.]]})
    return imp_info_cls, clust1


def get_Zadd(host_structure, add_position):
    """
    get Zatom for adding position
    """
    from aiida.common.constants import elements as PeriodicTableElements

    _atomic_numbers = {data['symbol']: num for num, data in PeriodicTableElements.items()}

    kind_name = host_structure.sites[add_position['ilayer']].kind_name
    symbol = host_structure.get_kind(kind_name).symbol
    Zadd = float(_atomic_numbers[symbol])
    return Zadd


def get_imp_cls_add(host_structure, add_position):
    """
    define auxiliary imp_info for adding position and generate rcls
    """
    Zadd = get_Zadd(host_structure, add_position)
    ilayer = add_position['ilayer']
    imp_info2 = orm.Dict({'ilayer_center': ilayer, 'Zimp': [Zadd], 'Rcut': 1e-5})
    # old version is too slow:
    # clust2 = get_scoef_single_imp(host_structure, imp_info2)
    # new version creates the array without calling the get_scoef_single_imp function:
    clust2 = np.array([[0., 0., 0., ilayer + 1, 0., 0.]])
    return imp_info2, clust2


def get_r_offset(clust1, clust2, host_structure, add_position):
    """
    find offset vector in rcls units
    """
    # calculate out-of plane vector from the ilayer indices of the two clusters
    r_out_of_plane = np.array([0., 0., 0.])
    # minus 1 because of conversion from fortran to python standard (counting starts at 0)
    ilayer1 = int(clust1[0, 3]) - 1
    ilayer2 = int(clust2[0, 3]) - 1
    if ilayer1 != ilayer2:
        pos1 = np.array(host_structure.sites[ilayer1].position)
        pos2 = np.array(host_structure.sites[ilayer2].position)
        r_out_of_plane = pos2 - pos1

    # convert from Ang to alat units (used internally in KKR)
    alat = get_alat_from_bravais(np.array(host_structure.cell), host_structure.pbc[2])
    r_out_of_plane /= alat

    # remove spurious offsets that might be there due to the choice of the unit cell positions
    r_out_of_plane = np.round(r_out_of_plane, 7)
    r_out_of_plane[:2] %= 1  # modulo 1 for x and y coordinate

    # calculate in-plane vector from da, db inputs
    da = add_position.get_dict().get('da', 0)
    db = add_position.get_dict().get('db', 0)
    cell = np.array(host_structure.cell)
    r_in_plane = da * cell[0] + db * cell[1]

    # convert from Ang to alat units (used internally in KKR)
    r_in_plane /= alat

    # combine to offset vector
    r_offset = r_out_of_plane + r_in_plane

    return r_offset


def offset_clust2(clust1, clust2, host_structure, add_position):
    """
    Compute and add offset to clust2
    """
    r_offset = get_r_offset(clust1, clust2, host_structure, add_position)

    clust2_offset = clust2.copy()
    clust2_offset[:, :3] += r_offset

    return clust2_offset


def get_imp_info_add_position(add_position, host_structure, imp_info):
    """
    Create combined impurity info node for the original
    imp cluster + an additional (STM tip) position
    """

    # extract host structure
    # host_structure = find_parent_structure(host_calc)

    # convert imp info to cls form
    imp_info_cls, clust1 = convert_to_imp_cls(host_structure, imp_info)

    # get imp cluster for adding position
    imp_info2, clust2 = get_imp_cls_add(host_structure, add_position)

    # shift clust2 by offset
    clust2_offset = offset_clust2(clust1, clust2, host_structure, add_position)

    # combine cluster information
    pos_exists_in_imp1, _ = pos_exists_already(clust1, clust2)
    if pos_exists_in_imp1:
        raise ValueError('Additional position exists already in impurity cluster.')
    cluster_combined, rimp_rel_combined, _, _ = combine_clusters(clust1, clust2_offset, False, debug=False)
    # combine the zimp arrays
    zimp_combined = imp_info['Zimp'] + imp_info2['Zimp']
    # now combine the imp info node
    imp_info_combined = orm.Dict({'imp_cls': cluster_combined, 'Zimp': zimp_combined, 'Rimp_rel': rimp_rel_combined})

    return imp_info_combined


@engine.calcfunction
def get_imp_info_add_position_cf(add_position, host_structure, imp_info):
    """
    Create a new impurty info node that combines the impurity cluster
    of an original calculation and an STM scanning position.
    """

    # then combine the imp info
    imp_info_combined = get_imp_info_add_position(add_position, host_structure, imp_info)

    return imp_info_combined


##############################################################################
# combine potentials


def extract_host_potential(add_position, host_calc):
    """
    Extract the potential of the position in the host that matches the additional position
    """

    # find ilayer from input node
    ilayer = add_position['ilayer']

    # get host calculation from remote
    #host_calc = host_remote.get_incoming(node_class=orm.CalcJobNode).first().node

    # read potential from host's retrieved node
    with host_calc.outputs.retrieved.open('out_potential') as _f:
        pot_txt = _f.readlines()
    iline_startpot = np.array([i for i, l in enumerate(pot_txt) if 'exc:' in l])

    # extract nspin from host calc
    nspin = host_calc.inputs.parameters['NSPIN']

    # get host's potential from ilayer
    pot_add = []
    for ispin in range(nspin):
        istart = iline_startpot[ilayer * nspin + ispin]
        iend = iline_startpot[ilayer * nspin + ispin + 1]
        pot_add += pot_txt[istart:iend]

    return pot_add


def add_host_potential_to_imp(add_position, host_calc, imp_potential_node):
    """
    combine host potential with impurity potential
    """
    # get add potential from host
    potname = f'host_pot:{add_position["ilayer"]}'
    if potname in imp_potential_node.extras:
        # reuse existing host position if we have found it previously
        pot_add = imp_potential_node.extras[potname]
    else:
        # get host postition and store as extra
        pot_add = extract_host_potential(add_position, host_calc)
        imp_potential_node.set_extra(potname, pot_add)

    # get impurity potential and convert to list
    pot_imp = imp_potential_node.get_content().split('\n')
    pot_imp = [line + '\n' for line in pot_imp if line != '']

    # glue potentials together and create SinglefileData
    pot_combined = pot_imp + pot_add

    return pot_combined


def create_combined_potential_node(add_position, host_calc, imp_potential_node):
    """
    Combine impurity potential with an additional potential from the host for
    the STM tip position (additional position)
    """
    import io

    # combine potential texts
    pot_combined = add_host_potential_to_imp(add_position, host_calc, imp_potential_node)

    # convert to byte string and put into SinglefilData node
    pot_combined_string = ''
    for line in pot_combined:
        pot_combined_string += line
    pot_combined_node = orm.SinglefileData(io.BytesIO(bytes(pot_combined_string, 'utf8')))

    return pot_combined_node


@engine.calcfunction
def create_combined_potential_node_cf(add_position, host_calc, imp_potential_node):
    """
    Calcfunction that combines the impurity potential with an addition potential site from the host
    """

    pot_combined_node = create_combined_potential_node(add_position, host_calc, imp_potential_node)

    return pot_combined_node


#####################################################################
# STM pathfinder


def STM_pathfinder(host_remote):
    """This function is used to help visualize the scanned positions
       and the symmetries that are present in the system

    inputs::
    host_remote : RemoteData : The Remote data contains all the information needed to create the path to scan

    outputs::
    struc_info : Dict  : Dictionary containing the structural information of the film
    matrices   : Array : Array containing the matrices that generate the symmetries of the system
    """

    def info_creation(structure):
        from ase.spacegroup import get_spacegroup
        # List of the Bravais vectors
        vec_list = structure.cell.tolist()

        # Find the Bravais vectors that are in plane vectors
        plane_vectors = {'plane_vectors': [], 'space_group': ''}
        for vec in vec_list:
            # Is this sufficient to find all the in-plane vectors?
            if vec[2] == 0:
                plane_vectors['plane_vectors'].append(vec[:2])

        space_symmetry = get_spacegroup(structure)
        plane_vectors['space_group'] = space_symmetry.no

        return plane_vectors

    def symmetry_finder(struc_info):
        from ase.spacegroup import Spacegroup
        # Here we get the symmetry operations that are possible
        symmetry_matrices = Spacegroup(struc_info['space_group'])

        # Reduce the dimensionality, we only want the 2D matrices
        matrices = []
        for element in symmetry_matrices.get_rotations():
            matrices.append(element[:2, :2])

        # Uniqueness of the elements must be preserved
        unique_matrices = []
        for matrix in matrices:
            if not any(np.array_equal(matrix, m) for m in unique_matrices):
                unique_matrices.append(matrix)

        return unique_matrices

    struc = find_parent_structure(host_remote)
    # clone the structure since it has already been saved in AiiDA and cannot be modified
    supp_struc = struc.clone()

    # If the structure is not periodic in every direction we force it to be.
    if not supp_struc.pbc[2]:
        # find film thickness
        zs = np.array([i.position[2] for i in supp_struc.sites])
        z = zs.max() - zs.min() + 5  # add 5 to have a unit cell larger than the considered film thickness
        # set third bravais vector along z direction
        cell = supp_struc.cell
        cell[2] = [0, 0, z]
        supp_struc.set_cell(cell)
        # change periodic boundary conditions to periodic
        supp_struc.pbc = (True, True, True)

    # ASE struc
    ase_struc = supp_struc.get_ase()

    # Structural informations are stored here
    struc_info = info_creation(ase_struc)

    # The structural informations are then used to find the symmetries of the system
    symm_matrices = symmetry_finder(struc_info)

    return struc_info, symm_matrices


@engine.calcfunction
def STM_pathfinder_cf(host_structure):
    """
    Calcfunction that gives back the structural information of the film, and the symmetries of the system
    """

    struc_info, symm_matrices = STM_pathfinder(host_structure)

    return struc_info, symm_matrices


###################################################################
# lattice generation (function of lattice plot)


def lattice_generation(x_len, y_len, rot, vec):
    import math
    """
    x_len : int  : value to create points between - x and x.
    y_len : int  : value to create points between - y and y.
    rot   : list : list of the rotation matrices given by the symmetry of the system.
    vec   : list : list containing the two Bravais vectors.
    """

    # Here we create a grid in made of points whic are the linear combination of the lattice vectors
    lattice_points = []

    for i in range(-x_len, x_len + 1):
        lattice_points_col = []
        for j in range(-y_len, y_len + 1):
            p = [i * x + j * y for x, y in zip(vec[0], vec[1])]
            lattice_points_col.append(p)
        lattice_points.append(lattice_points_col)

    # Eliminiatio of the symmetrical sites
    points_to_eliminate = []

    for i in range(-x_len, x_len + 1):
        for j in range(-y_len, y_len + 1):
            if ((lattice_points[i][j][0] > 0 or math.isclose(lattice_points[i][j][0], 0, abs_tol=1e-3)) and
                (lattice_points[i][j][1] > 0 or math.isclose(lattice_points[i][j][1], 0, abs_tol=1e-3))):
                for element in rot[1:]:
                    point = np.dot(element, lattice_points[i][j])
                    if point[0] >= 0 and point[1] >= 0:
                        continue
                    else:
                        points_to_eliminate.append(point)

    points_to_scan = []

    for i in range(-x_len, x_len + 1):
        for j in range(-y_len, y_len + 1):
            eliminate = False
            for elem in points_to_eliminate:
                # Since there can be some numerical error in the dot product we use the isclose function
                if (
                    math.isclose(elem[0], lattice_points[i][j][0], abs_tol=1e-4) and
                    math.isclose(elem[1], lattice_points[i][j][1], abs_tol=1e-4)
                ):
                    eliminate = True
            if not eliminate:
                points_to_scan.append(lattice_points[i][j])

    return points_to_eliminate, points_to_scan


###############################################################################################
# lattice plot


def lattice_plot(plane_vectors, symm_vec, symm_matrices, grid_length_x, grid_length_y):
    #from aiida_kkr.tools.tools_STM_scan import lattice_generation
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    origin = np.array([[0, 0], [0, 0]])
    # Generation of the points to plot
    unused, used = lattice_generation(grid_length_x, grid_length_y, symm_matrices, plane_vectors)

    # Plotting of the points
    for element in unused:
        plt.scatter(element[0], element[1], marker='s', s=130, c='#33638DFF')

    for element in used:
        plt.scatter(element[0], element[1], marker='D', s=130, c='#FDE725FF')

    # Plot of the crystal symmetry directions, tag must be activated.
    if symm_vec:
        import numpy.linalg as lin

        for element in symm_matrices:
            eig_val, eig_vec = lin.eig(element)

        for element in eig_vec:
            plt.quiver(
                *origin, element[0], element[1], alpha=1, color='#B8DE29FF', angles='xy', scale_units='xy', scale=1
            )

    # Plot of the Bravais lattice
    for element in plane_vectors:
        plt.quiver(*origin, element[0], element[1], color='#3CBB75FF', angles='xy', scale_units='xy', scale=1)

    legend_elements = [
        Line2D([0], [0], color='#33638DFF', lw=2, label='Unscanned Sites', marker='s'),
        Line2D([0], [0], color='#FDE725FF', lw=2, label='Scanned Sites', marker='D'),
        Line2D([0], [0], color='#3CBB75FF', lw=2, label='Bravais lattice'),
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(0.75, -0.15))

    plt.title('Lattice plot and symmetry directions')
    plt.ylabel('y direction')
    plt.xlabel('x direction')
    #plt.xticks(np.arange(-grid_length, grid_length, float(plane_vectors[0][0])))
    #plt.set_cmap(cmap)
    plt.grid(linestyle='--')
    plt.show()


##########################################################
# find linear combination coefficients


def find_linear_combination_coefficients(plane_vectors, vectors):
    from operator import itemgetter
    """This helper function takes the planar vectors and a list of vectors
       and return the coefficients in the base of the planar vectors"""

    # Formulate the system of equations Ax = b
    A = np.vstack((plane_vectors[0], plane_vectors[1])).T

    # Solve the system of equations using least squares method
    data = []
    for element in vectors:
        b = element
        # We use the least square mean error procedure to evaulate the units of da and db
        # lstsq returns: coeff, residue, rank, and singular value
        # We only need the coefficient.
        data.append(np.linalg.lstsq(A, b, rcond=None))

    indices = []
    for element in data:
        supp = []
        for elem in element[0]:
            # Here we round to an integer, this is because of the numerical error
            # which is present inside the calculation.
            supp.append(round(elem))
        indices.append(supp)

    # Before returning the indices, we reorder them first from the lowest to the highest valued
    # on the x axis and then from the lowest to the highest on the y axis.

    indices = sorted(indices, key=itemgetter(0, 1))

    return indices
