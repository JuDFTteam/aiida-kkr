# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools doing STM-like scans around impurity clusters
"""

import numpy as np
from aiida import orm, engine
from aiida_kkr.tools import find_parent_structure
from aiida.orm import CalcJobNode
from aiida_kkr.tools.imp_cluster_tools import get_scoef_single_imp
from aiida_kkr.tools.imp_cluster_tools import pos_exists_already, combine_clusters
from masci_tools.io.common_functions import get_alat_from_bravais

__copyright__ = (u'Copyright (c), 2023, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.6'
__contributors__ = (u'Philipp Rüßmann', u'Raffaele Aliberti')

##############################################################################
# combine impurty clusters


def convert_to_imp_cls(host_structure, imp_info):
    """
    convert imp info to rcls form
    """
    if 'imp_cls' in imp_info.get_dict():
        clust1 = np.array(imp_info['imp_cls'])
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
    #clust2 = get_scoef_single_imp(host_structure, imp_info2)
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

    # # remove spurious offsets that might be there due to the choice of the unit cell positions
    # r_out_of_plane = np.round(r_out_of_plane, 7)
    # r_out_of_plane[:2] %= 1  # modulo 1 for x and y coordinate

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
        # If the position exists already we simply skip the addition of the new scanning position
        return None
        #raise ValueError('Additional position exists already in impurity cluster.')
    cluster_combined, rimp_rel_combined, _, _ = combine_clusters(clust1, clust2_offset, False, debug=False)
    # combine the zimp arrays
    zimp1 = imp_info['Zimp']
    if not isinstance(zimp1, list):
        # convert to list if necessary
        zimp1 = [zimp1]
    zimp_combined = zimp1 + imp_info2['Zimp']
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


##############################################################################
# Helper function generating the point group symmetries in the system

def pointgrp(writesymfile=False):
    import numpy as np
    """
    Helper function contining the representation of the point group symmetry matrices as expressed
    in the KKR code. 
    """

    rotmat = np.zeros((64, 3, 3))
    rotname = [""] * 64

    rthree = np.sqrt(3.0) / 2.0
    half = 0.5

    rotmat[0, 0, 0] = 1.0
    rotmat[0, 1, 1] = 1.0
    rotmat[0, 2, 2] = 1.0
    rotname[0] = 'E'
    
    rotmat[1, 0, 1] = 1.0
    rotmat[1, 1, 2] = -1.0
    rotmat[1, 2, 0] = -1.0
    rotname[1] = 'C3alfa'

    rotmat[2, 0, 1] = -1.0
    rotmat[2, 1, 2] = -1.0
    rotmat[2, 2, 0] = 1.0
    rotname[2] = 'C3beta'

    rotmat[3, 0, 1] = -1.0
    rotmat[3, 1, 2] = 1.0
    rotmat[3, 2, 0] = -1.0
    rotname[3] = 'C3gamma'

    rotmat[4, 0, 1] = 1.0
    rotmat[4, 1, 2] = 1.0
    rotmat[4, 2, 0] = 1.0
    rotname[4] = 'C3delta'

    rotmat[5, 0, 2] = -1.0
    rotmat[5, 1, 0] = 1.0
    rotmat[5, 2, 1] = -1.0
    rotname[5] = 'C3alfa-1'

    rotmat[6, 0, 2] = 1.0
    rotmat[6, 1, 0] = -1.0
    rotmat[6, 2, 1] = -1.0
    rotname[6] = 'C3beta-1'

    rotmat[7, 0, 2] = -1.0
    rotmat[7, 1, 0] = -1.0
    rotmat[7, 2, 1] = 1.0
    rotname[7] = 'C3gamma-1'

    rotmat[8, 0, 2] = 1.0
    rotmat[8, 1, 0] = 1.0
    rotmat[8, 2, 1] = 1.0
    rotname[8] = 'C3delta-1'

    rotmat[9, 0, 0] = 1.0
    rotmat[9, 1, 1] = -1.0
    rotmat[9, 2, 2] = -1.0
    rotname[9] = 'C2x'

    rotmat[10, 0, 0] = -1.0
    rotmat[10, 1, 1] = 1.0
    rotmat[10, 2, 2] = -1.0
    rotname[10] = 'C2y'

    rotmat[11, 0, 0] = -1.0
    rotmat[11, 1, 1] = -1.0
    rotmat[11, 2, 2] = 1.0
    rotname[11] = 'C2z'

    rotmat[12, 0, 0] = 1.0
    rotmat[12, 1, 2] = 1.0
    rotmat[12, 2, 1] = -1.0
    rotname[12] = 'C4x'

    rotmat[13, 0, 2] = -1.0
    rotmat[13, 1, 1] = 1.0
    rotmat[13, 2, 0] = 1.0
    rotname[13] = 'C4y'

    rotmat[14, 0, 1] = 1.0
    rotmat[14, 1, 0] = -1.0
    rotmat[14, 2, 2] = 1.0
    rotname[14] = 'C4z'

    rotmat[15, 0, 0] = 1.0
    rotmat[15, 1, 2] = -1.0
    rotmat[15, 2, 1] = 1.0
    rotname[15] = 'C4x-1'

    rotmat[16, 0, 2] = 1.0
    rotmat[16, 1, 1] = 1.0
    rotmat[16, 2, 0] = -1.0
    rotname[16] = 'C4y-1'

    rotmat[17, 0, 1] = -1.0
    rotmat[17, 1, 0] = 1.0
    rotmat[17, 2, 2] = 1.0
    rotname[17] = 'C4z-1'

    rotmat[18, 0, 1] = 1.0
    rotmat[18, 1, 0] = 1.0
    rotmat[18, 2, 2] = -1.0
    rotname[18] = 'C2a'

    rotmat[19, 0, 1] = -1.0
    rotmat[19, 1, 0] = -1.0
    rotmat[19, 2, 2] = -1.0
    rotname[19] = 'C2b'

    rotmat[20, 0, 2] = 1.0
    rotmat[20, 1, 1] = -1.0
    rotmat[20, 2, 0] = 1.0
    rotname[20] = 'C2c'

    rotmat[21, 0, 2] = -1.0
    rotmat[21, 1, 1] = -1.0
    rotmat[21, 2, 0] = -1.0
    rotname[21] = 'C2d'

    rotmat[22, 0, 0] = -1.0
    rotmat[22, 1, 2] = 1.0
    rotmat[22, 2, 1] = 1.0
    rotname[22] = 'C2e'

    rotmat[23, 0, 0] = -1.0
    rotmat[23, 1, 2] = -1.0
    rotmat[23, 2, 1] = -1.0
    rotname[23] = 'C2f'

    for i1 in range(24):
        rotmat[i1+24] = -rotmat[i1]
        rotname[i1+24] = 'I' + rotname[i1]
    
    matrices = zip(rotname, rotmat)
    return list(matrices)


##############################################################################
# Parser of the symmetry contained in the host calculation

def symmetry_parser(host_calc):
    
    """
    Function used to get the relevant information regarding the symmetries of the sample directly from the claculation 
    node of the host calculation.
    
    Inputs :: 
    
    host_calc : CalcJobNode : Calculation node hosting the structural information regarding the sample
    
    Outputs :: 
    
    list_vec : list containing the normalized real space vectors constituting the geometry of the structure
    list_mat : list contining the rotatation matrices representing the symmetry of the system
    """
    
    sym_mat = pointgrp()
    
    with host_calc.outputs.retrieved.open('output.0.txt') as _f:
        read = _f.readlines()

    # Initialize variables to store symmetry information
    symmetry_operations = []
    
    # Flags to identify relevant sections
    in_symmetry_section = False
    start_read = False
    
    # Process the file line by line
    for line in read:
        line = line.strip()  # Remove leading and trailing whitespace
        if "3D symmetries" in line:
            # Start reading symmetry section
            in_symmetry_section = True
        if in_symmetry_section and line.startswith("------------------------------------------------------------"):
            start_read = True
        if in_symmetry_section and start_read and not line.startswith("------------------------------------------------------------"):
            # Found the line containing symmetry operations
            symmetry_operations.append(line.split())
        if in_symmetry_section and start_read and symmetry_operations and line.startswith("------------------------------------------------------------"):
            break # Exit the loop after every instance of the symmetries have been found
            
    sym_ops = [item for sublist in symmetry_operations for item in sublist]
            
    lattice_vectors = []
    
    in_lattice_section = False
    start_read = False
    #        
    for line in read:
        line = line.strip()  # Remove leading and trailing whitespace
        if "normalised (ALAT)" in line:
            # Start reading symmetry section
            in_lattice_section = True
        if in_lattice_section and line.startswith("----------------------                ----------------------"):
            start_read = True
        if in_lattice_section and start_read and not line.startswith("----------------------                ----------------------"):
            # Found the line containing symmetry operations
            lattice_vectors.append(line.split())
        if in_lattice_section and start_read and lattice_vectors and line.startswith("----------------------                ----------------------"):
            break # Exit the loop after every instance of the symmetries have been found
            
    # parsing of the data
    list_vec = []
    list_mat = []
            
    # Retrieve the vectors needed
    for l_vec in lattice_vectors:
        x = l_vec[1] ;y = l_vec[2]
        list_vec.append([float(x), float(y)])
    
    # Retrieve the matrix corresponding to the sysymmetry_operations
    for sym in sym_ops:
        for mat in sym_mat: 
            
            if sym == mat[0]:
                
                list_mat.append(mat[1][0:2, 0:2]) #only take the plane rotation
        
    # Output the extracted symmetry operations
    return list_vec, list_mat


##############################################################################
# STM pathfinder


#def STM_pathfinder(host_remote):
#    """
#    Calcfunction that gives back the structural information of the film, and the symmetries of the system
#
#    inputs ::
#
#           host_remote : Remote_data : node containing the remote data of the host material
#
#    return ::
#
#           plane_vectors   : list : list containing the 2 in plane vectors that span the surface.
#           unique_matrices : list : list of matrices, contains the 2x2 matrices that constitue the symmetry operations of the system.
#    """
#
#    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer, SymmOp
#
#    struc = find_parent_structure(host_remote)
#    # clone the structure since it has already been saved in AiiDA and cannot be modified
#    supp_struc = struc.clone()
#
#    # If the structure is not periodic in every direction we force it to be.
#    supp_struc.pbc = (True, True, True)
#
#    # Pymatgen struc
#    py_struc = supp_struc.get_pymatgen()
#
#    struc_dict = py_struc.as_dict()
#    # Find the Bravais vectors that are in-plane vectors (assumes 2D structure)
#    plane_vectors = {'plane_vectors': [], 'space_group': ''}
#    for vec in struc_dict['lattice']['matrix']:
#        # Is this sufficient to find all the in-plane vectors?
#        if vec[2] == 0 or (struc.pbc[2] and (vec[0] + vec[1]) > 0):
#            plane_vectors['plane_vectors'].append(vec[:2])
#    # finally check if setting of plane_vectors worked
#    if 'plane_vectors' not in plane_vectors:
#        raise ValueError('Could not set "plane_vectors" in STM_pathfinder')
#
#    # Here we get the symmetry operations that are possible
#    symmetry_matrices = SpacegroupAnalyzer(py_struc).get_point_group_operations(cartesian=True)
#
#    plane_vectors['space_group'] = SpacegroupAnalyzer(py_struc).get_symmetry_dataset()['number']
#
#    # Here we get the symmetry rotations
#
#    supp_mat = []
#
#    # Get the affine representation of the matrices
#    for symmops in range(len(symmetry_matrices)):
#        supp_mat.append(np.array(SymmOp.as_dict(symmetry_matrices[symmops])['matrix'][:3]))
#
#    # Get only the rotation matrices and correct the numerical error
#    rot_mat = []
#
#    # Take only the matrices for the in-plane rotations
#    for elements in range(len(supp_mat)):
#        rot_mat.append(supp_mat[elements][0:2, 0:2])
#
#    # Sometimes it will happen that some rotation matrices projected to the 2D space will have the same representation.
#    # here we only take the unique ones
#    unique_matrices = []
#    for matrix in rot_mat:
#        if not any(np.array_equal(matrix, m) for m in unique_matrices):
#            unique_matrices.append(matrix)
#
#    # Round off the numerical error
#    for elements in range(len(unique_matrices)):
#        for rows in range(len(unique_matrices[elements])):
#            for cols in range(len(unique_matrices[elements][rows])):
#                unique_matrices[elements][rows][cols] = round(unique_matrices[elements][rows][cols])
#
#    return plane_vectors, unique_matrices


#@engine.calcfunction
#def STM_pathfinder_cf(host_structure):
#    """
#    Calcfunction that gives back the structural information of the film, and the symmetries of the system
#    """
#
#    struc_info, symm_matrices = STM_pathfinder(host_structure)
#
#    return struc_info, symm_matrices


##############################################################################
# lattice generation (function of lattice plot)


def lattice_generation(rot, vec, x_start, y_start, xmax, ymax):
    """

        inputs ::

        x_len  : int  : value to create points between - x and x.
        y_len  : int  : value to create points between - y and y.
        rot    : list : list of the rotation matrices given by the symmetry of the system.
        vec    : list : list containing the two Bravais vectors.
        start_x: int  : starting value for the lattice generation in the x direction
        start_y: int  : starting value for the lattice generation in the y direction

        return ::

        points_to_eliminate : list : list of list containing the (x,y) positions to NOT to
                                     be scanned (unsorted)
        points_to_scan      : list : list of list containing the (x,y) positions to BE
                                     scanned (unsorted)
        """
    
    # Here we create a grid  made of points which are the linear combination of the lattice vectors
    x_len = xmax # 2
    y_len = ymax # 2  # maybe there is a way to make this more efficient... USE the lattice vectors! check the longest and use it

    x_interval = [i for i in range(-x_len, x_len+1)]

    y_interval = [i for i in range(-y_len, y_len+1)]

    points_to_scan = []
    lattice_points = []
    points_to_eliminate = []

    #Generat the lattice point and check if they fit in scanning area that it's wanted

    for i in x_interval:
        for j in y_interval:
            p = [i * x + j * y for x, y in zip(vec[0], vec[1])]
            if p[0] < xmax and p[0] > -xmax and p[1] < ymax and p[1] > -ymax:
                lattice_points.append(p)
    
    # sort the lattice points based on their y value. This is not necessary but makes the visualization nicer
    lattice_points = sorted(lattice_points, key=lambda y:y[1])

    for points in lattice_points:
        
        # First check if only the identity exists, in that case every point need to be scanned. 
        if len(rot) == 1: 
            points_to_scan = lattice_points
        else:
            for sym in rot[1:]:

                sym_point = np.dot(sym.tolist(), points).tolist() # Generate the symmetrical point
                
                if points not in points_to_eliminate and points not in points_to_scan:
                    points_to_scan.append(points)
                if sym_point not in points_to_eliminate and sym_point not in points_to_scan:
                    points_to_eliminate.append(sym_point)
                
            
    return points_to_eliminate, points_to_scan


##############################################################################
# lattice plot


def lattice_plot(plane_vectors, symm_vec, symm_matrices, grid_length_x, grid_length_y, **kwargs):
    """
        Helper tool to plot the position that will be scanned in the submission of the
        kkr_STM_wc workchain

        inputs ::

        plane_vectors : list : list containing the Bravais vector of the 2D lattice
        symm_vec      : bool : Toggle to show or not the Bravais vectors in the plot
        symm_matrices : list : list of the point-group symmetries matrices of the system
        grid_length_x : int  : scanning distance in the x direction
        grid_length_y : int  : scanning distance in the y direction

        return ::

        None

    """
    
    cused = kwargs.get('cused', '#FDE725FF')
    cunused = kwargs.get('cunused', '#33638DFF') 
    clattice = kwargs.get('clattice', '#3CBB75FF')

    #from aiida_kkr.tools.tools_STM_scan import lattice_generation
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    origin = np.array([[0, 0], [0, 0]])
    # Generation of the points to plot
    unused, used = lattice_generation(symm_matrices, plane_vectors, 0, 0, grid_length_x, grid_length_y)

    # Plotting of the points
    for element in unused:
        plt.scatter(element[0], element[1], marker='s', s=130, c=cunused)

    for element in used:
        plt.scatter(element[0], element[1], marker='D', s=130, c=cused)

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
        plt.quiver(*origin, element[0], element[1], color=clattice, angles='xy', scale_units='xy', scale=1)

    legend_elements = [
        Line2D([0], [0], color=cunused, lw=2, label='Unscanned Sites', marker='s'),
        Line2D([0], [0], color=cused, lw=2, label='Scanned Sites', marker='D'),
        Line2D([0], [0], color=clattice, lw=2, label='Bravais lattice'),
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(0.75, -0.15))

    plt.title('Lattice plot and symmetry directions')
    plt.ylabel('y direction')
    plt.xlabel('x direction')
    #plt.xticks(np.arange(-grid_length, grid_length, float(plane_vectors[0][0])))
    #plt.set_cmap(cmap)
    plt.grid(linestyle='--')
    plt.show()


##############################################################################
# find linear combination coefficients


def find_linear_combination_coefficients(plane_vectors, vectors):
    from operator import itemgetter
    """This helper function takes the planar vectors and a list of vectors
       and return the coefficients in the base of the planar vectors

       inputs ::

       plane vectors: list : list of list of the form [[a_x, a_y], [b_x, b_y]]

       returns ::

       indices : list : list of list of the form [[int_1, int_1]...[int_n, int_n]]
                        the integers refers to how many times that specific vectors is
                        present in the linear combination r = int_x * a + int_2 * b (sorted
                        list in the end)

       """

    # Formulate the system of equations Ax = b
    A = np.vstack((plane_vectors[0], plane_vectors[1])).T
    # inverse of A matrix
    Ainv = np.matrix(A)**(-1)
    indices = []
    # loop over flattened list
    for vec in np.array(vectors).reshape(-1, 2):
        # get indices from x = A^{-1}.b
        index = np.array((Ainv * np.matrix(vec).transpose()).transpose()).reshape(2)
        indices.append(index)
    # make sure to have integer values
    indices = np.round(np.array(indices), 0)

    isort = indices[:, 0] + indices[:, 1] / (10 * max(indices[:, 0]))
    indices = indices[isort.argsort()]

    return indices

##############################################################################
# Paser tool for retrieving data from an AiiDA Group of STM calculations

def STM_real_space_parser(group_STM, energy_pos, N0, _DEBUG_=False):
    from masci_tools.util.constants import BOHR_A
    from time import time
    
    """ 
    Function for parsing the data from a group of calculations for the STM 
    This function returns three lists: a position list, containing two lists: one for the x and one for the y positions
    a list containing the rho and another containing rhe mu values for the retrieved data. 
    
    Inputs:
        group_STM  (AiiDAGroup): A group containing the data we want to parse out.
        energy_pos (int) : The energy value that is being investigated. The value refers at the position of such value
                           in the list containing the data in the inputs. 
        N0 (int) : The numbers of atoms included in the original cluster that must be escluded. 
    """
    
    # First retrieve the lattice constant of the system. 
    
    try:
        ret  = group_STMnodes[0].called[0].called[0].called[1].outputs.retrieved
    except:
        print('WARNING: error while reading the group, it is possible that no node is contained here')
    
    
    with ret.open('inputcard') as _f:
        read = _f.readlines()
    
    for line in read:
        if 'ALATBASIS=' in line:
            # Split the line into words or space-separated values
            parts = line.split()
            
            # Assuming the numerical value is after 'ALTABASIS', get the next element
            if len(parts) > 1:
                alat = parts[1]  # The value is immediately attached to the nam. 
                break  # there is only one such parameter, break the loop after is found. 
    
    alat_ang = (alat * BOHR_A)
        
    if _DEBUG_:
        t0 = time() # Show time only in debugging procedure.
        
    # Create the lists containing the positions and the values of the collected data.
    all_pos = [[], []]
    all_dat_summed_rho = []
    all_dat_summed_mu = []
    for node in tqdm(group_STM.nodes):
        
        with node.called[0].called[0].called[1].outputs.retrieved.open("kkrflex_atominfo") as _f:
            pos = np.loadtxt(_f, skiprows=3) * alat_ang # Retrieve the positions of the atoms in the impurity cluster. 
        
        try:
            dat = node.outputs.STM_dos_data_lmdos.get_y()[0][1]
        except:
            continue
    
        for i in [-1,1]:
            for j in [-1,1]:
                
                all_pos[0] += list(i*pos[N0:,0])
                all_pos[1] += list(j*pos[N0:,1])
                all_dat_summed_rho += list(abs((dat[::2, energy_pos]+dat[1::2, energy_pos])[N0:])) #Both spin channels are considered here
                all_dat_summed_mu  += list(abs((dat[::2, energy_pos]-dat[1::2, energy_pos])[N0:]))
                
    
    all_pos = np.array(all_pos)
    all_dat_summed_rho = np.array(all_dat_summed_rho)
    all_dat_summed_mu= np.array(all_dat_summed_mu)
    
    
    if _DEBUG_:
        print(time()-t0)
    
    return all_pos, all_dat_summed_rho, all_dat_summed_mu

##############################################################################
# Real space plotting function
def STM_real_space_plot(positions, data, R0=0, R1=10, _DEBUG_=False, **kwargs):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    from time import time
    
    """
    Function for the real space plotting of the system
    
    Inputs:
    positions (List) : List containing 2 lists having thex and y positions
    data (List) : List containing the actual data to be plotted
    R0 (int) : Internal radius to be excluded (In Angstrom)
    R1 (int) : Outer radius from which to take the mean to normalise the plot (In Angstrom)
    
    kwargs ={
    's' = 125,                          #size of the plotted voronoi cells 
    's1' = 105 ,                        #size of the marker of the  escluded values,
    'lw' = 0 ,                          #line width 
    'cmap' = 'seismic',                 #type of map for the plotting
    'label' = '$\Delta n$ (states/ev)', #name of the plot
    'marker' = 'h',                     #shape of the marker for the plotting 
    'linthresh' = 0.05,                 #resolution of the plotting scale
    'figsize' = 10,                     #Size of the figure
    'fontsize' = '20',                  #fontsize of the plot title
    'fontsize_label' = 25,              #fontsize of the labels for the x and y axis
            }
    
    """
    
    # Extraction of the values for the plotting 
    s = kwargs.get('s', 125)
    s1 = kwargs.get('s1', 105)
    lw = kwargs.get('lw', 0)
    cmap = kwargs.get('cmap', 'seismic')
    label = kwargs.get('label' , '$\Delta n$ (states/ev)')
    markers = kwargs.get('markers', 'h')
    linthresh = kwargs.get('linthresh', 0.05)
    figsize = kwargs.get('figsize', 10)
    fontsize = kwargs.get('fontsize', 20)
    fontsize_label = kwargs.get('fontsize_label', 25)
            
    
    if _DEBUG_:
        t0 = time()
        
    all_dat_aux = data
    all_dat = data[np.sqrt(positions[0]**2+positions[1]**2)>R1] 
    
    all_dat_aux[np.sqrt(positions[0]**2+all_pospositions[1]**2)<R0] = np.NaN # Set to NaN those values that we don't want to see.
    plt.figure(figsize=(figsize,figsize))
    plt.scatter(positions[0], positions[1], c=all_dat_aux-np.nanmean(all_dat), cmap=cmap
                           , s=s, norm=colors.SymLogNorm(linthresh=linthresh), lw=lw, marker=marker) 

    cl = plt.gci().get_clim()
    cl = max(abs(cl[0]), cl[1])
    plt.clim(-cl, cl)
    
    cbar = plt.colorbar(orientation='vertical', aspect = 25, shrink = 1, pad=0.01)
    for t in cbar.ax.get_yticklabels():
         t.set_fontsize(20)
    cbar.set_label(label = label, fontsize=fontsize)
    
    nan_mask = np.isnan(all_dat_aux)
    nan_positions = np.argwhere(nan_mask)
    
    for ps in nan_positions:
            x, y = all_pos[0][ps[0]], all_pos[1][ps[0]]  # Get the x, y position of NaN
            plt.scatter(x, y, marker=marker, s = s1, color ='k')
    plt.xlabel('x ($\AA$)', fontsize = fontsize_label)
    plt.ylabel('y ($\AA$)', fontsize = fontsize_label)
    
    if _DEBUG_:
        print(time()-t0)
    
    plt.show()
    
##############################################################################
# Plotting for the FT of the real space image of a STM scanning

def FT_QPI(positions, data, length, R0=10, R1=20, _DEBUG_=False, **kwargs):
    import matplotlib.colors as colors 
    from scipy.interpolate import griddata
    import matplotlib.patches as patches
    from time import time

    """
    Function for the plotting of the Fourier transformed image of the real space STM imgage
    
    Inputs : 
    
    Inputs:
    positions (List) : List containing 2 lists having thex and y positions
    data (List) : List containing the actual data to be plotted
    length (int) : Length of the vector to divide the BZ zone
    res_points (int) : number of points to resolve
    R0 (int) : Internal radius to be excluded (In Angstrom)
    R1 (int) : Outer radius from which to take the mean to normalise the plot (In Angstrom)
    
    kwargs = {
            
            'xlim' = 2,                  #Lim for the x axis
            'ylim' = 2,                  #Lim for the y axis
            'cmap' = 'viridis',          #Color map that is used
            'method' = 'cubic',          #Interpolation method to use for the generation of the grid 
            'figsize' = 10,              #Size of the figure
            'fontsize' = 25,             #Fontsize of the x and y label
            'res_points' = 100,          #Number of points to resolve
            'tick_begin' = 4 ,           #Where to start the ticks
            'tick_spacing' = 0.25,       #Spacing between the ticks
            'interpolation' = 'quadric', #Interpolation method for the plotting of the FT 
            }
    
    """
    
    if _DEBUG_:
        t0 = time()
    
    xlim = kwargs.get('xlim', 2)
    ylim = kwargs.get('ylim', 2)
    cmap = kwargs.get('cmap','viridis')      
    method = kwargs.get('method','cubic')
    figsize = kwargs.get('figsize', 10)
    fontsize = kwargs.get('fontsize',25)  
    res_points = kwargs.get('res_points',100)
    tick_begin = kwargs.get('tick_begin',4)           
    tick_spacing = kwargs.get('tick_spacing',0.25)
    interpolation = kwargs.get('interpolation','quadric')
    
    # Generate the full set of points using the symmetry of the system before doing the interpolation
    grid_x, grid_y = np.mgrid[-length:length:(res_points*1j), -length:length:(res_points*1j)]
    
    # Use the position points, and then convert them to a 2D array 
    aux_pos = positions.copy()
    p = np.stack((aux_pos[0], aux_pos[1]), axis=-1)
            
    #Reduce the dimensionality of the data sample
    aux_data = data.copy()
    
    background_mean = np.nanmean(aux_data[np.sqrt(aux_pos[0]**2+aux_pos[1]**2)>=R1])
    
    mean = np.mean(aux_data)
    norm_data = aux_data-background_mean
    
    norm_data[np.sqrt(aux_pos[0]**2+aux_pos[1]**2)<=R0] = 0
    
    #Fourier transform of the data
    grid_z = griddata(p, norm_data, (grid_x, grid_y), method=method)
    ft = np.fft.fftshift(np.fft.fft2(grid_z))
    
    # resolution for the first BZ 
    plt.figure(figsize=(figsize,figsize))
    k_res = (np.pi/length)*res_points
    plt.imshow(np.abs(ft), extent=(-k_res, k_res, -k_res, k_res), cmap=cmap, interpolation=interpolation)#norm=colors.SymLogNorm(linthresh=0.00000001))
    plt.xlabel('$k_{x} (\AA^{-1})$',fontsize=fontsize)
    plt.ylabel('$k_{y} (\AA^{-1})$',fontsize=fontsize)


    # Adjust padding for better visualization
    plt.gca().tick_params(axis='both', which='major', pad=10)
    plt.xticks(np.arange(-tick_begin, tick_begin+0.5, tick_spacing))
    plt.yticks(np.arange(-tick_begin, tick_begin+0.5, tick_spacing))

    plt.xlim(-xlim, xlim)
    plt.ylim(-ylim, ylim)
    plt.colorbar(orientation='vertical', aspect = 25, shrink = 0.8, pad=0.01, label = 'Intensity')
    
    if _DEBUG_:
        print(time()-t0)
    
    plt.show()
        

##############################################################################
