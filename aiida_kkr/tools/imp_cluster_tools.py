# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools for manipulating the impurity info node
"""

import numpy as np
from aiida.engine import calcfunction
from aiida.orm import Dict
from aiida_kkr.tools.tools_kkrimp import create_scoef_array
from masci_tools.io.common_functions import get_alat_from_bravais

__copyright__ = (u'Copyright (c), 2023, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.1.0'
__contributors__ = (u'Philipp Rüßmann')

# activate debug writeout
debug = False


def get_zimp(impurity_info):
    """
    extract zimp from impurity_info node and return as list
    """
    zimp = impurity_info.get_dict().get('Zimp')

    try:
        zimp = list(zimp)
    except TypeError:
        zimp = [zimp]

    return zimp


def get_scoef_single_imp(host_structure, impinfo_node):
    """
    Create scoef cluster from host_structure and impurity info nodes

    :para impinfo_node: Impurity info node (should have at least Rcut and ilayer_center in dict)
    :type impinfo_node: aiida.orm.Dict
    :para host_structure: structure of host crystal into which impurities are embedded
    :type host_structure: aiida.orm.StructureData

    :return: scoef array (positions [x,y,z], layer index, distance to first position in imp cluster)
    :type: numpy.array
    """
    impinfo = impinfo_node.get_dict()
    Rcut = impinfo.get('Rcut', None)
    hcut = impinfo.get('hcut', -1.)
    cylinder_orient = impinfo.get('cylinder_orient', [0., 0., 1.])
    ilayer_center = impinfo.get('ilayer_center', 0)

    clust = create_scoef_array(host_structure, Rcut, hcut, cylinder_orient, ilayer_center)
    # sort after distance
    clust = clust[(clust[:, -1]).argsort()]

    return clust


def get_inplane_neighbor(host_structure, i_neighbor, r_out_of_plane):
    """
    create in-plane neighbor
    """
    # if offset of second impurity is given in units of the host structure bravais vectors:
    cell = np.array(host_structure.cell)

    # find and sort list of nearest neighbors
    dist = []
    icount = 0
    N_neighbor_search = 5  # define box in which neighbors are searched
    for j in list(range(N_neighbor_search + 1)
                  ) + list(range(-N_neighbor_search, 0)[::-1]):  # use this to have better ordering (nicer numbers)
        for i in list(range(N_neighbor_search + 1)) + list(range(-N_neighbor_search, 0)[::-1]):
            # create position
            r = cell[0] * i + cell[1] * j + r_out_of_plane
            d = np.sqrt(np.sum(r**2))
            # check if already in list
            add_pair = False
            if icount > 1:
                if abs(np.array(dist)[:, -1] - d).min() > 10**-5:
                    add_pair = True
            else:
                add_pair = True
            # add distance
            if add_pair:
                dist.append([i, j, r[0], r[1], r[2], d])
            icount += 1

    alat = get_alat_from_bravais(cell, host_structure.pbc[2])
    dist = np.array(dist)
    dist[:, 2:] = dist[:, 2:] / alat
    dist_sort_list = dist[:, -1].argsort()

    # find in-plane neighbor
    r_offset = dist[dist_sort_list[i_neighbor]][2:5]  # element 0 is no offset

    return r_offset


def pos_exists_already(pos_list_old, pos_new, debug=False):
    """
    check if pos_new is in pos_list_old
    """
    sort_ref = np.array(pos_list_old)
    dists = np.sqrt(np.sum((sort_ref - np.array(pos_new))**2, axis=1))
    mask = (dists < 10**-5)

    if debug:
        print('check pos', dists.min())

    if dists.min() < 10**-5:
        return True, [i for i in range(len(mask)) if mask[i]]
    else:
        return False, None


def combine_clusters(clust1, clust2_offset, single_single, debug=False):
    """
    combine impurity clusters and remove doubles in there

    :return cluster_combined:
    :return rimp_rel_combined: relative position of impurities in cluster
    :return kickout_list: licst of positions doubled in clust 2
    :return i_neighbor_inplane: position removed from cluster 1 (doubled with impurity position)
    """

    # now combine cluster1 and cluster2
    cluster_combined = list(clust1.copy())

    # check if imp position of cluster 2 is inside and remove that position
    # this ensures that imp2 is not kicked out
    found_double, i_removed_from_1 = pos_exists_already(clust1[:, :3], clust2_offset[0, :3], debug)
    if debug:
        print('i_removed_from_1:', i_removed_from_1, found_double)

    # remove doubled position from impcls1 if there is any
    if i_removed_from_1 is not None:
        #removed = cluster_combined.pop(i_removed_from_1)
        cluster_combined = [clust1[i] for i in range(len(clust1)) if i not in i_removed_from_1]

    # add the rest of the imp cluster of imp 2 if position does not exist already in cluster of imp 1
    kickout_list = []
    for ipos, pos_add in enumerate(clust2_offset):
        if ipos == 0 or not pos_exists_already(np.array(cluster_combined)[:, :3], pos_add[:3], debug)[0]:
            cluster_combined.append(pos_add)
        else:
            kickout_list.append(ipos)
    cluster_combined = np.array(cluster_combined)

    # fix distances for second half
    cluster_combined[:, -1] = np.sqrt(np.sum(cluster_combined[:, :3]**2, axis=1))

    # construct Rimp_rel list
    rimp_rel_combined = [clust1[0, :3]] + [clust2_offset[0, :3]]

    return cluster_combined, rimp_rel_combined, kickout_list, i_removed_from_1


def create_combined_imp_info(
    host_structure, impinfo1, impinfo2, offset_imp2, imps_info_in_exact_cluster, single_single, debug=False
):
    """
    create impurity clusters from impinfo nodes and combine these putting the second
    impurity to the i_neighbor_inplane-th in-plane neighbor
    """
    #TODO: Try it with the self object in the in the input

    imps_info_in_exact_cluster = imps_info_in_exact_cluster.get_dict()
    if single_single:
        zimp1 = imps_info_in_exact_cluster['Zimps'][0]
    zimp2 = imps_info_in_exact_cluster['Zimps'][-1]

    if 'imp_cls' in impinfo1.get_dict():
        clust1 = np.array(impinfo1['imp_cls'])

    elif single_single:
        # create cluster of imp1
        clust1 = get_scoef_single_imp(host_structure, impinfo1)

    # do the same for imp2
    clust2 = get_scoef_single_imp(host_structure, impinfo2)

    # set zimp in scoef file (not used by the code but makes it easier to read the files / debug)
    if single_single:
        clust1[0][4] = zimp1  # pylint: disable=possibly-used-before-assignment
    clust2[0][4] = zimp2
    #if debug:
    #    print('cls1:', clust1)
    #    print('cls2:', clust2)

    if 'r_offset' in offset_imp2.get_dict():
        # use offset given in input
        r_offset = offset_imp2['r_offset']
    else:
        # find offset taking into account the possible out-of-plane vector if the imps are in different layers
        r_out_of_plane = np.array([0, 0, 0])
        center_imp = imps_info_in_exact_cluster['ilayers'][0]
        layer2 = imps_info_in_exact_cluster['ilayers'][-1]
        if center_imp != layer2:
            pos1 = np.array(host_structure.sites[center_imp].position)
            pos2 = np.array(host_structure.sites[layer2].position)
            r_out_of_plane = pos2 - pos1
        i_neighbor_inplane = imps_info_in_exact_cluster['offset_imps'][-1]
        r_offset = get_inplane_neighbor(host_structure, i_neighbor_inplane, r_out_of_plane)
    if debug:
        print('r_offset:', r_offset)

    # add offset to cluster 2
    clust2_offset = clust2.copy()
    clust2_offset[:, :3] += r_offset

    cluster_combined, rimp_rel_combined, kickout_list, i_removed_from_1 = combine_clusters(
        clust1, clust2_offset, single_single, debug
    )

    if 'Rimp_rel' in impinfo1.get_dict():
        rimp_rel_combined = impinfo1['Rimp_rel'] + rimp_rel_combined[1:]

    if debug:
        #print('cls_combined:', cluster_combined)
        print('rimp_rel_combined:', rimp_rel_combined)
        print('kickout_list:', kickout_list)
        print('i_removed_from_1:', i_removed_from_1)

    zimp_combined = imps_info_in_exact_cluster['Zimps']

    # create new imp_info node with imp_cls, Rimp_rel and Zimp definig the cluster and impurity location
    imp_info_combined = Dict({'imp_cls': cluster_combined, 'Zimp': zimp_combined, 'Rimp_rel': rimp_rel_combined})

    # kickout info (used later in cfreation of combined potential)
    kickout_info = Dict({
        'i_removed_from_1': i_removed_from_1,
        'kickout_list': kickout_list,
        'Ncls1': len(clust1),
        'Ncls2': len(clust2),
        'Ncls_combined': len(cluster_combined)
    })

    return {'imp_info_combined': imp_info_combined, 'kickout_info': kickout_info}


@calcfunction
def create_combined_imp_info_cf(
    host_structure, impinfo1, impinfo2, offset_imp2, imps_info_in_exact_cluster, single_single
):
    """
    create impurity clusters from impinfo nodes and combine these putting the second
    impurity to the i_neighbor_inplane-th in-plane neighbor
    """

    return create_combined_imp_info(
        host_structure, impinfo1, impinfo2, offset_imp2, imps_info_in_exact_cluster, single_single
    )
