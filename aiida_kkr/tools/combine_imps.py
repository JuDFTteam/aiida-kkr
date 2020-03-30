# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools for the combine_imps_wc workchain
"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
from aiida.engine import calcfunction
from aiida.orm import Dict, SinglefileData, load_node
from aiida.common import InputValidationError
from aiida.common.folders import SandboxFolder
from aiida_kkr.tools.tools_kkrimp import modify_potential, create_scoef_array
from aiida_kkr.calculations import VoronoiCalculation
from masci_tools.io.common_functions import get_alat_from_bravais
from six.moves import range


__copyright__ = (u"Copyright (c), 2020, Forschungszentrum Jülich GmbH, "
                 "IAS-1/PGI-1, Germany. All rights reserved.")
__license__ = "MIT license, see LICENSE.txt file"
__version__ = "0.1.0"
__contributors__ = (u"Philipp Rüßmann")


def get_host_structure(impurity_workflow):
    """
    extract host structure from impurity
    """
    #TODO extract host parent no from input but take into account calculation of host GF from inside kkrimp full workflow

    if 'remote_data_gf' in impurity_workflow.inputs:
        host_parent = impurity_workflow.inputs.remote_data_gf
    else:
        host_parent = impurity_workflow.inputs.remote_data_host 
    host_structure, _ = VoronoiCalculation.find_parent_structure(host_parent)

    return host_structure


def get_nspin(imp_scf_workflow):
    """
    get nspin value from imp_scf_workflow output
    """

    uuid_last_params = imp_scf_workflow.outputs.workflow_info.get_dict().get('last_params_nodeinfo').get('uuid')
    nspin = load_node(uuid_last_params).get_dict()['NSPIN']

    return nspin


# combine clusters calcfunction

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
    clust = clust[(clust[:,-1]).argsort()]
    
    return clust

def get_inplane_neighbor(host_structure, i_neighbor):
    """
    create in-plane neighbor
    """
    # if offset of second impurity is given in units of the host structure bravais vectors:
    cell = np.array(host_structure.cell)

    # find and sort list of nearest neighbors
    dist =[]
    icount = 0
    N_neighbor_search = 5 # define box in which neighbors are searched
    for j in list(range(N_neighbor_search+1))+list(range(-N_neighbor_search,0)[::-1]): # use this to have better ordering (nicer numbers)
        for i in list(range(N_neighbor_search+1))+list(range(-N_neighbor_search,0)[::-1]):
            # create position
            r = cell[0]*i + cell[1]*j
            d = np.sqrt(np.sum(r**2))
            # check if already in list
            add_pair = False
            if icount>1:
                if abs(np.array(dist)[:,-1]-d).min()>10**-5:
                    add_pair = True
            else:
                add_pair = True
            # add distance
            if add_pair:
                dist.append([i, j, r[0], r[1], r[2], d])
            icount += 1

    alat = get_alat_from_bravais(cell, host_structure.pbc[2])
    dist = np.array(dist)
    dist[:,2:] = dist[:,2:] / alat
    dist_sort_list = dist[:,-1].argsort()

    # find in-plane neighbor
    r_offset =dist[dist_sort_list[i_neighbor]][2:5] # element 0 is no offset

    return r_offset

def pos_exists_already(pos_list_old, pos_new):
    """
    check if pos_new is in pos_list_old
    """
    sort_ref = np.array(pos_list_old)
    dists = np.sqrt(np.sum((sort_ref-np.array(pos_new))**2, axis=1))
    
    if dists.min() < 10**-5:
        return True, dists.argsort()[dists < 10**-5]
    else:
        return False, None

def combine_clusters(clust1, clust2_offset):
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
    pos_doubled, index = pos_exists_already(clust1[:,:3], clust2_offset[0,:3])
    if pos_doubled:
        i_removed_from_1 = index
    else:
        i_removed_from_1 = None
    # remove doubled position from impcls1 if there is any
    if i_removed_from_1 is not None:
        #removed = cluster_combined.pop(i_removed_from_1)
        cluster_combined = [clust1[i] for i in range(len(clust1)) if i not in i_removed_from_1]

    # add the rest of the imp cluster of imp 2 if position does not exist already in cluster of imp 1
    kickout_list = []
    ipos = 0
    for pos_add in clust2_offset:
        if ipos == 0 or not pos_exists_already(np.array(cluster_combined)[:,:3], pos_add[:3])[0]:
            cluster_combined.append(pos_add)
        else:
            kickout_list.append(ipos)
        ipos += 1
    cluster_combined = np.array(cluster_combined)

    # fix distances for second half
    cluster_combined[:,-1] = np.sqrt(np.sum(cluster_combined[:,:3]**2, axis=1))

    # construct Rimp_rel list
    rimp_rel_combined = [clust1[0,:3]] + [clust2_offset[0,:3]]

    return cluster_combined, rimp_rel_combined, kickout_list, i_removed_from_1

#@calcfunction
def create_combined_imp_info_cf(host_structure, impinfo1, impinfo2, offset_imp2):
    """
    create impurity clusters from impinfo nodes and combine these putting the second
    impurity to the i_neighbor_inplane-th in-plane neighbor
    """

    i_neighbor_inplane = offset_imp2['index']

    zimp1 = get_zimp(impinfo1)
    zimp2 = get_zimp(impinfo2)

    # combine Zimp lists
    zimp_combined = zimp1 + zimp2

    # create cluster of imp1
    clust1 = get_scoef_single_imp(host_structure, impinfo1)

    # do the same for imp2
    clust2 = get_scoef_single_imp(host_structure, impinfo2)

    # set zimp in scoef file (not used by the code but makes it easier to read the files / debug)
    clust1[0][4] = zimp1[0]
    clust2[0][4] = zimp2[0]

    # find offset
    #TODO allow also out-of-plane neighbors
    r_offset = get_inplane_neighbor(host_structure, i_neighbor_inplane)

    # add offset to cluster 2
    clust2_offset = clust2.copy()
    clust2_offset[:, :3] += r_offset

    cluster_combined, rimp_rel_combined, kickout_list, i_removed_from_1 = combine_clusters(clust1, clust2_offset)

    # create new imp_info node with imp_cls, Rimp_rel and Zimp definig the cluster and impurity location
    imp_info_combined = Dict(dict={'imp_cls': cluster_combined, 'Zimp': zimp_combined, 'Rimp_rel': rimp_rel_combined})
    
    # kickout info (used later in cfreation of combined potential)
    kickout_info = Dict(dict={'i_removed_from_1': i_removed_from_1, 'kickout_list': kickout_list, 
                              'Ncls1': len(clust1), 'Ncls2': len(clust2), 'Ncls_combined': len(cluster_combined)}
                        )

    return {'imp_info_combined': imp_info_combined, 'kickout_info': kickout_info}




# combine potentials calcfunction

@calcfunction
def combine_potentials_cf(kickout_info, pot_imp1, pot_imp2, nspin_node):

    # unpack kickout info
    kickout_list = kickout_info['kickout_list']
    i_removed_from_1 = kickout_info['i_removed_from_1']
    Ncls1 = kickout_info['Ncls1']
    Ncls2 = kickout_info['Ncls2']
    Ncls_combined = kickout_info['Ncls_combined']
    nspin = nspin_node.value

    # create neworder_pot list
    neworder_pot = list(range(Ncls1))
    if i_removed_from_1 is not None:
        removed = neworder_pot.pop(i_removed_from_1)
        print(('removed: ', removed))

    # add dummy lines which are replace with pot 2
    N0 = len(neworder_pot)
    N_add = Ncls_combined-N0
    replacepos = [i for i in range(N0+1, N0+N_add+1)]
    neworder_pot += replacepos

    # prepare index of pot2 without kciked out positions
    index_pot2 = [i for i in list(range(Ncls2)) if i not in kickout_list]

    # create replacelist (mapping which positions of neworder_pos are taken from pot2 instead)
    replacelist_pot2 = [(replacepos[i]-1, index_pot2[i]) for i in range(len(replacepos))]

    # take care of spin doubling for NSPIN==2
    if nspin>1:
        neworder_pot = np.array([[2*i, 2*i+1] for i in neworder_pot]).reshape(-1)
        replacelist_pot2 = np.array([[(2*i[0], 2*i[1]), (2*i[0]+1, 2*i[1]+1)] for i in replacelist_pot2]).reshape(-1, 2)


    # now combine potentials
    with SandboxFolder() as tempfolder:
        with tempfolder.open('potential_combined', 'w') as out_pot_fhandle:
            with pot_imp1.open(pot_imp1.filename) as pot1_filehande:
                with pot_imp2.open(pot_imp2.filename) as pot2_filehande:
                    # use neworder_potential function
                    modify_potential().neworder_potential(pot1_filehande, out_pot_fhandle, neworder_pot, potfile_2=pot2_filehande,
                                                        replace_from_pot2=replacelist_pot2)

            # store output potential to SinglefileData
            output_potential_sfd_node = SinglefileData(file=tempfolder.open('potential_combined', u'rb'))
            # add label and description
            output_potential_sfd_node.label = 'combined_potentials'
            output_potential_sfd_node.description = 'combined potential of imps {} and {}'.format(pot_imp1.uuid, pot_imp2.uuid)

    # return the combined potential
    return output_potential_sfd_node