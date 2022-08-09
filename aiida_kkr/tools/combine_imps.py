# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools for the combine_imps_wc workchain
"""
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import tarfile
from aiida.engine import calcfunction
from aiida.orm import Dict, SinglefileData, load_node, Bool
from aiida.common import InputValidationError
from aiida.common.folders import SandboxFolder
from aiida_kkr.tools.tools_kkrimp import modify_potential, create_scoef_array
from aiida_kkr.calculations import VoronoiCalculation, KkrimpCalculation
from aiida_kkr.workflows import kkr_imp_sub_wc
from masci_tools.io.common_functions import get_alat_from_bravais
from six.moves import range

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.3.2'
__contributors__ = (u'Philipp Rüßmann')

# activate debug writeout
debug = False


def get_host_structure(impurity_workflow_or_calc):
    """
    extract host structure from impurity
    """
    #TODO extract host parent no from input but take into account calculation of host GF from inside kkrimp full workflow
    print(
        f'This is line in the combine impurity tool files at:: /opt/aiida-kkr/aiida_kkr/tools for deburging the line',
        end=' '
    )
    print(f'impurity_workflow_or_calc: {impurity_workflow_or_calc}')
    if impurity_workflow_or_calc.process_class == KkrimpCalculation:
        host_parent = impurity_workflow_or_calc.inputs.host_Greenfunction_folder
        # Here 'impurity_workflow_or_calc.process_class== combine_imps_wc' occurs circular import with this present module
    elif impurity_workflow_or_calc.process_class.__name__ == 'combine_imps_wc':
        imp_sub_wc = impurity_workflow_or_calc.get_outgoing(node_class=kkr_imp_sub_wc).first().node
        kkr_imp_calc = imp_sub_wc.get_outgoing(node_class=KkrimpCalculation).all()[-1].node
        host_parent = kkr_imp_calc.inputs.host_Greenfunction_folder
    elif impurity_workflow_or_calc.process_class == kkr_imp_sub_wc:
        kkr_imp_calc = impurity_workflow_or_calc.get_outgoing(node_class=KkrimpCalculation).all()[-1].node
        host_parent = kkr_imp_calc.inputs.host_Greenfunction_folder
    elif 'remote_data' in impurity_workflow_or_calc.inputs:
        # this is the case if impurity_workflow_or_calc workflow is kkr_imp_sub
        host_parent = impurity_workflow_or_calc.inputs.remote_data
    elif 'remote_data_gf' in impurity_workflow_or_calc.inputs:
        host_parent = impurity_workflow_or_calc.inputs.remote_data_gf
    else:
        host_parent = impurity_workflow_or_calc.inputs.remote_data_host
    host_structure, _ = VoronoiCalculation.find_parent_structure(host_parent)

    return host_structure


def get_nspin(imp_scf_workflow):
    """
    get nspin value from imp_scf_workflow output
    """

    uuid_last_params = imp_scf_workflow.outputs.workflow_info.get_dict().get('last_params_nodeinfo').get('uuid')
    nspin = load_node(uuid_last_params).get_dict()['NSPIN']

    return nspin


@calcfunction
def make_potfile_sfd(**kwargs):
    """
    Make single file data from output potential in retrieved folder.
    Extract output potential from tarball first if necessary.

    kwargs should be a dict with a single entry:
      kwargs = {'linklabel': retrieved}
    """

    if len(kwargs.keys()) != 1:
        raise IOError('kwargs input should only have a single key value pair!')

    for key in kwargs.keys():
        retrieved = kwargs[key]

    with SandboxFolder() as tempfolder:
        # find path of tempfolder
        with tempfolder.open('.dummy', 'w') as dummyfile:
            tempfolder_path = dummyfile.name
            tempfolder_path = tempfolder_path.replace('.dummy', '')

        # extract output potential here
        tar_filenames = []
        if KkrimpCalculation._FILENAME_TAR in retrieved.list_object_names():
            # get path of tarfile
            with retrieved.open(KkrimpCalculation._FILENAME_TAR) as tf:
                tfpath = tf.name
            # extract file from tarfile of retrieved to tempfolder
            with tarfile.open(tfpath) as tf:
                tar_filenames = [ifile.name for ifile in tf.getmembers()]
                filename = KkrimpCalculation._OUT_POTENTIAL
                if filename in tar_filenames:
                    tf.extract(filename, tempfolder_path)  # extract to tempfolder

        # store as SingleFileData
        with tempfolder.open(KkrimpCalculation._OUT_POTENTIAL, 'rb') as potfile:
            potfile_sfd = SinglefileData(file=potfile)

        return potfile_sfd


def extract_potfile_from_retrieved(retrieved):
    """
    get output potential single file data from retrieved files or reuse existing
    """

    # check if retrieved has already a single file data child with given link label
    children = [res.node for res in retrieved.get_outgoing(link_label_filter='create_potfile_sfd').all()]
    if len(children) > 0:
        potfile_sfd = children[0].outputs.result
        print('take existing node')
    else:
        # create a new single file data node from output using calcfunction for data provenance
        potfile_sfd = make_potfile_sfd(create_potfile_sfd=retrieved)
        print('create node')
    return potfile_sfd


def get_nspin_and_pot(imp):
    """
    extract nspin value and impurty potential single file data
    """
    from aiida_kkr.workflows.kkr_imp_sub import kkr_imp_sub_wc

    if imp.process_class == KkrimpCalculation:
        nspin = imp.outputs.output_parameters['nspin']
        pot_imp = extract_potfile_from_retrieved(imp.outputs.retrieved)
    else:
        # find KKRimp scf sub workflows
        if imp.process_class == kkr_imp_sub_wc:
            imp_sub = imp
        else:
            imp_sub = imp.get_outgoing(node_class=kkr_imp_sub_wc).first().node

        nspin = get_nspin(imp_sub)

        # extract potential
        pot_imp = imp_sub.outputs.host_imp_pot

    return nspin, pot_imp


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
    _, i_removed_from_1 = pos_exists_already(clust1[:, :3], clust2_offset[0, :3], debug)
    if debug:
        print('i_removed_from_1:', i_removed_from_1)

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
        clust1[0][4] = zimp1
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
    imp_info_combined = Dict(dict={'imp_cls': cluster_combined, 'Zimp': zimp_combined, 'Rimp_rel': rimp_rel_combined})

    # kickout info (used later in cfreation of combined potential)
    kickout_info = Dict(
        dict={
            'i_removed_from_1': i_removed_from_1,
            'kickout_list': kickout_list,
            'Ncls1': len(clust1),
            'Ncls2': len(clust2),
            'Ncls_combined': len(cluster_combined)
        }
    )

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


# combine potentials calcfunction
def combine_potentials(kickout_info, pot_imp1, pot_imp2, nspin_node):

    # unpack kickout info
    kickout_list = kickout_info['kickout_list']
    i_removed_from_1 = kickout_info['i_removed_from_1']
    Ncls1 = kickout_info['Ncls1']
    Ncls2 = kickout_info['Ncls2']
    Ncls_combined = kickout_info['Ncls_combined']
    nspin = nspin_node.value
    if debug:
        print('kickout_list:', kickout_list)
        print('i_removed:', i_removed_from_1)
        print('params;', nspin, Ncls1, Ncls2, Ncls_combined)

    # create neworder_pot list
    neworder_pot = list(range(Ncls1))
    if i_removed_from_1 is not None:
        neworder_pot = [neworder_pot[i] for i in range(len(neworder_pot)) if i not in i_removed_from_1]
    if debug:
        print('neworder_pot:', neworder_pot)

    # add dummy lines which are replace with pot 2
    if i_removed_from_1 is not None:
        N0 = len(neworder_pot)
        N_add = Ncls_combined - N0
        replacepos = [i for i in range(N0 + 1, N0 + N_add + 1)]  # range(N0, N0+N_add)]
    else:
        N0 = len(neworder_pot)
        N_add = Ncls_combined - N0
        replacepos = [i for i in range(N0, N0 + N_add)]  # range(N0, N0+N_add)]

    neworder_pot += replacepos

    # prepare index of pot2 without kciked out positions
    index_pot2 = [i for i in list(range(Ncls2)) if i not in kickout_list]

    if debug:
        print('replacepos:', replacepos)
        print('index_pot2:', index_pot2)

    # create replacelist (mapping which positions of neworder_pos are taken from pot2 instead)
    if i_removed_from_1 is not None:
        replacelist_pot2 = [(replacepos[i] - 1, index_pot2[i]) for i in range(len(replacepos))]
    else:
        replacelist_pot2 = [(replacepos[i], index_pot2[i]) for i in range(len(replacepos))]

    # take care of spin doubling for NSPIN==2
    if nspin > 1:
        neworder_pot = np.array([[2 * i, 2 * i + 1] for i in neworder_pot]).reshape(-1)
        replacelist_pot2 = np.array([[(2 * i[0], 2 * i[1]), (2 * i[0] + 1, 2 * i[1] + 1)] for i in replacelist_pot2]
                                    ).reshape(-1, 2)

    # now combine potentials
    with SandboxFolder() as tempfolder:
        with tempfolder.open('potential_combined', 'w') as out_pot_fhandle:
            with pot_imp1.open(pot_imp1.filename) as pot1_filehande:
                with pot_imp2.open(pot_imp2.filename) as pot2_filehande:
                    # use neworder_potential function
                    modify_potential().neworder_potential(
                        pot1_filehande,
                        out_pot_fhandle,
                        neworder_pot,
                        potfile_2=pot2_filehande,
                        replace_from_pot2=replacelist_pot2,
                        debug=debug
                    )

            # store output potential to SinglefileData
            output_potential_sfd_node = SinglefileData(file=tempfolder.open('potential_combined', u'rb'))
            # add label and description
            output_potential_sfd_node.label = 'combined_potentials'
            output_potential_sfd_node.description = f'combined potential of imps {pot_imp1.uuid} and {pot_imp2.uuid}'

    # return the combined potential
    return output_potential_sfd_node


@calcfunction
def combine_potentials_cf(kickout_info, pot_imp1, pot_imp2, nspin_node):
    return combine_potentials(kickout_info, pot_imp1, pot_imp2, nspin_node)


def get_ldaumatrices(retrieved):
    """
    Take retrieved folder of KkrimpCalculation and extract ldaupot file
    If it exists we extract the LDAU matrices 'wldau', 'uldau' and 'phi'
    """
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
                txt_dict_ldaumats['wldau'] = txt[iwldau + 2:iuldau]
                txt_dict_ldaumats['uldau'] = txt[iuldau + 2:iphi]
                txt_dict_ldaumats['phi'] = txt[iphi + 2:]

    return has_ldaupot_file, txt_dict_ldaumats


@calcfunction
def combine_settings_ldau(**kwargs):
    """
    combine LDA+U settings using information from kickout_info to correct the atom index of the second impurity
    """

    imp1_has_ldau, has_old_ldaupot1 = False, False
    imp2_has_ldau, has_old_ldaupot2 = False, False

    if 'settings_LDAU1' in kwargs:
        imp1_has_ldau = True
        settings_LDAU1 = kwargs['settings_LDAU1'].get_dict()
        # get initial matrices from retrieved if given in input
        if 'retrieved1' in kwargs:
            has_old_ldaupot1, txts_ldaumat1 = get_ldaumatrices(kwargs['retrieved1'])

    if 'settings_LDAU2' in kwargs:
        imp2_has_ldau = True
        settings_LDAU2 = kwargs['settings_LDAU2'].get_dict()
        # get initial matrices from retrieved if given in input
        if 'retrieved2' in kwargs:
            has_old_ldaupot2, txts_ldaumat2 = get_ldaumatrices(kwargs['retrieved2'])

    if 'kickout_info' in kwargs:
        kickout_info = kwargs['kickout_info'].get_dict()
    else:
        raise KeyError('Need to have kickout_info key value pair in input.')

    # now combine LDAU settings
    settings_LDAU_combined = {}

    if has_old_ldaupot1 or has_old_ldaupot2:
        settings_LDAU_combined['initial_matrices'] = {}

    if imp1_has_ldau:
        for k, v in settings_LDAU1.items():
            if 'iatom' in k:
                iatom = int(k.split('=')[1])
                # TODO: implement something for the case when LDAU is not only on the impurity site at iatom==0
                settings_LDAU_combined[f'iatom={iatom}'] = v
                if has_old_ldaupot1:
                    settings_LDAU_combined['initial_matrices'][f'iatom={iatom}'] = txts_ldaumat1

    if imp2_has_ldau:
        for k, v in settings_LDAU2.items():
            if 'iatom' in k:
                iatom = int(k.split('=')[1])
                if kickout_info['i_removed_from_1'] is not None:
                    noffset = kickout_info['Ncls1'] - len(kickout_info['i_removed_from_1'])
                else:
                    noffset = kickout_info['Ncls1']
                settings_LDAU_combined[f'iatom={iatom+noffset}'] = v
                if has_old_ldaupot2:
                    settings_LDAU_combined['initial_matrices'][f'iatom={iatom+noffset}'] = txts_ldaumat2

    return Dict(dict=settings_LDAU_combined)
