# -*- coding: utf-8 -*-
"""
This module contains helper functions and tools for the combine_imps_wc workchain
"""
import numpy as np
import tarfile
from aiida.engine import calcfunction
from aiida.orm import Dict, SinglefileData, load_node, Bool
from aiida.common import InputValidationError
from aiida.common.folders import SandboxFolder
from aiida_kkr.tools.tools_kkrimp import modify_potential
from aiida_kkr.calculations import VoronoiCalculation, KkrimpCalculation
from aiida_kkr.workflows import kkr_imp_sub_wc
from .imp_cluster_tools import (
    pos_exists_already, get_inplane_neighbor, get_scoef_single_imp, get_zimp, combine_clusters,
    create_combined_imp_info, create_combined_imp_info_cf
)
from aiida_kkr.tools import get_ldaumatrices, get_LDAU_initmatrices_dict

__copyright__ = (u'Copyright (c), 2020, Forschungszentrum Jülich GmbH, '
                 'IAS-1/PGI-1, Germany. All rights reserved.')
__license__ = 'MIT license, see LICENSE.txt file'
__version__ = '0.3.4'
__contributors__ = (u'Philipp Rüßmann, David Antognini Silva')

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

    from aiida.plugins import DataFactory
    import tempfile
    import os

    # Create a SinglefileData node
    SinglefileData = DataFactory('singlefile')
    out_potential_content = retrieved.get_object_content('out_potential')

    # Create a temporary file
    temp_dir = tempfile.gettempdir()
    temp_file = os.path.join(temp_dir, 'out_potential')
    with open(temp_file, 'w') as f:
        f.write(out_potential_content)

    # Create a SinglefileData node with the temporary file
    potfile_sfd = SinglefileData(temp_file)

    # Remove the temporary file
    os.remove(temp_file)

    return potfile_sfd


def extract_potfile_from_retrieved(retrieved):
    """
    get output potential single file data from retrieved files or reuse existing
    """

    # check if retrieved has already a single file data child with given link label
    children = [res.node for res in retrieved.get_outgoing(link_label_filter='create_potfile_sfd').all()]
    if len(children) > 0 and 'result' in children[0].outputs:
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
            with tempfolder.open('potential_combined', u'rb') as _f:
                output_potential_sfd_node = SinglefileData(file=_f)
            # add label and description
            output_potential_sfd_node.label = 'combined_potentials'
            output_potential_sfd_node.description = f'combined potential of imps {pot_imp1.uuid} and {pot_imp2.uuid}'

    # return the combined potential
    return output_potential_sfd_node


@calcfunction
def combine_potentials_cf(kickout_info, pot_imp1, pot_imp2, nspin_node):
    return combine_potentials(kickout_info, pot_imp1, pot_imp2, nspin_node)


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

    if imp1_has_ldau:
        for k, v in settings_LDAU1.items():  # pylint: disable=used-before-assignment
            if 'iatom' in k:
                iatom = int(k.split('=')[1])
                settings_LDAU_combined[f'iatom={iatom}'] = v

    if imp2_has_ldau:
        for k, v in settings_LDAU2.items():  # pylint: disable=used-before-assignment
            if 'iatom' in k:
                iatom = int(k.split('=')[1])
                if kickout_info['i_removed_from_1'] is not None:
                    noffset = kickout_info['Ncls1'] - len(kickout_info['i_removed_from_1'])
                else:
                    noffset = kickout_info['Ncls1']
                settings_LDAU_combined[f'iatom={iatom+noffset}'] = v

    if has_old_ldaupot1:
        settings_LDAU_combined['initial_matrices'] = {}
        settings_LDAU_combined['initial_matrices'] = get_LDAU_initmatrices_dict(txts_ldaumat1)  # pylint: disable=used-before-assignment

    if has_old_ldaupot2:
        # If the dictionary already exists, update it with the new values
        if 'initial_matrices' in settings_LDAU_combined:
            settings_LDAU_combined['initial_matrices'].update(get_LDAU_initmatrices_dict(txts_ldaumat2, noffset))
        else:
            # If the dictionary does not exist, create it with the new values
            settings_LDAU_combined['initial_matrices'] = {}
            settings_LDAU_combined['initial_matrices'] = get_LDAU_initmatrices_dict(txts_ldaumat2, noffset)

    return Dict(settings_LDAU_combined)
