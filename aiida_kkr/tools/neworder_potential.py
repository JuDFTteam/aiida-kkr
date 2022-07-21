# -*- coding: utf-8 -*-
"""
Helper calcfunction that allws to modify a potential file
"""

from aiida.engine import calcfunction
from aiida.common.exceptions import InputValidationError
from aiida_kkr.tools.context import open_context_to_stack
from contextlib import ExitStack


@calcfunction
def neworder_potential_wf(settings_node, parent_calc_folder, **kwargs):
    """
    Workfunction to create database structure for aiida_kkr.tools.modify_potential.neworder_potential function
    A temporary file is written in a Sandbox folder on the computer specified via
    the input computer node before the output potential is stored as SinglefileData
    in the Database.

    :param settings_node: settings for the neworder_potential function (Dict)
    :param parent_calc_folder: parent calculation remote folder node where the input
        potential is retreived from (RemoteData)
    :param parent_calc_folder2: *optional*, parent calculation remote folder node where
        the second input potential is retreived from in case 'pot2' and 'replace_newpos'
        are also set in settings_node (RemoteData)
    :param debug: *optional*, contol wether or not debug information is written out (aiida.orm.Bool)

    :returns: output_potential node (SinglefileData)

    .. note::

        The settings_node dictionary needs to be of the following form::

            settings_dict = {'neworder': [list of intended order in output potential]}

        Optional entries are::

            'out_pot': '<filename_output_potential>'  name of the output potential file, defaults to 'potential_neworder' if not specified
            'pot1': '<filename_input_potential>'      if not given we will try to find it from the type of the parent remote folder
            'pot2': '<filename_second_input_file>'
            'replace_newpos': [[position in neworder list which is replace with potential from pot2, position in pot2 that is chosen for replacement]]
            'switch_spins': [indices of atom for which spins are exchanged] (indices refer to position in neworder input list)
            'label': 'label_for_output_node'
            'description': 'longer_description_for_output_node'
    """
    import os
    from aiida_kkr.tools.tools_kkrimp import modify_potential
    from aiida.common.folders import SandboxFolder
    from aiida.common.exceptions import UniquenessError
    from aiida.orm import CalcJobNode, Dict, RemoteData, SinglefileData

    if 'debug' in list(kwargs.keys()):
        debug = kwargs.get('debug').value
    else:
        debug = False

    if 'parent_calc_folder2' in list(kwargs.keys()):
        parent_calc_folder2 = kwargs.get('parent_calc_folder2', None)
    else:
        parent_calc_folder2 = None

    # check input consistency
    if not isinstance(settings_node, Dict):
        raise InputValidationError('settings_node needs to be a valid aiida Dict node')
    if not isinstance(parent_calc_folder, RemoteData):
        raise InputValidationError('parent_calc_folder needs to be a valid aiida RemoteData node')
    if parent_calc_folder2 is not None and not isinstance(parent_calc_folder2, RemoteData):
        raise InputValidationError('parent_calc_folder2 needs to be a valid aiida RemoteData node')

    settings_dict = settings_node.get_dict()
    pot1 = settings_dict.get('pot1', None)
    if pot1 is None:
        # try to extract the potential name from the type of the parent_calc_folder
        try:
            pot1 = extract_potname_from_remote(parent_calc_folder)
        except ValueError:
            raise InputValidationError(
                'settings_node_dict needs to have key "pot1" containing the filename of the input potential'
            )
    out_pot = settings_dict.get('out_pot', 'potential_neworder')
    neworder = settings_dict.get('neworder', None)
    if neworder is None:
        raise InputValidationError(
            'settings_node_dict needs to have key "neworder" containing the list of new positions'
        )
    pot2 = settings_dict.get('pot2', None)
    if pot2 is None and parent_calc_folder2 is not None:
        # try to extract the potential name from the type of the parent_calc_folder
        try:
            pot2 = extract_potname_from_remote(parent_calc_folder2)
        except ValueError:
            raise InputValidationError(
                'settings_node_dict needs to have key "pot2" containing the filename of the input potential'
            )
    replace_newpos = settings_dict.get('replace_newpos', None)
    switch_spins = settings_dict.get('switch_spins', [])

    # Create Sandbox folder for generation of output potential file
    # and construct output potential
    with SandboxFolder() as tempfolder:
        # Get abolute paths of input files from parent calc and filename
        parent_calcs = parent_calc_folder.get_incoming(node_class=CalcJobNode).all()
        n_parents = len(parent_calcs)
        if n_parents != 1:
            raise UniquenessError(
                'Input RemoteData is child of {} '
                'calculation{}, while it should have a single parent'
                ''.format(n_parents, '' if n_parents == 0 else 's')
            )
        parent_calc = parent_calcs[0].node

        # extract nspin from parent calc's input parameter node
        nspin = parent_calc.inputs.parameters.get_dict().get('NSPIN')
        neworder_spin = []
        ii = 0
        for iatom in neworder:
            spins = range(nspin)
            # change spin order if needed
            if ii in switch_spins:
                spins = spins[::-1]
            for ispin in spins:
                neworder_spin.append(iatom * nspin + ispin)
            ii += 1
        neworder = neworder_spin

        # Copy optional files?
        if pot2 is not None and parent_calc_folder2 is not None:
            parent_calcs = parent_calc_folder2.get_incoming(node_class=CalcJobNode).all()
            n_parents = len(parent_calcs)
            if n_parents != 1:
                raise UniquenessError(
                    'Input RemoteData of parent_calc_folder2 is child of {} '
                    'calculation{}, while it should have a single parent'
                    ''.format(n_parents, '' if n_parents == 0 else 's')
                )
            else:
                parent_calc2 = parent_calcs[0].node
            if pot2 not in parent_calc2.outputs.retrieved.list_object_names():
                raise InputValidationError(
                    'neworder_potential_wf: pot2 does not exist', pot2,
                    parent_calc.outputs.retrieved.list_object_names()
                )

        # open files in context (if name is None the fhandle will just be None and then ignored)
        with ExitStack() as stack:
            out_pot_fhandle = open_context_to_stack(stack, tempfolder, out_pot, u'w')
            pot1_fhandle = open_context_to_stack(stack, parent_calc.outputs.retrieved, pot1)
            pot2_fhandle = None
            if pot2 is not None:
                pot2_fhandle = open_context_to_stack(stack, parent_calc2.outputs.retrieved, pot2)
            # run neworder_potential function
            modify_potential().neworder_potential(
                pot1_fhandle,
                out_pot_fhandle,
                neworder,
                potfile_2=pot2_fhandle,
                replace_from_pot2=replace_newpos,
                debug=debug
            )

        # store output potential to SinglefileData
        with tempfolder.open(out_pot, u'rb') as fhandle:
            output_potential_sfd_node = SinglefileData(file=fhandle)

        lbl = settings_dict.get('label', None)
        if lbl is not None:
            output_potential_sfd_node.label = lbl
        desc = settings_dict.get('description', None)
        if desc is not None:
            output_potential_sfd_node.description = desc

        # TODO create shapefun sfd node accordingly
        """
        out_shape_path =

        output_shapefun_sfd_node = SinglefileData(file=out_shape_path)

        lbl2 = settings_dict.get('label_shape', None)
        if lbl2 is None and lbl is not None:
            lbl2 = lbl
        if lbl2 is not None:
            output_shapefun_sfd_node.label = lbl2
        desc2 = settings_dict.get('description_shape', None)
        if desc2 is None and desc is not None:
            desc2 = desc
        if desc2 is not None:
            output_shapefun_sfd_node.description = desc2

        return output_potential_sfd_node, output_shapefun_sfd_node
        """
        return output_potential_sfd_node


def extract_potname_from_remote(parent_calc_folder):
    """
    extract the bname of the output potential from a RemoteData folder
    """
    from aiida_kkr.calculations import KkrCalculation, VoronoiCalculation
    from aiida.orm import CalcJobNode

    pot_name = None
    # extract list of parents (can only extract the parent calculation
    # if there is only a single incoming link to follow)
    parents = parent_calc_folder.get_incoming(node_class=CalcJobNode)
    if len(list(parents)) == 1:
        parent = parents.first().node
        # now extract the pot_name dependeing on the parent calculation's type
        if parent.process_class == KkrCalculation:
            pot_name = KkrCalculation._OUT_POTENTIAL
        elif parent.process_class == VoronoiCalculation:
            pot_name = VoronoiCalculation._OUT_POTENTIAL_voronoi

    # return the potential name or raise an error if nothing was found
    if pot_name is not None:
        return pot_name
    else:
        raise ValueError('Could not extract a potential name')
