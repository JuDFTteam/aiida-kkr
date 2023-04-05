# -*- coding: utf-8 -*-
"""
Here we specify some defaults for cli commands
"""


# Structures
def get_si_bulk_structure():
    """Return a `StructureData` representing bulk silicon.

    The database will first be queried for the existence of a bulk silicon crystal. If this is not the case, one is
    created and stored. This function should be used as a default for CLI options that require a `StructureData` node.
    This way new users can launch the command without having to construct or import a structure first. This is the
    reason that we hardcode a bulk silicon crystal to be returned. More flexibility is not required for this purpose.

    :return: a `StructureData` representing bulk silicon
    """
    from ase.spacegroup import crystal
    from aiida.orm import QueryBuilder, StructureData

    # Filters that will match any elemental Silicon structure with 2 or less sites in total
    filters = {
        'attributes.sites': {
            'of_length': 2
        },
        'attributes.kinds': {
            'of_length': 1
        },
        'attributes.kinds.0.symbols.0': 'Si'
    }

    builder = QueryBuilder().append(StructureData, filters=filters)
    results = builder.first()  # pylint: disable=no-value-for-parameter

    if not results:
        alat = 5.43
        ase_structure = crystal(
            'Si',
            [(0, 0, 0)],
            spacegroup=227,
            cellpar=[alat, alat, alat, 90, 90, 90],
            primitive_cell=True,
        )
        structure = StructureData(ase=ase_structure)
        structure.store()
    else:
        structure = results[0]  # pylint: disable=unsubscriptable-object

    return structure.uuid


def get_cu_bulk_structure():
    """Return a `StructureData` representing bulk copper.

    The database will first be queried for the existence of a bulk Cu crystal. If this is not the case, one is
    created and stored. This function should be used as a default for CLI options that require a `StructureData` node.
    This way new users can launch the command without having to construct or import a structure first. This is the
    reason that we hardcode a bulk silicon crystal to be returned. More flexibility is not required for this purpose.

    :return: a `StructureData` representing bulk silicon
    """
    from ase.spacegroup import crystal
    from aiida.orm import QueryBuilder, StructureData

    # Filters that will match any elemental Copper structure with 1 sites in total
    filters = {
        'attributes.sites': {
            'of_length': 1
        },
        'attributes.kinds': {
            'of_length': 1
        },
        'attributes.kinds.0.symbols.0': 'Cu'
    }

    builder = QueryBuilder().append(StructureData, filters=filters)
    results = builder.first()  # pylint: disable=no-value-for-parameter

    if not results:
        alat = 3.615
        ase_structure = crystal(
            'Cu',
            [(0, 0, 0)],
            spacegroup=225,
            cellpar=[alat, alat, alat, 90, 90, 90],
            primitive_cell=True,
        )
        structure = StructureData(ase=ase_structure)
        structure.store()
    else:
        structure = results[0]  # pylint: disable=unsubscriptable-object

    return structure.uuid


# Parameters
def get_kkrpara_defaults():
    """
    Create or load kkrparams default values
    """
    from aiida.orm import QueryBuilder, Dict
    from aiida_kkr.tools import kkrparams

    # Filters that will match the Dict node with kkr_default_params
    filters = {'extras.aiida-kkr_cmdline_info': 'kkr_default_params'}

    # now query for this
    builder = QueryBuilder().append(Dict, filters=filters)
    results = builder.first()  # pylint: disable=no-value-for-parameter

    # create kkr params node or take from database if found
    if not results:
        defaults, version = kkrparams.get_KKRcalc_parameter_defaults()
        params = Dict(defaults)
        params.label = 'Default KKR parameters'
        params.description = f'Version of the KKR params: {version}'
        params.extras['aiida-kkr_cmdline_info'] = 'kkr_default_params'  # pylint: disable=no-member
        params.store()
    else:
        params = results[0]  # pylint: disable=unsubscriptable-object

    # return params node
    return params.uuid


# Codes
def get_voro():
    """Return a `Code` node of the latest added inpgen executable in the database."""
    return get_last_code('kkr.voro')


def get_kkr():
    """Return a `Code` node of the latest added inpgen executable in the database."""
    return get_last_code('kkr.kkr')


def get_kkrimp():
    """Return a `Code` node of the latest added inpgen executable in the database."""
    return get_last_code('kkr.kkrimp')


def get_last_code(entry_point_name):
    """Return a `Code` node of the latest code executable of the given entry_point_name in the database.

    The database will be queried for the existence of a inpgen node.
    If this is not exists and NotExistent error is raised.


    :param entry_point_name: string
    :return: the uuid of a inpgen `Code` node
    :raise: aiida.common.exceptions.NotExistent
    """
    from aiida.orm import QueryBuilder, Code
    from aiida.common.exceptions import NotExistent

    filters = {'attributes.input_plugin': {'==': entry_point_name}}

    builder = QueryBuilder().append(Code, filters=filters)
    builder.order_by({Code: {'ctime': 'desc'}})
    results = builder.first()  # pylint: disable=no-value-for-parameter

    if not results:
        raise NotExistent(f'ERROR: Could not find any Code in the database with entry point: {entry_point_name}!')
    else:
        code = results[0]  # pylint: disable=unsubscriptable-object

    return code.uuid
