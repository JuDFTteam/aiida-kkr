"""
Helper tools to find parent calculations or structures.
"""

from aiida.common.exceptions import NotExistent


def get_struc(parent_calc):
    """
    Get structure from a parent_folder (result of a calculation, typically a remote folder)
    """
    return parent_calc.inputs.structure


def has_struc(parent_folder):
    """
    Check if parent_folder has structure information in its input
    """
    return 'structure' in parent_folder.base.links.get_incoming().all_link_labels()


def get_remote(parent_folder):
    """
    get remote_folder from input if parent_folder is not already a remote folder
    """
    parent_folder_tmp0 = parent_folder
    try:
        parent_folder_tmp = parent_folder_tmp0.base.links.get_incoming().get_node_by_label('remote_folder')
    except NotExistent:
        try:
            # check if GFhost_folder is there, this is the case for a KkrimpCalculation
            parent_folder_tmp = parent_folder_tmp0.base.links.get_incoming().get_node_by_label('GFhost_folder')
        except NotExistent:
            parent_folder_tmp = parent_folder_tmp0
    return parent_folder_tmp


def get_parent(input_folder):
    """
    Get the parent folder of the calculation. If no parent was found return input folder
    """
    input_folder_tmp0 = input_folder

    # first option: parent_calc_folder (KkrimpCalculation)
    try:
        parent_folder_tmp = input_folder_tmp0.base.links.get_incoming().get_node_by_label('parent_calc_folder')
        return_input = False
    except NotExistent:
        return_input = True

    # second option: parent_folder (KkrCalculation)
    try:
        parent_folder_tmp = input_folder_tmp0.base.links.get_incoming().get_node_by_label('parent_folder')
        return_input = False
    except NotExistent:
        return_input = return_input & True

    # third option: parent_KKR option (special mode of VoronoiCalculation)
    try:
        parent_folder_tmp = input_folder_tmp0.base.links.get_incoming().get_node_by_label('parent_KKR')
        return_input = False
    except NotExistent:
        return_input = return_input & True

    if return_input:
        parent_folder_tmp = input_folder_tmp0

    return parent_folder_tmp


def find_parent_structure(parent_folder):
    """
    Find the Structure node recuresively in chain of parent calculations (structure node is input to voronoi calculation)
    """
    iiter = 0
    Nmaxiter = 1000
    parent_folder_tmp = get_remote(parent_folder)
    while not has_struc(parent_folder_tmp) and iiter < Nmaxiter:
        parent_folder_tmp = get_remote(get_parent(parent_folder_tmp))
        iiter += 1
        if iiter % 200 == 0:
            print(
                'Warning: find_parent_structure takes quite long (already searched {} ancestors). Stop after {}'.format(
                    iiter, Nmaxiter
                )
            )
    if has_struc(parent_folder_tmp):
        struc = get_struc(parent_folder_tmp)
        return struc, parent_folder_tmp
    else:
        raise ValueError('structure not found')


def get_calc_from_remote(calc_remote):
    """
    Get the parent calculation from a RemoteData
    """
    from aiida.orm import CalcJobNode
    from aiida.orm import RemoteData

    if not isinstance(calc_remote, RemoteData):
        raise ValueError('input node is not a RemoteData folder')

    parents = calc_remote.base.links.get_incoming(node_class=CalcJobNode).all()
    if len(parents) != 1:
        raise ValueError('Parent is not unique!')

    return parents[0].node
