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
__version__ = '0.1.0'
__contributors__ = (u'Philipp Rüßmann')

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
    imp_info2 = orm.Dict({'ilayer_center': add_position['ilayer'], 'Zimp': [Zadd], 'Rcut': 1e-5})
    clust2 = get_scoef_single_imp(host_structure, imp_info2)
    return imp_info2, clust2


def get_r_offset(clust1, clust2, host_structure, add_position):
    """
    find offset vector in rcls units
    """
    # calculate out-of plane vector from the ilayer indices of the two clusters
    r_out_of_plane = np.array([0., 0., 0.])
    ilayer1 = int(clust1[0, 3])
    ilayer2 = int(clust2[0, 3])
    if ilayer1 != ilayer2:
        pos1 = np.array(host_structure.sites[ilayer1].position)
        pos2 = np.array(host_structure.sites[ilayer2].position)
        r_out_of_plane = pos2 - pos1

    # calculate in-plane vector from da, db inputs
    da = add_position.get_dict().get('da', 0)
    db = add_position.get_dict().get('db', 0)
    cell = np.array(host_structure.cell)
    r_in_plane = da * cell[0] + db * cell[1]

    # combine to offset vector
    r_offset = r_out_of_plane + r_in_plane

    # convert from Ang to alat units (used internally in KKR)
    alat = get_alat_from_bravais(np.array(host_structure.cell), host_structure.pbc[2])
    r_offset /= alat

    return r_offset


def offset_clust2(clust1, clust2, host_structure, add_position):
    """
    Compute and add offset to clust2
    """
    r_offset = get_r_offset(clust1, clust2, host_structure, add_position)

    clust2_offset = clust2.copy()
    clust2_offset[:, :3] += r_offset

    return clust2_offset


def get_imp_info_add_position(host_calc, imp_info, add_position):
    """
    Create combined impurity info node for the original
    imp cluster + an additional (STM tip) position
    """

    # extract host structure
    host_structure = find_parent_structure(host_calc)

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
def get_imp_info_add_position_cf(host_remote, imp_info, add_position):
    """
    Create a new impurty info node that combines the impurity cluster
    of an original calculation and an STM scanning position.
    """

    # then combine the imp info
    imp_info_combined = get_imp_info_add_position(host_remote, imp_info, add_position)

    return imp_info_combined


##############################################################################
# combine potentials


def extract_host_potential(add_position, host_remote):
    """
    Extract the potential of the position in the host that matches the additional position
    """

    # find ilayer from input node
    ilayer = add_position['ilayer']

    # get host calculation from remote
    host_calc = host_remote.get_incoming(node_class=orm.CalcJobNode).first().node

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


def add_host_potential_to_imp(add_position, host_remote, imp_potential_node):
    """
    combine host potential with impurity potential
    """
    # get add potential from host
    pot_add = extract_host_potential(add_position, host_remote)

    # get impurity potential and convert to list
    pot_imp = imp_potential_node.get_content().split('\n')
    pot_imp = [line + '\n' for line in pot_imp if line != '']

    # glue potentials together and create SinglefileData
    pot_combined = pot_imp + pot_add

    return pot_combined


def create_combined_potential_node(add_position, host_remote, imp_potential_node):
    """
    Combine impurity potential with an additional potential from the host for
    the STM tip position (additional position)
    """
    import io

    # combine potential texts
    pot_combined = add_host_potential_to_imp(add_position, host_remote, imp_potential_node)

    # convert to byte string and put into SinglefilData node
    pot_combined_string = ''
    for line in pot_combined:
        pot_combined_string += line
    pot_combined_node = orm.SinglefileData(io.BytesIO(bytes(pot_combined_string, 'utf8')))

    return pot_combined_node


@engine.calcfunction
def create_combined_potential_node_cf(add_position, host_remote, imp_potential_node):
    """
    Calcfunction that combines the impurity potential with an addition potential site from the host
    """

    pot_combined_node = create_combined_potential_node(add_position, host_remote, imp_potential_node)

    return pot_combined_node
