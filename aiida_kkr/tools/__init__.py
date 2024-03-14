"""
tools provided by aiida-kkr plugin
"""

# import all tools here to expose them in `aiida_kkr.tools` directly
from .common_workfunctions import (
    update_params_wf, prepare_VCA_structure_wf, prepare_2Dcalc_wf, test_and_get_codenode, get_inputs_kkr,
    get_inputs_kkrimporter, get_inputs_voronoi, get_inputs_kkrimp, get_parent_paranode,
    generate_inputcard_from_structure, check_2Dinput_consistency, structure_from_params, vca_check
)
from .find_cluster_radius import find_cluster_radius
from .plot_kkr import plot_kkr
from .parse_dos import parse_dosfiles
from .tools_kkrimp import modify_potential, rotate_onto_z, find_neighbors, make_scoef
# make the most important things from masci-tools importable here
from masci_tools.io.kkr_params import kkrparams
from .multi_imps_data_extract import MultiImpuritiesData
from .kick_out_core_states import *
from .neworder_potential import *
from .find_parent import get_calc_from_remote, get_remote, get_parent
from .bdg_tools import get_anomalous_density_data
from .ldau import *


# expose structure finder from VoronoiCalculation
def find_parent_structure(calc_or_remote, return_voro=False):
    """
    Traverse the provenance graph upwards to find the input structure
    from a KkrCalculation or VoronoiCalculation node or their remote_folder
    nodes, respectively.

    :param calc_or_remote: CalcJobNode or RemoteData node of VoronoiCalculation or KkrCalculation
    :return struc: parent StructureData node
    """
    from .find_parent import find_parent_structure

    struc, voro_calc = find_parent_structure(calc_or_remote)

    if return_voro:
        return struc, voro_calc
    else:
        return struc


def search_kkrparams(search_string):
    """
    Search keywords and their description of kkrparams. Useful to find the correct name of a parameter.

    :param search_string: string which is searched (case insensitive) in parameter names and descriptions of kkrparams
    """
    kkrparams().get_description(search=search_string)
