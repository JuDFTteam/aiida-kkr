"""
tools provided by aiida-kkr plugin
"""

# import all tools here to expose them in `aiida_kkr.tools` directly
from .common_workfunctions import (
    update_params_wf, prepare_VCA_structure_wf, prepare_2Dcalc_wf, test_and_get_codenode, get_inputs_kkr,
    get_inputs_kkrimporter, get_inputs_voronoi, get_inputs_kkrimp, get_parent_paranode,
    generate_inputcard_from_structure, check_2Dinput_consistency, structure_from_params, neworder_potential_wf,
    vca_check, kick_out_corestates_wf, find_cluster_radius
)
from .plot_kkr import plot_kkr
from .tools_kkrimp import modify_potential, rotate_onto_z, find_neighbors, make_scoef
# make the most important things from masci-tools importable here
from masci_tools.io.kkr_params import kkrparams


# expose structure finder from VoronoiCalculation
def find_parent_structure(calc_or_remote, return_voro=False):
    """
    Traverse the provenance graph upwards to find the input structure
    from a KkrCalculation or VoronoiCalculation node or their remote_folder
    nodes, respectively.

    :param calc_or_remote: CalcJobNode or RemoteData node of VoronoiCalculation or KkrCalculation
    :return struc: parent StructureData node
    """
    from aiida_kkr.calculations.voro import VoronoiCalculation

    struc, voro_calc = VoronoiCalculation.find_parent_structure(calc_or_remote)

    return struc
