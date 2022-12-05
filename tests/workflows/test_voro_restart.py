#!/usr/bin/env python

import pytest
from aiida.manage.tests.pytest_fixtures import aiida_profile
from ..dbsetup import queuename
from ..conftest import voronoi_local_code, data_dir
from aiida_testing.export_cache._fixtures import run_with_cache


# tests
def test_voronoi_base_restart(aiida_profile, voronoi_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example
    """
    from aiida.orm import Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.restart.restart_voro import VoronoiCalculationBaseWorkChain

    # create a sample StructureData
    # with an additional Cu atom to artificially increase the number of atoms in the cluster
    alat = 3.61  # lattice constant in Angstroem
    bravais = [[0.5 * alat, 0.5 * alat, 0], [0.5 * alat, 0, 0.5 * alat], [0, 0.5 * alat,
                                                                          0.5 * alat]]  # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0, 0, 0], symbols='Cu')
    Cu.append_atom(position=[0, 0.1, 0], symbols='Cu')

    # create Dict input node using kkrparams class from masci-tools
    # use large RLCUSTZ which is then reduced in the BaseWorkChain workchain
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=4.5)
    ParaNode = Dict(params.get_dict())

    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
    builder = VoronoiCalculationBaseWorkChain.get_builder()
    builder.calc.code = voronoi_local_code
    builder.calc.metadata.options = options
    builder.calc.parameters = ParaNode
    builder.calc.structure = Cu

    # run test
    out, node = run_with_cache(builder, data_dir=data_dir)
    print(out, node)
    # check if two voronoi runs were done (one failed, but was handled)
    assert len(node.called) == 2
    # check if cluster radius has right size
    assert 'output_parameters' in out
    out_dict = out['output_parameters'].get_dict()
    assert out_dict['cluster_info_group']['cluster_info_atoms'][1]['sites'] == 530
