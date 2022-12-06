#!/usr/bin/env python

import pytest
from aiida.manage.tests.pytest_fixtures import aiida_profile
from ..dbsetup import queuename
from ..conftest import kkrhost_local_code, data_dir, import_with_migration
from aiida_testing.export_cache._fixtures import run_with_cache


# tests
def test_kkr_base_restart(aiida_profile, kkrhost_local_code, run_with_cache):
    """
    simple example that fixes the R_LOG error
    """
    from aiida.orm import load_node, Dict
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.restart.restart_kkr import KkrCalculationBaseWorkChain

    # load necessary files from db_dump files
    import_with_migration('files/db_dump_vorocalc.tar.gz')

    # first load parent voronoi calculation
    voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

    # extract and update KKR parameter (add missing values)
    params = kkrparams(**voro_calc.inputs.parameters.get_dict())
    params.set_multiple_values(RMAX=7., GMAX=65., DECOUPLE_SPIN_CHEBY=True, USE_CHEBYCHEV_SOLVER=True, R_LOG=0.01)
    params_node = Dict(params.get_dict())

    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
    builder = KkrCalculationBaseWorkChain.get_builder()
    builder.calc.code = kkrhost_local_code
    builder.calc.metadata.options = options
    builder.calc.parameters = params_node
    builder.calc.parent_folder = voro_calc.outputs.remote_folder

    # run test
    out, node = run_with_cache(builder, data_dir=data_dir)
    # check if two voronoi runs were done (one failed, but was handled)
    assert len(node.called) == 2
    # check if cluster radius has right size
    assert 'output_parameters' in out
    out_dict = out['output_parameters'].get_dict()
    assert out_dict['convergence_group']['rms'] < 0.18
