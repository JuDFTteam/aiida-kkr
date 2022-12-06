#!/usr/bin/env python

import pytest
from aiida.manage.tests.pytest_fixtures import aiida_profile
from ..dbsetup import queuename
from ..conftest import kkrimp_local_code, data_dir, import_with_migration
from aiida_testing.export_cache._fixtures import run_with_cache


# tests
def test_kkr_base_restart(aiida_profile, kkrimp_local_code, run_with_cache):
    """
    simple example that fixes the R_LOG error
    """
    from aiida.orm import Code, load_node, Dict, load_group, WorkChainNode, CalcJobNode
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.restart.restart_kkrimp import KkrimpCalculationBaseWorkChain

    # # first load parent voronoi calculation
    # group_pk = import_with_migration('data_dir/KkrCalculationBaseWorkChain-nodes-706d766ad8897a8e0a5dc48d3b07320d.tar.gz')
    # import_group = load_group(group_pk)
    # kkr_restart_wf = [i for i in import_group.nodes if isinstance(i, WorkChainNode)][0]
    # kkrhost_parent = [i for i in kkr_restart_wf.get_outgoing(node_class=CalcJobNode).all() if i.link_label=='iteration_02'][0]

    # # now create a SingleFileData node containing the impurity starting potential
    # from aiida_kkr.tools import neworder_potential_wf
    # settings_dict = {'neworder': [0]]}
    # settings = Dict(settings_dict)
    # startpot_imp_sfd = neworder_potential_wf(
    #     settings_node=settings, parent_calc_folder=kkrhost_parent.outputs.remote_folder
    # )

    # # set 1 simple mixing step
    # kkrimp_params = kkrparams(params_type='kkrimp')
    # kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
    # ParamsKKRimp = Dict(kkrimp_params.get_dict())

    # # create new KKRimp calculation
    # options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
    # builder = KkrimpCalculationBaseWorkChain.get_builder()
    # builder.calc.code = kkrimp_local_code
    # builder.calc.host_Greenfunction_folder = GF_host_calc.outputs.remote_folder
    # builder.calc.impurity_potential = startpot_imp_sfd
    # builder.calc.metadata.options = options
    # builder.calc.parameters = ParamsKKRimp

    # # run test
    # out, node = run_with_cache(builder, data_dir=data_dir)
    # # check if two voronoi runs were done (one failed, but was handled)
    # assert len(node.called) == 2
    # # check if cluster radius has right size
    # assert 'output_parameters' in out
    # out_dict = out['output_parameters'].get_dict()
    # assert out_dict['convergence_group']['rms'] < 0.18
