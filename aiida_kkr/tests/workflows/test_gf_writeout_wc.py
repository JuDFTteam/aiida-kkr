#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache
from ..conftest import voronoi_local_code, kkrhost_local_code
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile


@pytest.mark.timeout(240, method='thread')
def test_kkrflex_writeout_wc(aiida_profile, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData, RemoteData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.gf_writeout import kkr_flex_wc
    from numpy import array
    import os

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd =kkr_flex_wc.get_wf_defaults()
    options = {'queue_name' : queuename, 'resources': {"num_machines": 1},'max_wallclock_seconds' : 5*60, 'custom_scheduler_commands' : '', 'withmpi' : False}
    options = Dict(dict=options)

    imp_info = Dict(dict={'Rcut':2.5533, 'ilayer_center': 0, 'Zimp':[29.]})

    label = 'GF_writeout Cu bulk'
    descr = 'GF_writeout workflow for Cu bulk'

    from aiida.tools.importexport import import_data
    import_data('files/db_dump_kkrcalc.tar.gz')
    kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

    # create process builder to set parameters
    builder = kkr_flex_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkr = kkrhost_local_code
    builder.options = options
    builder.remote_data = kkr_calc_remote
    builder.impurity_info = imp_info

    # now run calculation
    out, node = run_with_cache(builder)
    print(out)

    n = out['workflow_info']
    n = n.get_dict()

    assert n.get('successful')
    assert n.get('list_of_errors') == []

    d = out['GF_host_remote']
    assert isinstance(d, RemoteData)

    kkrflex_calc = load_node(n.get('pk_flexcalc'))
    kkrflex_retrieved = kkrflex_calc.outputs.retrieved
    for name in 'tmat green atominfo intercell_cmoms intercell_ref'.split():
        assert 'kkrflex_'+name in kkrflex_retrieved.list_object_names()


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   test_kkrflex_writeout_wc()
