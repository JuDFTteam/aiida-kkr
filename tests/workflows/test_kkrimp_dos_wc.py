#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from ..dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
from aiida_kkr.conftest import voronoi_local_code, kkrhost_local_code, data_dir
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test
from ..conftest import import_with_migration


@pytest.mark.timeout(300, method='thread')
def test_dos_startpot_wc(clear_database_before_test, kkrimp_local_code, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.kkr_imp_dos import kkr_imp_dos_wc
    from numpy import array

    # import precomputed GF host writeout
    import_with_migration('files/db_dump_kkrflex_create.tar.gz')
    GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

    wfd = kkr_imp_dos_wc.get_wf_defaults()
    wfd['clean_impcalc_retrieved'] = False  # deactivate cleaning of unused data to regain cachability
    print(wfd)

    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    options = Dict(dict=options)

    # now create a SingleFileData node containing the impurity starting potential
    from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
    from numpy import loadtxt
    neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:, 3] - 1]
    settings_dict = {'pot1': 'out_potential', 'out_pot': 'potential_imp', 'neworder': neworder_pot1}
    settings = Dict(dict=settings_dict)

    from aiida.manage.caching import enable_caching
    with enable_caching():  # should enable caching globally in this python interpreter
        startpot_imp_sfd = neworder_potential_wf(
            settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder
        )

    label = 'kkrimp_dos Cu host_in_host'
    descr = 'kkrimp_dos workflow for Cu bulk'

    imp_info = GF_host_calc.inputs.impurity_info.get_dict()
    imp_info['Rcut'] = 2.5533
    print(imp_info)

    # create process builder to set parameters
    builder = kkr_imp_dos_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.options = options
    builder.kkr = kkrhost_local_code
    builder.kkrimp = kkrimp_local_code
    builder.imp_pot_sfd = startpot_imp_sfd
    builder.wf_parameters = Dict(dict=wfd)
    builder.impurity_info = Dict(dict=imp_info)
    builder.host_remote = GF_host_calc.outputs.remote_folder

    # now run calculation
    out, node = run_with_cache(builder, data_dir=data_dir)
    print(node)
    print(out)
    from aiida.orm import WorkChainNode
    for i in node.get_outgoing(node_class=WorkChainNode).all():
        print(i.node, list(i.node.outputs))

    assert 'last_calc_info' in list(out.keys())
    assert 'last_calc_output_parameters' in list(out.keys())
    assert 'workflow_info' in list(out.keys())
    assert 'dos_data' in list(out.keys())
    assert 'dos_data_interpol' in list(out.keys())
    assert len(out['dos_data_interpol'].get_y()) == 5
    assert len(out['dos_data_interpol'].get_y()[0]) == 3
    assert len(out['dos_data_interpol'].get_y()[0][0]) == 20


"""
@pytest.mark.timeout(300, method='thread')
def test_dos_reuse_gf_writeout(self):
    pass


@pytest.mark.timeout(300, method='thread')
def test_dos_from_kkrimp_sub(self):
    pass


@pytest.mark.timeout(300, method='thread')
def test_dos_from_kkrimp_full(self):
    pass
"""

#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    test_dos_startpot_wc()
    #test_dos_reuse_gf_writeout()
    #test_dos_from_kkrimp_sub()
    #test_dos_from_kkrimp_full()
