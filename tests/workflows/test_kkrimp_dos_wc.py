#!/usr/bin/env python

import pytest
from ..dbsetup import *
from ..conftest import voronoi_local_code, kkrhost_local_code, data_dir, import_with_migration
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test


@pytest.mark.timeout(300, method='thread')
def test_dos_startpot_wc(clear_database_before_test, kkrimp_local_code, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow for impurity host-in-host
    """
    from numpy import array, loadtxt
    from masci_tools.io.kkr_params import kkrparams
    from aiida.orm import Code, load_node, Dict, StructureData, load_group
    from aiida.orm.querybuilder import QueryBuilder
    from aiida.manage.caching import enable_caching
    from aiida_kkr.tools import neworder_potential_wf
    from aiida_kkr.workflows.kkr_imp_dos import kkr_imp_dos_wc
    from aiida_kkr.calculations import KkrCalculation

    # import precomputed GF host writeout
    o = import_with_migration('data_dir/kkr_flex_wc-nodes-94b44d25a7af0e584343419207f3a7e1.tar.gz')
    print('import', o)
    # gf_workflow = load_node('0adfe62d-05aa-4c79-90ca-86b1753612a2')
    gf_workflow = [i for i in load_group(o).nodes if i.label == 'GF_writeout Cu bulk'][0]
    GF_host_calc = gf_workflow.get_outgoing(node_class=KkrCalculation).first().node
    print('GF_host_calc', GF_host_calc)

    # set workflow settings
    wfd = kkr_imp_dos_wc.get_wf_defaults()
    wfd['clean_impcalc_retrieved'] = False  # deactivate cleaning of unused data to regain cachability
    wfd['retrieve_kkrflex'] = True
    print(wfd)

    # set computer options
    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    options = Dict(options)

    # now create a SingleFileData node containing the impurity starting potential
    with GF_host_calc.outputs.retrieved.open('scoef') as _f:
        neworder_pot1 = [int(i) for i in loadtxt(_f, skiprows=1)[:, 3] - 1]
    settings_dict = {'pot1': 'out_potential', 'out_pot': 'potential_imp', 'neworder': neworder_pot1}
    settings = Dict(settings_dict)

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
    builder.wf_parameters = Dict(wfd)
    builder.impurity_info = Dict(imp_info)
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
