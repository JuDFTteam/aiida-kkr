#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
import pytest
from ..dbsetup import *
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
from ..conftest import voronoi_local_code, kkrhost_local_code
from ..conftest import import_with_migration, data_dir
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test


@pytest.mark.timeout(240, method='thread')
def test_dos_wc_Cu(clear_database_before_test, kkrhost_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida import get_version
    from aiida.orm import Code, load_node, Dict, StructureData
    from aiida.orm import Computer
    from aiida.orm.querybuilder import QueryBuilder
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.calculations import KkrCalculation
    from aiida_kkr.workflows.dos import kkr_dos_wc
    from numpy import array

    print(f'AiiDA version: {get_version()}')

    # Then set up the structure
    alat = 6.83  # in a_Bohr
    abohr = 0.52917721067  # conversion factor to Angstroem units
    bravais = array([[0.5, 0.5, 0.0], [0.5, 0.0, 0.5], [0.0, 0.5, 0.5]])  # bravais vectors
    a = 0.5 * alat * abohr
    Cu = StructureData(cell=[[a, a, 0.0], [a, 0.0, a], [0.0, a, a]])
    Cu.append_atom(position=[0.0, 0.0, 0.0], symbols='Cu')

    Cu.store()
    print(Cu)

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd = kkr_dos_wc.get_wf_defaults()
    wfd['kmesh'] = [10, 10, 10]
    wfd['nepts'] = 10
    params_dos = Dict(dict=wfd)

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

    label = 'dos Cu bulk'
    descr = 'DOS workflow for Cu bulk'

    import_with_migration('files/db_dump_kkrcalc.tar.gz')
    kkr_calc_remote = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566').outputs.remote_folder

    # create process builder to set parameters
    builder = kkr_dos_wc.get_builder()
    builder.metadata.description = descr
    builder.metadata.label = label
    builder.kkr = kkrhost_local_code
    builder.wf_parameters = params_dos
    builder.options = options
    builder.remote_data = kkr_calc_remote

    # now run calculation
    out, node = run_with_cache(builder, data_dir=data_dir)

    # check outcome
    print('check outputs', node.exit_status, out)
    print(node.get_outgoing(node_class=KkrCalculation).first().node.outputs.retrieved.list_object_names())
    n = out['results_wf']
    n = n.get_dict()
    assert n.get('successful')
    assert n.get('list_of_errors') == []
    print(n)
    assert n.get('dos_params').get('nepts') == 10

    d = out['dos_data']
    x = d.get_x()
    y = d.get_y()

    #print(x[1][0])
    #print(y[0][1][0])
    assert sum(
        abs(
            x[1][0] - array([
                -1.00000441e+01, -8.33337744e+00, -6.66671074e+00, -5.00004407e+00, -3.33337737e+00, -1.66671073e+00,
                -4.40824453e-05, 1.66662256e+00, 3.33328921e+00, 4.99995599e+00
            ])
        )
    ) < 1e-6
    assert sum(
        abs(
            y[0][1][0] - array([
                0.00509658, 0.09912667, 0.33763353, 0.44783853, 3.48078323, 3.31767819, 0.7819404, 0.12646486,
                0.24972391, 0.19468126
            ])
        )
    ) < 1e-6


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    test_dos_wc_Cu()
