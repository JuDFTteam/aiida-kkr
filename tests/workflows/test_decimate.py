#!/usr/bin/env python

import numpy as np
from aiida.orm import load_node, Dict, load_group, KpointsData
from aiida.engine import run_get_node
from aiida_kkr.workflows import kkr_decimation_wc
from ..conftest import import_with_migration
from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir, import_with_migration


def get_builder(dosmode, kkrhost_local_code, voronoi_local_code):
    # import slab structure
    import_with_migration(test_dir / 'files/deci_data.aiida')
    scf_Au9 = load_node('f7012a27-98f8-4c5c-ba29-b1ab2109c427')

    # set computer options
    options = {
        'queue_name': '',
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }

    # create process builder
    builder = kkr_decimation_wc.get_builder()
    builder.kkr = kkrhost_local_code
    builder.voronoi = voronoi_local_code
    builder.options = Dict(options)
    builder.remote_data = scf_Au9.outputs.last_RemoteData
    # use DOS mode for faster test
    settings = kkr_decimation_wc.get_wf_defaults()
    settings['nkz'] = 5
    settings['nprinc'] = 3
    settings['nplayer'] = 2
    settings['dosmode'] = dosmode
    settings['dos_params'] = {'emin_EF': -5.0, 'emax_EF': 3.0, 'nepts': 1, 'tempr': 200, 'kmesh': [5, 5, 5]}
    builder.wf_parameters = Dict(dict=settings)
    # decrease cluster radius for faster results
    builder.calc_parameters = Dict(dict(RCLUSTZ=2.0))

    return builder


def test_decimate(clear_database_before_test, kkrhost_local_code, voronoi_local_code, enable_archive_cache):
    """
    test for decimation workflow
    """

    # run test
    with enable_archive_cache(data_dir / 'decimate.aiida'):
        out, node = run_get_node(get_builder(False, kkrhost_local_code, voronoi_local_code))

    print((out, node))

    # check outcome
    for key in [
        'structure_decimate',
        'structure_substrate',
        'out_params_calc_deci_out',
        'out_params_calc_decimate',
        'out_remote_calc_decimate',
        'out_retrieved_calc_decimate',
    ]:
        assert key in out

    print(out['out_retrieved_calc_decimate'].list_object_names())
    from pprint import pprint
    with out['out_retrieved_calc_decimate'].open('inputcard') as _f:
        pprint(_f.readlines())
    with out['out_retrieved_calc_decimate'].open('out_kkr') as _f:
        pprint(_f.readlines())
    assert min(out['out_params_calc_decimate']['convergence_group']['charge_neutrality_all_iterations']) < 0.01


def test_decimate_dos(
    clear_database_before_test, kkrhost_local_code, voronoi_local_code, enable_archive_cache, ndarrays_regression
):
    """
    test for decimation workflow
    """

    # run test
    with enable_archive_cache(data_dir / 'decimate_dos.aiida'):
        out, node = run_get_node(get_builder(True, kkrhost_local_code, voronoi_local_code))

    print((out, node))

    # check outcome
    assert 'lmdos.06.2.dat' in out['out_retrieved_calc_decimate'].list_object_names()
    with out['out_retrieved_calc_decimate'].open('lmdos.06.2.dat') as _f:
        dos = np.loadtxt(_f)
    ndarrays_regression.check({'dos': dos})


def test_decimate_bandstruc(
    clear_database_before_test, kkrhost_local_code, voronoi_local_code, enable_archive_cache, ndarrays_regression
):
    """
    test for decimation workflow
    """

    # run test
    builder = get_builder(False, kkrhost_local_code, voronoi_local_code)
    kpts = KpointsData()
    kpts.set_cell([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    kpts.set_kpoints([[0., 0., 0.], [0., 1., 0.]], cartesian=True)
    builder.kpoints = kpts
    with enable_archive_cache(data_dir / 'decimate_bandstruc.aiida'):
        out, node = run_get_node(builder)

    print((out, node))

    # check outcome
    print(out['out_retrieved_calc_decimate'].list_object_names())
    from pprint import pprint
    with out['out_retrieved_calc_decimate'].open('inputcard') as _f:
        pprint(_f.readlines())
    with out['out_retrieved_calc_decimate'].open('out_kkr') as _f:
        pprint(_f.readlines())

    # check outcome
    assert 'qdos.06.2.dat' in out['out_retrieved_calc_decimate'].list_object_names()
    with out['out_retrieved_calc_decimate'].open('qdos.06.2.dat') as _f:
        qdos = np.loadtxt(_f)
    ndarrays_regression.check({'qdos': qdos})
