#!/usr/bin/env python

import pytest
from ..dbsetup import *
from ..conftest import voronoi_local_code, kkrhost_local_code, kkrimp_local_code, data_dir, import_with_migration
from aiida.engine import run_get_node
from aiida import orm


@pytest.mark.timeout(900, method='thread')
def test_imp_BdG_wc(
    clear_database_before_test, kkrimp_local_code, kkrhost_local_code, voronoi_local_code, enable_archive_cache,
    ndarrays_regression
):
    """
    simple Cu noSOC host-in-host, use imp_BdG workchain
    """
    from aiida_kkr.workflows import kkrimp_BdG_wc

    # import host calculation with SOC solver
    o = import_with_migration(data_dir / 'jij_soc.aiida')
    host_parent_calc = [i for i in orm.load_group(o).nodes if i.label == 'jij_calc_z'][0]

    # set up process builder
    builder = kkrimp_BdG_wc.get_builder()

    # codes and computer options
    builder.kkr = kkrhost_local_code
    builder.kkrimp = kkrimp_local_code
    builder.voronoi = voronoi_local_code

    options = orm.Dict({
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 12 * 60,
        'withmpi': False,
        'custom_scheduler_commands': 'ulimit -s unlimited'
    })
    builder.options = options

    # impurity setup
    imp_info = orm.Dict({'Rcut': 3.5, 'ilayer_center': 0, 'Zimp': [47.]})
    builder.impurity_info = imp_info

    # setting for impurity scf
    builder.imp_scf.remote_data_host = host_parent_calc.outputs.remote_folder
    builder.imp_scf.wf_parameters = orm.Dict({
        'kkr_runmax': 1,
        'nsteps': 1,
        'retrieve_kkrflex': True,
        'do_final_cleanup': False,
        'accuracy_params':
        {  # lower accuracy to have faster runtime
            'RADIUS_LOGPANELS': 0.6,
            'NPAN_LOG': 3,
            'NPAN_EQ': 7,
            'NCHEB': 6
        },
    })
    builder.imp_scf.scf.params_overwrite = orm.Dict({'TOL_ALAT_CHECK': 1e-8})
    builder.imp_scf.gf_writeout.params_kkr_overwrite = orm.Dict({'RUNOPT': []})

    # activate BdG DOS step
    builder.calc_DOS = orm.Bool(True)
    builder.dos.wf_parameters = orm.Dict({
        'retrieve_kkrflex': True,
        'dos_params': {
            'nepts': 3,
            'tempr': 200.0,
            'emin': -1.0,
            'emax': 1.0,
            'kmesh': [10, 10, 10]
        },
        'clean_impcalc_retrieved': False
    })
    # Activate BdG for DOS gf writeout step
    builder.dos.gf_writeout.params_kkr_overwrite = orm.Dict({
        'RUNOPT': [],
        '<USE_BDG>': True,
        '<DELTA_BDG>': 1e-2,
        'NPAN_LOG': 3,
        'NPAN_EQ': 7,
        'NCHEB': 6,
    })

    # now run calculation
    with enable_archive_cache(data_dir / 'imp_BdG_wc.aiida'):
        out, node = run_get_node(builder)
    print(node)
    print(out)

    check_dict = {'x': node.outputs.dos_data.get_x()[1][0], 'y': node.outputs.dos_data.get_y()[0][1]}
    ndarrays_regression.check(check_dict)
