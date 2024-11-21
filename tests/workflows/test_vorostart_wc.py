#!/usr/bin/env python

import pytest
from aiida.engine import run_get_node
from ..dbsetup import *
from ..conftest import voronoi_local_code, kkrhost_local_code, test_dir, data_dir, import_with_migration


@pytest.mark.timeout(500, method='thread')
def test_kkr_startpot_wc_Cu(clear_database_before_test, voronoi_local_code, kkrhost_local_code, enable_archive_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.voro_start import kkr_startpot_wc
    from numpy import array

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
    wfd = kkr_startpot_wc.get_wf_defaults()
    wfd['check_dos'] = True
    wfd['natom_in_cls_min'] = 20
    wfd['num_rerun'] = 2
    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }
    params_vorostart = Dict(wfd)

    # Finally we use the kkrparams class to prepare a valid set of KKR parameters that are stored as a Dict object for the use in aiida
    ParaNode = Dict(kkrparams(LMAX=2, NSPIN=1, RCLUSTZ=1.9).get_dict())

    # create process builder to set parameters
    builder = kkr_startpot_wc.get_builder()
    builder.calc_parameters = ParaNode
    builder.metadata.description = 'voronoi startpot workflow for Cu bulk'
    builder.metadata.label = 'startpot for Cu bulk'
    builder.voronoi = voronoi_local_code
    builder.structure = Cu
    builder.wf_parameters = params_vorostart
    builder.options = Dict(options)
    builder.kkr = kkrhost_local_code

    with enable_archive_cache(data_dir / 'kkr_startpot_wc_Cu.aiida'):
        out, node = run_get_node(builder)
    print('outputs:', node, out)

    # check output
    n = out['results_vorostart_wc']
    n = n.get_dict()
    print('results_dict:', n)
    assert n.get('successful')
    assert n.get('last_voro_ok')
    assert n.get('list_of_errors') == []
    assert abs(n.get('starting_fermi_energy') - 0.409241) < 10**-14

    kkrcalc = [i for i in node.called_descendants if i.process_label == 'KkrCalculation'][0]
    vorocalc = [i for i in node.called_descendants if i.process_label == 'VoronoiCalculation'][0]
    print('hashes of computed voro and kkr calcs', vorocalc.get_hash(), kkrcalc.get_hash())
    print('was cached (voro/kkr)?', vorocalc.get_cache_source(), kkrcalc.get_cache_source())
    print('caching info voro', vorocalc._get_objects_to_hash())
    print('caching info kkr ', kkrcalc._get_objects_to_hash())

    # find imported and new calculations
    from aiida_kkr.calculations import VoronoiCalculation, KkrCalculation
    from aiida.orm import QueryBuilder

    qub = QueryBuilder()
    qub.append(VoronoiCalculation)
    voro_calcs = qub.all()
    print('\n\nvoronoi calculations:', voro_calcs)
    for calc in voro_calcs:
        code = calc[0].inputs.code
        print('hash', calc[0].get_hash())
        print('was cached?', calc[0].get_cache_source())
        print('caching info', calc[0]._get_objects_to_hash())
        print('code info:', code.get_hash())
        print('code info:', code._get_objects_to_hash())
        print('code info:', code._hash_ignored_attributes)
        print('code attrs:', code.attributes)

    qub = QueryBuilder()
    qub.append(KkrCalculation)
    kkr_calcs = qub.all()
    print('\n\nkkr calculations:', kkr_calcs)


@pytest.mark.timeout(500, method='thread')
def test_kkr_startpot_parent_KKR(clear_database_before_test, voronoi_local_code, enable_archive_cache):
    """
    simple Cu noSOC, FP, lmax2 full example using scf workflow
    """
    from aiida.orm import Code, load_node, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.voro_start import kkr_startpot_wc
    from numpy import array

    # here we create a parameter node for the workflow input (workflow specific parameter) and adjust the convergence criterion.
    wfd = kkr_startpot_wc.get_wf_defaults()
    wfd['check_dos'] = False
    wfd['natom_in_cls_min'] = 20
    wfd['num_rerun'] = 2
    options = {
        'queue_name': queuename,
        'resources': {
            'num_machines': 1
        },
        'max_wallclock_seconds': 5 * 60,
        'withmpi': False,
        'custom_scheduler_commands': ''
    }

    # load necessary files from db_dump files
    import_with_migration('files/db_dump_kkrcalc.tar.gz')

    # first load parent voronoi calculation
    kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

    # extract KKR parameter and remote_data folder
    params_kkr_parent = kkr_calc.inputs.parameters
    parent_calc_remote = kkr_calc.outputs.remote_folder

    # increase lmax value
    params = kkrparams(params_type='voronoi', **params_kkr_parent)
    params.set_multiple_values(LMAX=3)

    # create process builder to set parameters
    builder = kkr_startpot_wc.get_builder()
    builder.calc_parameters = Dict(params.get_dict())
    builder.metadata.description = 'voronoi startpot workflow with parent_KKR input'
    builder.metadata.label = 'startpot for increased lmax'
    builder.voronoi = voronoi_local_code
    builder.wf_parameters = Dict(wfd)
    builder.options = Dict(options)
    builder.parent_KKR = parent_calc_remote

    with enable_archive_cache(data_dir / 'kkr_startpot_parent_KKR.aiida'):
        out, node = run_get_node(builder)
    print('outputs:', node, out)

    # check output
    n = out['results_vorostart_wc']
    n = n.get_dict()
    print('results_dict:', n)
    assert n.get('successful', False)


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    #Test = Test_vorostart_workflow()
    #Test.test_vorostart_wc_Cu()
    test_vorostart_wc_Cu()
