#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
import pytest
import pathlib
from aiida.manage.tests.pytest_fixtures import clear_database, clear_database_after_test, clear_database_before_test
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint, absolute_archive_path
from ..dbsetup import *
from ..conftest import voronoi_local_code, test_dir, data_dir, import_with_migration

kkr_codename = 'kkrhost'

#TODO
# implement missing tests:
# * test_vca_structure
# * test_overwrite_alat_input
# * test_voronoi_after_kkr
# * test_overwrite_potential


# tests
def test_voronoi_dry_run(aiida_profile, voronoi_local_code):
    """
    simple Cu noSOC, FP, lmax2 full example
    """
    from aiida.orm import Code, Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # create StructureData instance for Cu
    alat = 3.61  # lattice constant in Angstroem
    bravais = [[0.5 * alat, 0.5 * alat, 0], [0.5 * alat, 0, 0.5 * alat], [0, 0.5 * alat,
                                                                          0.5 * alat]]  # Bravais matrix in Ang. units
    Cu = StructureData(cell=bravais)
    Cu.append_atom(position=[0, 0, 0], symbols='Cu')

    # create Dict input node using kkrparams class from masci-tools
    params = kkrparams(params_type='voronoi')
    params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    ParaNode = Dict(dict=params.get_dict())

    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = ParaNode
    builder.structure = Cu
    builder.metadata.dry_run = True
    from aiida.engine import run
    run(builder)


def test_voronoi_cached(clear_database_before_test, voronoi_local_code, run_with_cache):
    """
    simple Cu noSOC, FP, lmax2 full example
    """
    import numpy as np
    from masci_tools.io.kkr_params import kkrparams
    from aiida.orm import Code, Dict, StructureData
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # create StructureData instance for Cu
    alat = 3.61  # lattice constant in Angstroem
    bravais = [[0.5 * alat, 0.5 * alat, 0.0], [0.5 * alat, 0.0, 0.5 * alat],
               [0.0, 0.5 * alat, 0.5 * alat]]  # Bravais matrix in Ang. units
    structure = StructureData(cell=np.round(bravais, 3))
    structure.append_atom(position=[0, 0, 0], symbols='Cu')

    # create Dict input node using kkrparams class from masci-tools
    kkr_params = kkrparams(params_type='voronoi')
    kkr_params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
    parameters = Dict(dict={k: v for k, v in kkr_params.items() if v})

    # computer options
    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}

    # set up builder
    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = parameters
    builder.structure = structure
    # now run calculation or use cached result
    print('data_dir:', data_dir)
    out, node = run_with_cache(builder, data_dir=data_dir)
    # check output
    print('out, node:', out, node)
    print('cache_source:', node.get_cache_source())
    print('hash', node.get_hash())
    print('_get_objects_to_hash', node._get_objects_to_hash())
    print('ignored attributes:', node._hash_ignored_attributes)
    print('===== code =====')
    print('hash:', voronoi_local_code.get_hash())
    print('objects to hash:', voronoi_local_code._get_objects_to_hash())
    print('ignored attributes:', voronoi_local_code._hash_ignored_attributes)
    print('===== structure =====')
    print('structure hash:', structure.get_hash())
    print('objects to hash:', structure._get_objects_to_hash())
    print('ignored attributes:', structure._hash_ignored_attributes)
    print('===== parameters =====')
    print('hash:', parameters.get_hash())
    print('objects to hash:', parameters._get_objects_to_hash())
    print('ignored attributes:', parameters._hash_ignored_attributes)
    assert node.get_cache_source() is not None


def test_vca_structure(aiida_profile, voronoi_local_code):
    """
    test for vca_structure behaviour
    """
    pass


def test_overwrite_alat_input(aiida_profile, voronoi_local_code):
    """
    test using 'use_alat_input' keyword in input parameters
    """
    pass


def test_voronoi_after_kkr(aiida_profile, voronoi_local_code, run_with_cache, nopytest=False):
    """
    test voronoi run from parent kkr calculation (e.g. to update to a higher lmax value)
    """
    from aiida.orm import Dict, load_node
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.calculations.voro import VoronoiCalculation

    # load necessary files from db_dump files
    import_with_migration(test_dir / 'files/db_dump_kkrcalc.tar.gz')

    # first load parent voronoi calculation
    kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

    # extract KKR parameter and remote_data folder
    params_kkr_parent = kkr_calc.inputs.parameters
    parent_calc_remote = kkr_calc.outputs.remote_folder

    # set computer options
    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}

    # increase LMAX value from previous run
    params = kkrparams(params_type='voronoi', **params_kkr_parent)
    params.set_multiple_values(LMAX=3)
    new_params = Dict(dict=params.get_dict())

    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = new_params
    builder.parent_KKR = parent_calc_remote

    # now run calculation (or use cached results)
    if not nopytest:
        out, node = run_with_cache(builder, data_dir=data_dir)
        print('cache_source:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)
    else:
        from aiida.engine import run_get_node
        out, node = run_get_node(builder)

    print(out, node)

    # extract output nodes
    out_dict = node.outputs.output_parameters
    print(out_dict.get_dict())
    ret = node.outputs.retrieved

    # check if LMAX was increased
    assert params_kkr_parent.get_dict().get('LMAX') < new_params.get_dict().get('LMAX')

    # check if parsing was successful
    assert out_dict.get_dict().get('parser_errors') == []

    # check if overwrite_potential file is present
    assert 'overwrite_potential' in ret.list_object_names()


def test_overwrite_potential(aiida_profile, voronoi_local_code):
    """
    test providing overwirte_potential input node which overwrites the starting potentai with the given input
    """
    pass
