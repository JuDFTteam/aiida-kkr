#!/usr/bin/env python

from builtins import object
import pytest
import pathlib
from aiida_test_cache.archive_cache import enable_archive_cache
from aiida.engine import run_get_node
from ..dbsetup import *
from ..conftest import voronoi_local_code, test_dir, data_dir, import_with_migration

kkr_codename = 'kkrhost'


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
    ParaNode = Dict(params.get_dict())

    options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = ParaNode
    builder.structure = Cu
    builder.metadata.dry_run = True
    from aiida.engine import run
    run(builder)


def test_voronoi_cached(aiida_profile_clean, voronoi_local_code, enable_archive_cache):
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
    parameters = Dict({k: v for k, v in kkr_params.items() if v})

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
    with enable_archive_cache(data_dir / 'voronoi_cached.aiida'):
        out, node = run_get_node(builder)

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


def test_voronoi_after_kkr(aiida_profile_clean, voronoi_local_code, enable_archive_cache, nopytest=False):
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
    new_params = Dict(params.get_dict())

    builder = VoronoiCalculation.get_builder()
    builder.code = voronoi_local_code
    builder.metadata.options = options
    builder.parameters = new_params
    builder.parent_KKR = parent_calc_remote

    # now run calculation (or use cached results)
    if not nopytest:
        with enable_archive_cache(data_dir / 'voronoi_after_kkr.aiida'):
            out, node = run_get_node(builder)
        print('hash:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)
    else:
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
