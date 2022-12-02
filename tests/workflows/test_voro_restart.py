#!/usr/bin/env python

import pytest
from aiida.manage.tests.pytest_fixtures import aiida_profile
from aiida_testing.export_cache._fixtures import export_cache
from ..dbsetup import queuename
from ..conftest import voronoi_local_code


# tests
def test_voronoi_base_restart(aiida_profile, voronoi_local_code):
    """
    simple Cu noSOC, FP, lmax2 full example
    """
    from aiida.orm import Dict, StructureData
    from masci_tools.io.kkr_params import kkrparams
    from aiida_kkr.workflows.restart.restart_voro import VoronoiCalculationBaseWorkChain

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
    builder = VoronoiCalculationBaseWorkChain.get_builder()
    builder.calc.code = voronoi_local_code
    builder.calc.metadata.options = options
    builder.calc.parameters = ParaNode
    builder.calc.structure = Cu
    from aiida.engine import run
    run(builder)
