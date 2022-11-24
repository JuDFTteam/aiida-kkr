#!/usr/bin/env python

from builtins import object
import pytest
from ..dbsetup import *
from ..conftest import import_with_migration
from aiida.engine import run_get_node


# tests
class Test_KKRnano_Calculation(object):
    """
    Tests for the KKRnanoCalculation
    """

    def test_kkrnano_from_voronoi(self, kkrhost_local_code, dry_run=True):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrnano import KKRnanoCalculation

        # load necessary files from db_dump files
        imported_nodes = import_with_migration('files/db_dump_KKRnano.aiida')

        # first load parent voronoi calculation
        voro_calc = load_node('c7bada4d-82ea-4d28-a4d4-55aaba918291')

        # set KKRnano input parameters
        voro_out = voro_calc.outputs.output_parameters.get_dict()
        voro_para = voro_calc.inputs.parameters.get_dict()
        params = {
            'bzdivide': [10, 10, 10],
            'emin': voro_out['emin'],
            'emax': 1.,
            'npnt1': 3,
            'npnt2': 10,
            'npnt3': 4,
            'npol': 7,
            'scfsteps': 1,
            'imix': 0,
            'mixing': 0.01,
            'rmax': 10.,
            'gmax': 100.,
            'nsra': 2,
            'kte': 1,
            'rclust_voronoi': voro_para['RCLUSTZ']
        }
        params_node = Dict(params)

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KKRnanoCalculation.get_builder()
        builder.code = kkrhost_local_code  # TODO replace with KKRnano code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.outputs.remote_folder
        builder.metadata.dry_run = dry_run
        out, node = run_get_node(builder)

    def test_kkrnano_from_kkrnano(self, kkrhost_local_code, dry_run=True):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrnano import KKRnanoCalculation

        # load necessary files from db_dump files
        imported_nodes = import_with_migration('files/db_dump_KKRnano.aiida')

        # first load parent voronoi calculation
        kkrnano_calc = load_node('da234580-b279-490b-972a-cdf9254c9278')

        # extract and update KKR parameter (add missing values)
        para_parent = kkrnano_calc.inputs.parameters.get_dict()
        params = {
            k: para_parent[k] for k in [
                'bzdivide', 'emin', 'emax', 'npnt1', 'npnt2', 'npnt3', 'npol', 'rmax', 'gmax', 'nsra', 'kte',
                'rclust_voronoi'
            ]
        }
        # maybe update some input parameters
        params['scfsteps'] = 3
        params['imix'] = 5
        params['mixing'] = 0.01
        params_node = Dict(params)

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KKRnanoCalculation.get_builder()
        builder.code = kkrhost_local_code  # TODO replace with KKRnano code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkrnano_calc.outputs.remote_folder
        builder.metadata.dry_run = dry_run
        out, node = run_get_node(builder)
