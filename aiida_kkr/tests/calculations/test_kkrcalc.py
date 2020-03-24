#!/usr/bin/env python

from __future__ import absolute_import
from __future__ import print_function
from builtins import object
import pytest
from aiida.engine import run, run_get_node
from aiida_kkr.tests.dbsetup import *
from ..conftest import kkrhost_local_code
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint
from aiida.manage.tests.pytest_fixtures import aiida_local_code_factory, aiida_localhost, temp_dir, aiida_profile, clear_database, clear_database_after_test

# some global settings
eps = 10**-14 # threshold for float comparison equivalence

#kkr_codename = 'kkrhost_intel19'
kkr_codename = 'kkrhost'

dry_run = True


# tests
class Test_kkr_calculation(object):
    """
    Tests for the kkr calculation
    """

    def test_kkr_from_voronoi(self, kkrhost_local_code):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_vorocalc.tar.gz', extras_mode_existing='nnl')

        # first load parent voronoi calculation
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inputs.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = Dict(dict=params.get_dict())

        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.outputs.remote_folder
        builder.metadata.dry_run = dry_run
        out, node = run_get_node(builder)


    def test_kkr_cached(self, aiida_profile, kkrhost_local_code, run_with_cache):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_vorocalc.tar.gz', extras_mode_existing='nnl')

        # first load parent voronoi calculation
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inputs.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = Dict(dict=params.get_dict())

        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.outputs.remote_folder
        builder.metadata.dry_run = False
        print(options)
        print(params_node)
        print(kkrhost_local_code)
        print(voro_calc)
        print(builder)
        out, node = run_with_cache(builder)
        print((node, out))
        print((node.get_cache_source()))
        assert node.get_cache_source() is not None
        out_dict = node.outputs.output_parameters.get_dict()
        from pprint import pprint 
        pprint(out_dict)


    def test_kkr_from_kkr(self, kkrhost_local_code):
        """
        continue KKR calculation after a previous KKR calculation instead of starting from voronoi
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.metadata.dry_run = dry_run
        out = run(builder)
        print(out)


    def test_kkrflex(self, kkrhost_local_code):
        """
        test kkrflex file writeout (GF writeout for impurity calculation)
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add KKRFLEX option)
        params_node = kkr_calc.inputs.parameters
        params = params_node.get_dict()
        params['RUNOPT'] = ['KKRFLEX']
        params_node = Dict(dict=params)

        # create an impurity_info node
        imp_info = Dict(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[29.]})

        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.impurity_info = imp_info
        builder.metadata.dry_run = dry_run
        out = run(builder)
        print(out)


    def test_kkr_qdos(self, kkrhost_local_code):
        """
        run bandstructure calculation
        """
        from aiida.orm import load_node, Dict, KpointsData
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # define k-path
        kpoints = KpointsData()
        kpoints.set_kpoints([[0,0,0],[0.1,0,0],[0.2,0,0],[0.3,0,0],[0.4,0,0]])
        kpoints.set_cell([[1.0,0,0],[0,1.0,0],[0,0,1.0]])

        # load necessary files from db_dump files
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        # construct process builder and run calc
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.kpoints = kpoints
        builder.metadata.dry_run = dry_run
        out = run(builder)
        print(out)


    def test_kkr_increased_lmax(self, kkrhost_local_code, run_with_cache):
        """
        run kkr calculation from output of previous calculation but with increased lmax
        (done with auxiliary voronoi calculation which is imported here).
        """
        from aiida.orm import load_node
        from aiida_kkr.calculations import KkrCalculation

        # import previous voronoi calc (ran with parent_KKR mode and increased LMAX in input params)
        from aiida.tools.importexport import import_data
        import_data('data_dir/VoronoiCalculation-nodes-8c7aed435f2140768f52c78b0b1b0629.tar.gz')
        voro_with_kkr_input = load_node('69441815-8d55-4412-baf6-1793665aba19')

        # extract KKR parameter from imported voronoi calc
        params_node = voro_with_kkr_input.inputs.parameters

        # construct process builder and run calc
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_with_kkr_input.outputs.remote_folder

        # now run or load from cached data
        out, node = run_with_cache(builder)
        print('cache_source:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)

        # inspect result
        out_dict = node.outputs.output_parameters.get_dict()
        assert len(out_dict['parser_errors']) < 2
        assert node.inputs.parameters.get_dict().get('LMAX') == 3
        # check if parent voronoi calculation had parent_KKR input
        from aiida.orm import RemoteData
        input_remote = node.get_incoming(node_class=RemoteData).first().node
        v = input_remote.get_incoming().first().node
        assert 'parent_KKR' in [i.link_label for i in v.get_incoming()]

