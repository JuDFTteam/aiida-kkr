#!/usr/bin/env python

from builtins import object
import pytest
from aiida.engine import run, run_get_node
from ..dbsetup import *
from ..conftest import kkrhost_local_code, test_dir, data_dir, import_with_migration

# some global settings
eps = 10**-14  # threshold for float comparison equivalence

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
        import_with_migration('files/db_dump_vorocalc.tar.gz')

        # first load parent voronoi calculation
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inputs.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = Dict(params.get_dict())

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.outputs.remote_folder
        builder.metadata.dry_run = dry_run
        out, node = run_get_node(builder)

    def test_kkr_cached(self, aiida_profile, kkrhost_local_code, enable_archive_cache):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        import_with_migration('files/db_dump_vorocalc.tar.gz')

        # first load parent voronoi calculation
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inputs.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = Dict(params.get_dict())

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
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
        with enable_archive_cache(data_dir / 'kkr_cached.aiida'):
            out, node = run_get_node(builder)
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
        import_with_migration('files/db_dump_kkrcalc.tar.gz')
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
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
        from aiida.orm import load_node, Dict, Bool
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        import_with_migration('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add KKRFLEX option)
        params_node = kkr_calc.inputs.parameters
        params = params_node.get_dict()
        params['RUNOPT'] = ['KKRFLEX']
        params_node = Dict(params)

        # create an impurity_info node
        imp_info = Dict({'Rcut': 1.01, 'ilayer_center': 0, 'Zimp': [29.]})

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.impurity_info = imp_info
        builder.retrieve_kkrflex = Bool(False)
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
        kpoints.set_kpoints([[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0], [0.4, 0, 0]])
        kpoints.set_cell([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])

        # load necessary files from db_dump files
        import_with_migration('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        # construct process builder and run calc
        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.kpoints = kpoints
        builder.metadata.dry_run = dry_run
        out = run(builder)
        print(out)

    def test_kkr_increased_lmax(self, clear_database_before_test, kkrhost_local_code, enable_archive_cache):
        """
        run kkr calculation from output of previous calculation but with increased lmax
        (done with auxiliary voronoi calculation which is imported here).
        """
        from aiida.orm import load_node, CalcJobNode, Dict
        from aiida_kkr.calculations import KkrCalculation, VoronoiCalculation
        from aiida_kkr.tools import kkrparams

        # import previous voronoi calc (ran with parent_KKR mode and increased LMAX in input params)
        import_with_migration('files/export_kkr_lmax_change.tar.gz')
        voro_with_kkr_input = load_node('4d92dc05-041f-422c-945d-dec1fb10301e')

        # extract KKR parameter from imported voronoi calc
        params_node = voro_with_kkr_input.inputs.parameters
        p = kkrparams(
            **{
                k: v
                for k, v in params_node.get_dict().items()
                if k not in ['<NEWVERSION_BDG>', '<DECOUPLE_SPINS_CHEBY>']
            }
        )

        # construct process builder and run calc
        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = Dict(p)
        builder.parent_folder = voro_with_kkr_input.outputs.remote_folder

        # now run or load from cached data
        with enable_archive_cache(data_dir / 'kkr_increased_lmax.aiida'):
            out, node = run_get_node(builder)
        print('cache_source:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)

        print('output files', node.outputs.retrieved.list_object_names())
        print('std.out')
        with node.outputs.retrieved.open('_scheduler-stdout.txt') as f:
            print(f.readlines())
        print('std.err')
        with node.outputs.retrieved.open('_scheduler-stderr.txt') as f:
            print(f.readlines())
        print('inputcard')
        with node.outputs.retrieved.open('inputcard') as f:
            print(f.readlines())
        print('out_kkr')
        with node.outputs.retrieved.open('out_kkr') as f:
            print(f.readlines())
        print('output.000.txt')
        with node.outputs.retrieved.open('output.000.txt') as f:
            print(f.readlines())

        # inspect result
        out_dict = node.outputs.output_parameters.get_dict()
        assert len(out_dict['parser_errors']) < 2
        assert node.inputs.parameters.get_dict().get('LMAX') == 3
        # check if parent voronoi calculation had parent_KKR input
        from aiida.orm import RemoteData
        input_remote = node.get_incoming(node_class=RemoteData).first().node
        v = input_remote.get_incoming().first().node
        assert 'parent_KKR' in [i.link_label for i in v.get_incoming()]

    def test_kkr_gf_writeout_full_impcls(self, kkrhost_local_code):
        """
        run kkr calculation from output of previous calculation but with increased lmax
        (done with auxiliary voronoi calculation which is imported here).
        """
        from aiida.orm import load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation

        # load necessary files from db_dump files
        import_with_migration('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add KKRFLEX option)
        params_node = kkr_calc.inputs.parameters
        params = params_node.get_dict()
        params['RUNOPT'] = ['KKRFLEX']
        params_node = Dict(params)

        # create an impurity_info node
        imp_info = Dict({'imp_cls': []})

        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = kkrhost_local_code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.impurity_info = imp_info
        builder.metadata.dry_run = dry_run
        out = run(builder)
        print(out)
