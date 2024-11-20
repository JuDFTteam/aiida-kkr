#!/usr/bin/env python

from builtins import object
import pytest
from aiida.engine import run_get_node
from ..dbsetup import *
from ..conftest import data_dir, import_with_migration


# tests
class Test_kkrimp_calculation(object):
    """
    Tests for the kkrimp calculation
    """

    def test_host_in_host(self, aiida_profile, kkrimp_local_code):
        """
        simple Cu noSOC, FP, lmax2
        """
        from aiida.orm import Code, load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrimp import KkrimpCalculation

        # first load parent voronoi calculation
        import_with_migration('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools import neworder_potential_wf
        from numpy import loadtxt
        with GF_host_calc.outputs.retrieved.open('scoef') as _f:
            neworder_pot1 = [int(i) for i in loadtxt(_f, skiprows=1)[:, 3] - 1]
        settings_dict = {'pot1': 'out_potential', 'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(settings_dict)
        startpot_imp_sfd = neworder_potential_wf(
            settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder
        )

        # set 1 simple mixing step
        kkrimp_params = kkrparams(params_type='kkrimp')
        kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
        ParamsKKRimp = Dict(kkrimp_params.get_dict())

        # create new KKRimp calculation
        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrimpCalculation.get_builder()
        builder.code = kkrimp_local_code
        builder.host_Greenfunction_folder = GF_host_calc.outputs.remote_folder
        builder.impurity_potential = startpot_imp_sfd
        builder.metadata.options = options
        builder.parameters = ParamsKKRimp
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)

    def test_cached(self, aiida_profile, kkrimp_local_code, enable_archive_cache):
        """
        simple Cu noSOC, FP, lmax2
        """
        from aiida.orm import Code, load_node, Dict
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrimp import KkrimpCalculation

        # first load from imported calculation
        import_with_migration('files/export_kkrimp_full.tar.gz')
        kkrimp_calc = load_node('81a7c87a-bddb-4ea0-9ce6-31b725ea02e2')
        startpot_imp_sfd = kkrimp_calc.inputs.impurity_potential
        host_Greenfunction_folder = kkrimp_calc.inputs.host_Greenfunction_folder

        # set 1 simple mixing step
        kkrimp_params = kkrparams(params_type='kkrimp')
        kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
        ParamsKKRimp = Dict(dict=kkrimp_params.get_dict())

        # create new KKRimp calculation
        options = {'resources': {'num_machines': 1, 'tot_num_mpiprocs': 1}, 'queue_name': queuename}
        builder = KkrimpCalculation.get_builder()
        builder.code = kkrimp_local_code
        builder.host_Greenfunction_folder = host_Greenfunction_folder
        builder.impurity_potential = startpot_imp_sfd
        builder.metadata.options = options
        builder.parameters = ParamsKKRimp

        with enable_archive_cache(data_dir / 'kkrimp_cached.aiida'):
            out, node = run_get_node(builder)
        print(out, node)
        print('cache_source:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)
        assert node.get_cache_source() is not None
        print(out['output_parameters'].get_dict())
        assert out['output_parameters']['parser_errors'] == []


#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    Test = Test_kkrimp_calculation()
    Test.test_host_in_host()
