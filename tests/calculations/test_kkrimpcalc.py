#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest
from ..dbsetup import *
from aiida.manage.tests.pytest_fixtures import aiida_profile, clear_database, clear_database_after_test
from aiida_testing.export_cache._fixtures import run_with_cache, export_cache, load_cache, hash_code_by_entrypoint, absolute_archive_path
from ..conftest import data_dir, import_with_migration


# tests
class Test_kkrimp_calculation(object):
    """
    Tests for the kkrimp calculation
    """

    def test_host_in_host(self, clear_database, kkrimp_local_code):
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
        settings = Dict(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(
            settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder
        )

        # set 1 simple mixing step
        kkrimp_params = kkrparams(params_type='kkrimp')
        kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
        ParamsKKRimp = Dict(dict=kkrimp_params.get_dict())

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


'''
    def test_cached(self, clear_database, kkrimp_local_code, run_with_cache):
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
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)
        # need to load startpot to use caching
        #startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder)
        #from aiida.tools.importexport import export
        #export([startpot_imp_sfd], outfile='files/export_startpot_kkrimp.zip', overwrite=True, include_comments=False)
        #print('exported startpot:', startpot_imp_sfd)
        import_with_migration('files/export_startpot_kkrimp.zip')
        startpot_imp_sfd = load_node('085ee1f1-da92-464f-a883-edeb5900de94')



        # set 1 simple mixing step
        kkrimp_params = kkrparams(params_type='kkrimp')
        kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
        ParamsKKRimp = Dict(dict=kkrimp_params.get_dict())

        # create new KKRimp calculation
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrimpCalculation.get_builder()
        builder.code = kkrimp_local_code
        builder.host_Greenfunction_folder = GF_host_calc.outputs.remote_folder
        builder.impurity_potential = startpot_imp_sfd
        builder.metadata.options = options
        builder.parameters = ParamsKKRimp
        out, node = run_with_cache(builder, data_dir=data_dir)
        print(out, node)
        print('cache_source:', node.get_hash())
        print('cache_source:', node.get_cache_source())
        print('code objects to hash:', node._get_objects_to_hash())
        print('ignored attributes:', node._hash_ignored_attributes)
        assert node.get_cache_source() is not None
        #print(out['output_parameters'].get_dict())
        assert out['output_parameters']['parser_errors'] == []
'''

#run test manually
if __name__ == '__main__':
    from aiida import load_profile
    load_profile()
    Test = Test_kkrimp_calculation()
    Test.test_host_in_host()
