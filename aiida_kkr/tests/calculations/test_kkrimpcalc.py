#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida.manage.tests.pytest_fixtures import aiida_profile, clear_database, clear_database_after_test

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
        from aiida.tools.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz', silent=True)
        GF_host_calc = load_node('baabef05-f418-4475-bba5-ef0ee3fd5ca6')

        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.outputs.retrieved.open('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = Dict(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.outputs.remote_folder)

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
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


#run test manually
if __name__=='__main__':
   from aiida import load_profile
   load_profile()
   Test = Test_kkrimp_calculation()
   Test.test_host_in_host()
