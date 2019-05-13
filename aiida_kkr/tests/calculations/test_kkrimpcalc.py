#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest
from aiida_kkr.tests.calculations.test_vorocalc import wait_for_it
from aiida_kkr.tests.dbsetup import *

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkrimp_calculation(object):
    """
    Tests for the kkrimp calculation
    """

    def test_host_in_host(self):
        """
        simple Cu noSOC, FP, lmax2
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrimp import KkrimpCalculation
        Dict = DataFactory('dict')

        # first load parent voronoi calculation
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('de9b5093-25e7-407e-939e-9282c4431343') #'9459b4ea-ead5-4268-aa29-1c5e18654d77')

        # prepare computer and code (needed so that
        prepare_code(kkrimp_codename, codelocation, computername, workdir)

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
        #kkrimp_code = Code.get_from_string(codename)
        kkrimp_code = Code.get_from_string(kkrimp_codename+'@'+computername)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrimpCalculation.get_builder()
        builder.code = kkrimp_code
        builder.host_Greenfunction_folder = GF_host_calc.outputs.remote_folder
        builder.impurity_potential = startpot_imp_sfd
        builder.metadata.options = options
        builder.parameters = ParamsKKRimp
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkrimp_calculation()
   Test.test_host_in_host()
