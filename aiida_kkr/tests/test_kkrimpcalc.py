#!/usr/bin/env python

from builtins import object
import pytest

# some global settings

codename = 'KKRimp@iff003'
queuename = 'th1_node'

from .test_vorocalc import wait_for_it

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
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkrimp import KkrimpCalculation
        ParameterData = DataFactory('parameter')

        # first load parent voronoi calculation       
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_kkrflex_create.tar.gz')
        GF_host_calc = load_node('de9b5093-25e7-407e-939e-9282c4431343') #'9459b4ea-ead5-4268-aa29-1c5e18654d77')
       
        # now create a SingleFileData node containing the impurity starting potential
        from aiida_kkr.tools.common_workfunctions import neworder_potential_wf
        from numpy import loadtxt
        neworder_pot1 = [int(i) for i in loadtxt(GF_host_calc.out.retrieved.get_abs_path('scoef'), skiprows=1)[:,3]-1]
        settings_dict = {'pot1': 'out_potential',  'out_pot': 'potential_imp', 'neworder': neworder_pot1}
        settings = ParameterData(dict=settings_dict)
        startpot_imp_sfd = neworder_potential_wf(settings_node=settings, parent_calc_folder=GF_host_calc.out.remote_folder)

        # set 1 simple mixing step
        kkrimp_params = kkrparams(params_type='kkrimp')
        kkrimp_params.set_multiple_values(SCFSTEPS=1, IMIX=0, MIXFAC=0.05)
        ParamsKKRimp = ParameterData(dict=kkrimp_params.get_dict())

        # create new KKRimp calculation
        kkrimp_code = Code.get_from_string(codename)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrimpCalculation.get_builder()
        builder.code = kkrimp_code
        builder.host_Greenfunction_folder = GF_host_calc.out.remote_folder
        builder.impurity_potential = startpot_imp_sfd
        builder.options = options
        builder.parameters = ParamsKKRimp
        builder.submit_test()

 
#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkrimp_calculation()
   Test.test_host_in_host()
