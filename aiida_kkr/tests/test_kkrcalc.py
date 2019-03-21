#!/usr/bin/env python

from builtins import object
import pytest

# some global settings

codename = 'KKRhost@iff003'
queuename = 'th1_node'
eps = 10**-14 # threshold for float comparison equivalence

from .test_vorocalc import wait_for_it

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkr_calculation(object):
    """
    Tests for the kkr calculation
    """
    
    def test_kkr_from_voronoi(self):
        """
        simple Cu noSOC, FP, lmax2 full example 
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        ParameterData = DataFactory('parameter')

        # load necessary files from db_dump files
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_vorocalc.tar.gz')
        import_data('files/db_dump_kkrcalc.tar.gz')

        # first load parent voronoi calculation       
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inp.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = ParameterData(dict=params.get_dict())

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(codename)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.out.remote_folder
        builder.submit_test()

    

    def test_kkr_from_kkr(self):
        """
        continue KKR calculation after a previous KKR calculation instead of starting from voronoi
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        ParameterData = DataFactory('parameter')

        # first load parent voronoi calculation       
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inp.parameters

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(codename)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.out.remote_folder
        builder.submit_test()
       
    
    def test_kkrflex(self):
        """
        test kkrflex file writeout (GF writeout for impurity calculation)
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        ParameterData = DataFactory('parameter')

        # first load parent voronoi calculation       
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add KKRFLEX option)
        params_node = kkr_calc.inp.parameters
        params = params_node.get_dict()
        params['RUNOPT'] = ['KKRFLEX']
        params_node = ParameterData(dict=params)
        
        # create an impurity_info node
        imp_info = ParameterData(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[29.]})

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(codename)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.out.remote_folder
        builder.impurity_info = imp_info
        builder.submit_test()
       
 
#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkr_calculation()
   Test.test_kkr_from_voronoi()
   Test.test_kkr_from_kkr()
   Test.test_kkrflex()
