#!/usr/bin/env python

from __future__ import absolute_import
from builtins import object
import pytest
from aiida_kkr.tests.calculations.test_vorocalc import wait_for_it
from aiida_kkr.tests.dbsetup import *

# some global settings
eps = 10**-14 # threshold for float comparison equivalence

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
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        Dict = DataFactory('dict')

        # load necessary files from db_dump files
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_vorocalc.tar.gz')
        import_data('files/db_dump_kkrcalc.tar.gz')

        # prepare computer and code (needed so that
        prepare_code(kkr_codename, codelocation, computername, workdir)

        # first load parent voronoi calculation
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inputs.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = Dict(dict=params.get_dict())

        # load code from database and create new voronoi calculation
        #code = Code.get_from_string(codename)
        code = Code.get_from_string(kkr_codename+'@'+computername)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = voro_calc.outputs.remote_folder
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


    def test_kkr_from_kkr(self):
        """
        continue KKR calculation after a previous KKR calculation instead of starting from voronoi
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        Dict = DataFactory('dict')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(kkr_codename+'@'+computername)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


    def test_kkrflex(self):
        """
        test kkrflex file writeout (GF writeout for impurity calculation)
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        Dict = DataFactory('dict')

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add KKRFLEX option)
        params_node = kkr_calc.inputs.parameters
        params = params_node.get_dict()
        params['RUNOPT'] = ['KKRFLEX']
        params_node = Dict(dict=params)

        # create an impurity_info node
        imp_info = Dict(dict={'Rcut':1.01, 'ilayer_center': 0, 'Zimp':[29.]})

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(kkr_codename+'@'+computername)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.impurity_info = imp_info
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


    def test_kkr_qdos(self):
        """
        run bandstructure calculation
        """
        from aiida.orm import Code, load_node
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.kkr import KkrCalculation
        Dict = DataFactory('dict')
        KpointsData = DataFactory('array.kpoints')

        # define k-path
        kpoints = KpointsData()
        kpoints.set_kpoints([[0,0,0],[0.1,0,0],[0.2,0,0],[0.3,0,0],[0.4,0,0]])

        # first load parent voronoi calculation
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inputs.parameters

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(kkr_codename+'@'+computername)
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = KkrCalculation.get_builder()
        builder.code = code
        builder.metadata.options = options
        builder.parameters = params_node
        builder.parent_folder = kkr_calc.outputs.remote_folder
        builder.kpoints = kpoints
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)


#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkr_calculation()
   Test.test_kkr_from_voronoi()
   Test.test_kkr_from_kkr()
   Test.test_kkrflex()
   Test.test_kkr_qdos()
