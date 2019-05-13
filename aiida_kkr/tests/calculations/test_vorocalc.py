#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from builtins import object
from aiida_kkr.tests.dbsetup import *
import pytest

#TODO
# implement missing tests:
# * test_vca_structure
# * test_overwrite_alat_input
# * test_voronoi_after_kkr
# * test_overwrite_potential

def wait_for_it(calc, maxwait=300, dT=10):
    """
    helper function used to wait until calculation reaches FINISHED state
    wait for maximally <maxwait> seconds and check the calculation's state every <dT> seconds
    """
    from time import sleep
    nsteps = maxwait/dT
    print('waiting for calculation to finish (maximally wait for {} seconds)'.format(maxwait))
    istep = 0
    calcstate = u'UNKNOWN'
    while istep < nsteps:
        print('checking status')
        sleep(dT)
        calcstate = calc.get_state()
        istep += 1
        if calcstate == u'FINISHED' or calcstate == u'FAILED':
            break

    if calcstate == u'FINISHED':
        print('calculation reached FINISHED state')
    elif calcstate == u'FAILED':
        print('calculation in FAILED state')
    else:
        print('maximum waiting time exhausted')


# tests
@pytest.mark.usefixtures("aiida_env")
class Test_voronoi_calculation(object):
    """
    Tests for the voronoi calculation
    """
    def test_startpot_Cu_simple(self):
        """
        simple Cu noSOC, FP, lmax2 full example
        """
        from aiida.orm import Code
        from aiida.plugins import DataFactory
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.calculations.voro import VoronoiCalculation

        Dict = DataFactory('dict')
        StructureData = DataFactory('structure')

        # create StructureData instance for Cu
        alat = 3.61 # lattice constant in Angstroem
        bravais = [[0.5*alat, 0.5*alat, 0], [0.5*alat, 0, 0.5*alat], [0, 0.5*alat, 0.5*alat]] # Bravais matrix in Ang. units
        Cu = StructureData(cell=bravais)
        Cu.append_atom(position=[0,0,0], symbols='Cu')

        # create Dict input node using kkrparams class from masci-tools
        params = kkrparams(params_type='voronoi')
        params.set_multiple_values(LMAX=2, NSPIN=1, RCLUSTZ=2.3)
        Dict = DataFactory('dict') # use DataFactory to get ParamerterData class
        ParaNode = Dict(dict=params.get_dict())

        # import computer etc from database dump
        from aiida.orm.importexport import import_data
        import_data('files/db_dump_vorocalc.tar.gz')

        # prepare computer and code (needed so that
        prepare_code(voro_codename, codelocation, computername, workdir)

        # load code from database and create new voronoi calculation
        #code = Code.get_from_string(codename)
        code = Code.get_from_string(voro_codename+'@'+computername)

        #code = Code.get_from_string('voronoi@localhost_new')
        options = {'resources': {'num_machines':1, 'tot_num_mpiprocs':1}, 'queue_name': queuename}
        builder = VoronoiCalculation.get_builder()
        builder.code = code
        builder.metadata.options = options
        builder.parameters = ParaNode
        builder.structure = Cu
        builder.metadata.dry_run = True
        from aiida.engine import run
        run(builder)

    def test_vca_structure(self):
        """
        test for vca_structure behaviour
        """
        pass

    def test_overwrite_alat_input(self):
        """
        test using 'use_alat_input' keyword in input parameters
        """
        pass

    def test_voronoi_after_kkr(self):
        """
        test voronoi run from parent kkr calculation (e.g. to update to a higher lmax value)
        """
        pass

    def test_overwrite_potential(self):
        """
        test providing overwirte_potential input node which overwrites the starting potentai with the given input
        """
        pass


#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_voronoi_calculation()
   Test.test_startpot_Cu_simple()
