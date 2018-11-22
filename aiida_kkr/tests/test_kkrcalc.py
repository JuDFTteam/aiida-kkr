#!/usr/bin/env python

import pytest

# some global settings

codename = 'KKRhost@iff003'
queuename = 'th1_node'
eps = 10**-14 # threshold for float comparison equivalence

from test_vorocalc import wait_for_it

# tests
@pytest.mark.usefixtures("aiida_env")
class Test_kkr_calculation():
    """
    Tests for the kkr calculation
    """
    
    def test_kkr_from_voronoi(self):
        """
        simple Cu noSOC, FP, lmax2 full example 
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        ParameterData = DataFactory('parameter')

        # first load parent voronoi calculation       
        voro_calc = load_node('559b9d9b-3525-402e-9b24-ecd8b801853c')

        # extract and update KKR parameter (add missing values)
        params = kkrparams(**voro_calc.inp.parameters.get_dict())
        params.set_multiple_values(RMAX=7., GMAX=65.)
        params_node = ParameterData(dict=params.get_dict())

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(codename)
       
        calc = code.new_calc()
        calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})
        calc.use_parameters(params_node)
        calc.set_queue_name(queuename)
        calc.use_parent_folder(voro_calc.out.remote_folder)
       
        # now store all nodes and submit calculation
        calc.store_all()
        calc.submit()

        # now wait for the calculation to finish
        wait_for_it(calc)

        # finally check some output
        print '\n\ncheck values ...\n-------------------------------------------------'
       
        test_ok = calc.get_state() == u'FINISHED'
        print 'calculation state', calc.get_state(), 'OK?', test_ok
        assert test_ok
       
        test_ok = calc.res.alat_internal == 4.82381975
        print 'alat internal units', calc.res.alat_internal, 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.charge_valence_states_per_atom[0]-10.965496)<eps
        print 'valence charge', calc.res.charge_valence_states_per_atom, 'OK?', test_ok
        assert test_ok
       
        test_ok = not calc.res.convergence_group.get('calculation_converged')
        print 'calculation_converged', calc.res.convergence_group.get('calculation_converged'), 'OK?', test_ok
        assert test_ok
       
        test_ok = calc.res.convergence_group.get('nsteps_exhausted')
        print 'nsteps_exhausted', calc.res.convergence_group.get('nsteps_exhausted'), 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.convergence_group.get('charge_neutrality')+0.034504)<eps
        print 'charge neutrality', calc.res.convergence_group.get('charge_neutrality'), 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.convergence_group.get('rms')-0.178)<eps
        print 'rms', calc.res.convergence_group.get('rms'), 'OK?', test_ok
        assert test_ok
       
        print '\ndone with checks\n'
    

    def test_kkr_from_kkr(self):
        """
        continue KKR calculation after a previous KKR calculation instead of starting from voronoi
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
        ParameterData = DataFactory('parameter')

        # first load parent voronoi calculation       
        kkr_calc = load_node('3058bd6c-de0b-400e-aff5-2331a5f5d566')

        # extract KKR parameter (add missing values)
        params_node = kkr_calc.inp.parameters

        # load code from database and create new voronoi calculation
        code = Code.get_from_string(codename)
       
        calc = code.new_calc()
        calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})
        calc.use_parameters(params_node)
        calc.set_queue_name(queuename)
        calc.use_parent_folder(kkr_calc.out.remote_folder)
       
        # now store all nodes and submit calculation
        calc.store_all()
        calc.submit()

        # now wait for the calculation to finish
        wait_for_it(calc)

        # finally check some output
        print '\n\ncheck values ...\n-------------------------------------------------'
       
        test_ok = calc.get_state() == u'FINISHED'
        print 'calculation state', calc.get_state(), 'OK?', test_ok
        assert test_ok
       
        test_ok = calc.res.alat_internal == 4.82381975
        print 'alat internal units', calc.res.alat_internal, 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.charge_valence_states_per_atom[0]-10.972076)<eps
        print 'valence charge', calc.res.charge_valence_states_per_atom, 'OK?', test_ok
        assert test_ok
       
        test_ok = not calc.res.convergence_group.get('calculation_converged')
        print 'calculation_converged', calc.res.convergence_group.get('calculation_converged'), 'OK?', test_ok
        assert test_ok
       
        test_ok = calc.res.convergence_group.get('nsteps_exhausted')
        print 'nsteps_exhausted', calc.res.convergence_group.get('nsteps_exhausted'), 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.convergence_group.get('charge_neutrality')+0.027924)<eps
        print 'charge neutrality', calc.res.convergence_group.get('charge_neutrality'), 'OK?', test_ok
        assert test_ok
       
        test_ok = abs(calc.res.convergence_group.get('rms')-0.17427)<eps
        print 'rms', calc.res.convergence_group.get('rms'), 'OK?', test_ok
        assert test_ok
       
        print '\ndone with checks\n'

    
    def test_qdos(self):
        """
        ...
        """
        pass

    
    def test_kkrflex(self):
        """
        test kkrflex file writeout (GF writeout for impurity calculation)
        """
        from aiida.orm import Code, load_node, DataFactory
        from masci_tools.io.kkr_params import kkrparams
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
       
        calc = code.new_calc()
        calc.set_resources({'num_machines':1, 'tot_num_mpiprocs':1})
        calc.use_parameters(params_node)
        calc.set_queue_name(queuename)
        calc.use_parent_folder(kkr_calc.out.remote_folder)
        calc.use_impurity_info(imp_info)
       
        # now store all nodes and submit calculation
        calc.store_all()
        calc.submit()

        # now wait for the calculation to finish
        wait_for_it(calc)

        # finally check some output
        print '\n\ncheck values ...\n-------------------------------------------------'
       
        test_ok = calc.get_state() == u'FINISHED'
        print 'calculation state', calc.get_state(), 'OK?', test_ok
        assert test_ok

        # check if kkrflex file have been retreived
        path = calc.out.retrieved.get_abs_path('')
        import os
        retrieved_files = os.listdir(path)
        for filename in 'kkrflex_green kkrflex_tmat kkrflex_atominfo kkrflex_intercell_ref kkrflex_intercell_cmoms'.split():
            print '{} file there?'.format(filename)
            assert filename in retrieved_files
            print 'OK'
       
        print '\ndone with checks\n'

    
    def test_vca(self):
        """
        ...
        """
        pass

    
    def test_use_alat_input(self):
        """
        ...
        """
        pass

    
    def test_set_ef_value(self):
        """
        ...
        """
        pass

    
    def test_dos(self):
        """
        ...
        """
        pass

    
    def test_jij(self):
        """
        ...
        """
        pass

 
#run test manually
if __name__=='__main__':
   from aiida import is_dbenv_loaded, load_dbenv
   if not is_dbenv_loaded():
      load_dbenv()
   Test = Test_kkr_calculation()
   Test.test_kkrflex()
