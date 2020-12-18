#!/usr/bin/env python
# coding: utf-8

from __future__ import absolute_import
from __future__ import print_function
import pytest
from aiida_kkr.tests.dbsetup import *
from aiida.engine import WorkChain, if_, ToContext, submit

@pytest.mark.usefixtures("aiida_env")
class Test_dbs_workflow():
    """
    Tests for the kkr_startpot workflow
    """

    @pytest.mark.timeout(240, method='thread')
    def test_bs_wc_Cu(self):
       
        from aiida import get_version
        from aiida.orm import Code, load_node, Dict, StructureData
        from aiida.plugins import DataFactory
        from aiida.orm import Computer
        from aiida.orm.querybuilder import QueryBuilder
        from masci_tools.io.kkr_params import kkrparams
        from aiida_kkr.workflows.bs import kkr_bs_wc
        from numpy import array
        import numpy as np   
        
        print('AiiDA version: {}'.format(get_version()))  
        Dict = DataFactory('dict')
        StructureData = DataFactory('structure')
        
        # preparation for computer and code 
        kkr_codename = 'kkr'
        computername = 'claix18_init'
        
        wfbs = kkr_bs_wc.get_wf_default()

        wfbs['NPT2'] = 12
        wfbs['EMAX'] = 5
        wfbs['EMIN'] = 10
        wfbs['RCLUSTZ'] = 2.3
        wfbs['TEMPR'] = 50.0  
        
        params_bs = Dict(dict=wfbs)
        
        # for runing in cluster or remote computer and should be rewiten for diferent computer
        options1 = {'max_wallclock_seconds': 36000,'resources': 
                  {'tot_num_mpiprocs': 48, 'num_machines': 1},
                  'custom_scheduler_commands': 
                  '#SBATCH --account=jara0191\n\nulimit -s unlimited; export OMP_STACKSIZE=2g',
                  'withmpi':
                  True}
        
        # for runing in local computer        
        options2 = {'queue_name' : queuename, 'resources': {"num_machines": 1}, 'max_wallclock_seconds' : 5*60, 'withmpi' : False, 'custom_scheduler_commands' : ''}
        options = Dict(dict=options1)
        
        # The scf-workflow needs also the voronoi and KKR codes to be able to run the calulations
        KKRCode = Code.get_from_string(kkr_codename+'@'+computername)
        
        label = 'bs calc Cu bulk'
        descr = ' testing bs workflow for Cu bulk'

        from aiida.tools.importexport import import_data

        import_data('files/db_dump_bs/db_dump_kkrcalc_bs.tar.gz', silent=True)
        kkr_calc_remote = load_node('d5782162-8393-4212-9340-c8ee8b725474').outputs.remote_folder
        
        builder = kkr_bs_wc.get_builder()
        builder.metadata.description = descr
        builder.metadata.label = label
        builder.kkr = KKRCode#kkrhost_local_code
        builder.wf_parameters = params_bs
        builder.options = options
        builder.remote_data = kkr_calc_remote

        from aiida.engine import run

        out = run(builder)
        n = out['results_wf']
        n = n.get_dict()
        assert n.get('successful')
        assert n.get('list_of_errors') == []
        d = out['BS_Data']
        kpts = d.get_array('Kpts')
        eng_points = d.get_array('energy_points')
        BlochSpectral = d.get_array('BlochSpectralFunction')
        
        # Loading the data from the previous calc (run manually) with the same config 
        test_BSF = np.loadtxt('files/db_dump_bs/test_spectral')
        test_Kpts = np.loadtxt('files/db_dump_bs/test_Kpts')
        test_eng = np.loadtxt('files/db_dump_bs/test_energy_points')
        # Define the error boundary
        dos_limt = abs(np.amin(BlochSpectral)*1E-2)
        eng_min = abs(np.amin(eng_points)*1E-7)

        differ_BSF = np.subtract(BlochSpectral,test_BSF)
        shape_differ_BSF = np.shape(differ_BSF)
         # shape(differ_BSF) =~ (kpts*eng) =~ (y*x) 
        for i in range(shape_differ_BSF[0]): 
            
            assert sum(abs(kpts[i,:] - test_Kpts[i,:] )) < 10**-7
            for j in range(shape_differ_BSF[1]):    
                # here to inspect the density validity
                assert dos_limt > abs(differ_BSF[i,j])
                if i==0 :
                    # here to inspect the energy validity     
                    assert eng_min > abs(eng_points[j] -test_eng[j])
                    

