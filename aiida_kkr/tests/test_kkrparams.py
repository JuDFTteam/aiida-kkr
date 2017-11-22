#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 16:43:31 2017

@author: ruess
"""

import pytest
from aiida_kkr.tools.kkr_params import kkrparams
    
    
#class Test_kkrparams_class():
#    """automatic tests for kkrparams class (defines KKR parameter and gives methods to fill, check and write input for KKR, voronoi etc.)"""
    
    
class Test_create_and_set_keys():
    def test_create_params_with_inital_values(self):
        p = kkrparams(RBASIS=[0,0,0], params_type='voronoi')
        assert type(p)==kkrparams
        assert p.values['<RBASIS>'] == [0,0,0]
    
    def test_default_values(self):
        p = kkrparams()
        assert p.values['EMIN'] is None
        
    def test_set_single_value(self):
        p = kkrparams()
        p.set_value('EMIN', 2)
        assert p.values['EMIN'] == 2.
        assert p.values['EMAX'] is None
        
    def test_set_multiple_values(self):
        p = kkrparams()
        p.set_multiple_values(EMIN=1, EMAX=2)
        assert p.values['EMIN']== 1.
        assert p.values['EMAX']== 2.

    
class Test_capture_wrong_input():
    def test_wrong_input_type(self):
        p = kkrparams()
        known_error = False
        try:
            p.set_value('EMIN', '2')
        except TypeError:
            known_error = True
        assert known_error
        
        known_error = False
        try:
            p.set_value('EMIN', False)
        except TypeError:
            known_error = True
        assert known_error
        
    def test_wrong_input_array_dimension(self):
        p = kkrparams()
        from numpy import array, sqrt
        bravais = array([[0.7071067812, -0.5, 0.0], [0.7071067812, 0.5, 0.0], [sqrt(2), 0.0, 0.866025404]])
        
        # atom positions in relative coordinates
        basis_vectors = []
        for iatom in range(6):
            tmp = array([0, 0, 0])+iatom*array([0.5, 0.5, bravais[2, 2]])
            tmp[0] = tmp[0]%1
            tmp[1] = tmp[1]%1
            print(iatom, tmp)
            basis_vectors.append(tmp)
        basis_vectors = array(basis_vectors)
        p.set_value('INTERFACE', True)
        p.set_value('<RBLEFT>', array([[1,1],[0,1]]))
        
    def test_input_consistency_check_fail(self):
        knownError = False
        try: 
            p = kkrparams(ZATOM=29., LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,1]], RMAX=7, GMAX=65, NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1)
            p.set_value('LDAU_PARA', [1,2])
            p._check_input_consistency()
        except TypeError:
           knownError = True
        assert knownError

    def test_inconsistency_bulk_mode_bravais(self):
        p = kkrparams(LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,0]], NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1, RMAX=7, GMAX=65, ZATOM=29.)
        knownError = False
        try:
            p.fill_keywords_to_inputfile()
        except ValueError:
            knownError = True
        assert knownError
    
class Test_get_info():
    def test_get_mandatory(self):
        p = kkrparams()
        manlist = p.get_all_mandatory()
        assert set(manlist)==set(['LMAX', 'NAEZ', 'BRAVAIS', 'RMAX', 'GMAX', 'NSPIN', '<RBASIS>', 'ALATBASIS', '<ZATOM>'])
        
    def test_get_set_values(self):
        p = kkrparams()
        setlist = p.get_set_values()
        assert setlist==[]
        
    def test_get_set_values2(self):
        from numpy import array
        p = kkrparams()
        p.set_multiple_values(EMIN=1, EMAX=2)
        setlist = p.get_set_values()
        assert set(array(setlist).flatten()) == set(array([['EMIN', 1.], ['EMAX', 2.]]).flatten())
        
    def test_get_description(self):
        p = kkrparams()
        desc = p.get_description('EMIN')
        assert desc=='Accuracy, Valence energy contour: Lower value (in Ryd) for the energy contour'
        
    def test_get_type(self):
        p = kkrparams()
        tlist = p.get_type('BRAVAIS')
        assert tlist == [float, float, float, float, float, float, float, float, float]
        
    def test_is_mandatory(self):
        p = kkrparams()
        man = p.is_mandatory('EMAX')
        assert (not man)

    
class Test_fill_inputfile():
    def test_fill_inputfile_minimal_Voronoi(self):
        p = kkrparams(ZATOM=29., LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,1]], RCLUSTZ=1.5, NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1)
        p.fill_keywords_to_inputfile(is_voro_calc=True)
        txt = open('inputcard').readlines()
        ref = ['ALATBASIS= 1.000000\n', 'BRAVAIS\n', '1.000000 0.000000 0.000000\n',
               '0.000000 1.000000 0.000000\n', '0.000000 0.000000 1.000000\n', 'NAEZ= 1\n', '<RBASIS>\n',
               '0.000000 0.000000 0.000000\n', '<ZATOM>\n', '29.000000\n', 'NSPIN= 2\n', 'LMAX= 2\n', 'RCLUSTZ= 1.500000\n']
        assert set(txt)==set(ref)
        
    def test_fill_inputfile_KKR(self):
        reffile = ['ALATBASIS= 1.000000\n', 'BRAVAIS\n', '1.000000 0.000000 0.000000\n', '<ZATOM>\n', '29.000000\n',
                   '0.000000 1.000000 0.000000\n', '0.000000 0.000000 1.000000\n', 'NAEZ= 1\n',
                   '<RBASIS>\n', '0.000000 0.000000 0.000000\n', 'NSPIN= 2\n', 'LMAX= 2\n', 
                   'RCLUSTZ= 1.500000\n', 'RMAX= 7.000000\n', 'GMAX= 65.000000\n']
        p = kkrparams(ZATOM=29., LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,1]], RMAX=7, GMAX=65, RCLUSTZ=1.5, NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1)
        p.fill_keywords_to_inputfile()
        txt = open('inputcard').readlines()
        assert set(txt)==set(reffile)
        
    def test_fill_inputfile_empty_check(self):
        p = kkrparams(LMAX=2, NAEZ=1)
        known_error = False
        try:
            p.fill_keywords_to_inputfile()
        except ValueError:
            known_error = True
        assert known_error
        
    def test_fill_inputfile_all_keys(self):  
        """Example filling all keys"""      
        from numpy import array, sqrt
       
        alat=5.416871386
        naez=6
        bravais=array([[0.7071067812, -0.5, 0.0], [0.7071067812, 0.5, 0.0], [sqrt(2), 0.0, 0.866025404]])
        lmax=2
        nspin=2
        nucl_numbers=[0,0,26,27,26,27,0,0]
        cpa_info = [naez+2, [1., 1., 0.98, 0.02, 0.98, 0.02, 1., 1.], [1, 2, 3, 3, 4, 4, 5, 6]]
        npol=4
        npt1, npt2, npt3 = 3, 10, 3
        tempr= 800
        basis_vectors = []
        for iatom in range(naez):
            tmp = array([0, 0, 0])+iatom*array([0.5, 0.5, bravais[2, 2]])
            tmp[0] = tmp[0]%1
            tmp[1] = tmp[1]%1
            print(iatom, tmp)
            basis_vectors.append(tmp)
        basis_vectors = array(basis_vectors)
        natyp = cpa_info[0]
        cpa_conc = cpa_info[1]
        cpa_sites = cpa_info[2]
        ins=1
        kshape=ins
        rmax, gmax= 7, 65
        rcls=1.5
        bzdivide=[10,10,0]
        emin=-0.4
        p = kkrparams()
        p.set_multiple_values(ZATOM=nucl_numbers, RBASIS=basis_vectors, BRAVAIS=bravais, NAEZ=naez, ALATBASIS=alat)
        p.set_multiple_values(NSPIN=nspin, LMAX=lmax, NPOL=npol, NPT1=npt1, NPT2=npt2, NPT3=npt3, TEMPR=tempr)
        p.set_multiple_values(RMAX=rmax, GMAX=gmax)
        p.set_multiple_values(RCLUSTZ=rcls, BZDIVIDE=bzdivide, EMIN=emin)
        p.set_multiple_values(INS=ins, KSHAPE=kshape)
        p.set_multiple_values(INTERFACE=True, NLBASIS=1, NRBASIS=1,
                              ZPERIODL=array([-0.5, -0.5, -bravais[2, 2]]), 
                              ZPERIODR=array([0.5, 0.5, bravais[2, 2]]), 
                              RBLEFT=basis_vectors[0]+array([-0.5, -0.5, -bravais[2, 2]]), 
                              RBRIGHT=basis_vectors[naez-1]+array([0.5, 0.5, bravais[2, 2]]))
        p.set_value('LINIPOL', True)
        p.set_value('XINIPOL', [1 for i in range(natyp)])
        p.set_value('HFIELD', 0.02)
        p.set_value('NSTEPS', 1)
        p.set_value('IMIX', 0)
        p.set_value('STRMIX', 0.01)
        p.set_value('CARTESIAN', False)
        p.set_multiple_values(KAOEZR=1, KAOEZL=1, FPRADIUS=[-1 for i in range(natyp)], RCLUSTXY=rcls, 
                              TKSEMI=800,EMAX=1, NPOLSEMI=0, N2SEMI=0, N1SEMI=0, N3SEMI=0, FSEMICORE=0, 
                              KVREL=1, NCHEB=7, VCONST=0, SOCSCL=[1 for i in range(natyp)], 
                              LAMBDA_XC=1, FCM=20, ITDBRY=20, KREADLDAU=0, RUNOPT=['LDAU', 'SEMICORE', 'IRGENDWAS FALSCHES'],
                              TESTOPT=['TSTOPTX0', 'TSTOPTX1', 'TSTOPTX2', 'TSTOPTX3', 'TSTOPTX4', 'TSTOPTX5', 'TSTOPTX6', 'TSTOPTX7', 'TSTOPTX8', 'TSTOPTXYZZZZZZ'], QBOUND=10**-3, 
                              NPAN_LOG=3, NPAN_EQ=4, CPAINFO=[10**-3, 20], LLOYD=0, EMUSEMI=0, ICST=2, 
                              TOLRDIF=0.01, BRYMIX=0.01, EBOTSEMI=0, NRIGHTHO=10, KEXCOR=2, NLEFTHOS=10,
                              R_LOG=0.4, LDAU_PARA=[1, 2, 0, 0, 0], NAT_LDAU=0, 
                              RMTREFL=2.3, RMTREFR=2.3, DELTAE=[10**-5, 0], RMTREF=[2.3 for i in range(natyp)])
        p.set_value('<SHAPE>', [1 for i in range(natyp)])
        p.set_multiple_values(NATYP=natyp, SITE=cpa_sites)
        p.set_value('<CPA-CONC>', cpa_conc)
        p.set_value('FILES', ['output.pot', ''])
        p.fill_keywords_to_inputfile(is_voro_calc=True)

#TODO: implement and test read_inputfile
"""
class Test_read_inputfile():
    def test_read_minimal_inputfile(self):
        p = kkrparams(LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,1]], RCLUSTZ=1.5, NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1)
        p.fill_keywords_to_inputfile(is_voro_calc=True)
        p2 = kkrparams(params_type='voronoi')
        p2.read_keywords_from_inputcard()
        assert p.values==p2.values
        assert set(p.get_all_mandatory())==set(p2.get_all_mandatory())
        
    def test_read_unsorted_inputfile(self):
        p = kkrparams(LMAX=2, NAEZ=1, BRAVAIS=[[1,0,0],[0,1,0],[0,0,1]], RCLUSTZ=1.5, NSPIN=2, RBASIS=[0,0,0], ALATBASIS=1, RMAX=7, GMAX=65)
        p.fill_keywords_to_inputfile(output='input.temp.txt')
        txt = open('input.temp.txt', 'r').readlines()
        tmp = txt[0]
        txt[0] = txt[1]; txt[1]=tmp
        open('input.temp_unsorted.txt', 'w').writelines(txt)
        p2 = kkrparams()
        p2.read_keywords_from_inputcard()
        assert p.values==p2.values
        assert set(p.get_all_mandatory())==set(p2.get_all_mandatory())
        
    def test_read_old_inputstyle(self):
        p = kkrparams()
        p.read_keywords_from_inputcard(inputcard='input_old')
        #TODO define input_old and assertion
        p2 = kkrparams(LMAX=2)
        assert p.values==p2.values
"""
        
        
