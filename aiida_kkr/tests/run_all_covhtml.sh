#!/usr/bin/env sh

#pytest -v --cov-report=html --cov=aiida_kkr --ignore=test_entrypoints.py --ignore=test_scf_wc_simple.py --ignore=test_common_workfunctions.py
#pytest -v --cov-report=html --cov=aiida_kkr --ignore=test_scf_wc_simple.py
#pytest -s --cov-report=html --cov=aiida_kkr --ignore=test_scf_wc_simple.py 
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_vorostart_workflow
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_kkrimp_full_workflow
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_scf_workflow

#pytest -s --cov-report=html --cov=aiida_kkr --ignore=jukkr -k Test_dos_workflow
#pytest -s --cov-report=html --cov=aiida_kkr --ignore=jukkr -k Test_gf_writeout_workflow -k Test_dos_workflow
#pytest --cov-report=html --cov=aiida_kkr --ignore=jukkr --ignore=test_scf_wc_simple.py --ignore=test_vorostart_wc.py --ignore=test_kkrimp_sub_wc.py --ignore=test_kkrimp_full_wc.py
pytest -v -s --cov-report=html --cov=aiida_kkr --ignore=jukkr
