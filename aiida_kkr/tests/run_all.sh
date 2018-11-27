#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
#pytest -sv
#pytest -v
#pytest --cov-report=term-missing --cov=aiida_kkr --ignore=test_scf_wc_simple.py
pytest --cov-report=term-missing --cov=aiida_kkr --ignore=test_scf_wc_simple.py --ignore=test_vorostart_wc.py --ignore=test_dos_wc.py --ignore=test_gf_writeout_wc.py --ignore=test_kkrimp_sub_wc.py --ignore=test_kkrimp_full_wc.py
