#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
#pytest -sv
#pytest -v
pytest --cov-report=term-missing --cov=aiida_kkr --ignore=test_scf_wc_simple.py
