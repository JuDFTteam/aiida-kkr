#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
#pytest -sv
#pytest -v
pytest --cov-report=html --cov=aiida_kkr
#pytest --cov-report=term-missing --cov=aiida_kkr
