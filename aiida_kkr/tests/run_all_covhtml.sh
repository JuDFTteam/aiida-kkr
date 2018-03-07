#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
pytest -v --cov-report=html --cov=aiida_kkr
