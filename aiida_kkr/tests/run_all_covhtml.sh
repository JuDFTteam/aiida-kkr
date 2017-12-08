#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
pytest --cov-report=html --cov=aiida_kkr
