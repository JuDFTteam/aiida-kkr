#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';
pytest -vs $@
