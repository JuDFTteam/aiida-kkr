#!/usr/bin/env sh

pytest -v -s --cov-report=html --cov=aiida_kkr --ignore=jukkr --mpl -p no:warnings
#pytest -v -s --cov-report=html --cov=aiida_kkr --ignore=jukkr --ignore=workflows --mpl -p no:warnings
