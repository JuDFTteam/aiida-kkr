#!/usr/bin/env sh
export AIIDA_PATH='.';
mkdir -p '.aiida';

#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr --mpl -p no:warnings

# tests without running actual calculations
pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=workflows --ignore=jukkr --mpl -p no:warnings

# test running full workflows, need compiled codes and execute them
pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_vorostart_workflow

# these tests fail at the moment because gfortran is too old
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_dos_workflow
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_gf_writeout_workflow
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_scf_workflow
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_eos_workflow

# these are deactivated becaus gfortran compilation of kkrimp does not work
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_kkrimp_scf_workflow
#pytest --cov-report=term-missing --cov-append --cov=aiida_kkr --ignore=jukkr -k Test_kkrimp_full_workflow
