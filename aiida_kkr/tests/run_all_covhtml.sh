#!/usr/bin/env sh
#export AIIDA_PATH='.';
#mkdir -p '.aiida';
#
##
#AIIDA_BACKEND="django"
#AIIDA_EMAIL="aiida@localhost"
#AIIDA_DB_HOST="ruess"
#AIIDA_DB_PORT="5432"
#AIIDA_DB_NAME="aiida_default"
#AIIDA_DB_USER="aiida"
#AIIDA_DB_PASS="password"
#AIIDA_REPO="/Users/ruess/sourcecodes/aiida/repositories/scratch_local_machine/"
#AIIDA_FIRST_NAME="AiiDa"
#AIIDA_LAST_NAME="In Docker"
#AIIDA_INSTITUTION="None"
##
#verdi setup \
#        --non-interactive \
#        --backend="${AIIDA_BACKEND}" \
#        --email="${AIIDA_EMAIL}" \
#        --db_host="${AIIDA_DB_HOST}" \
#        --db_port="${AIIDA_DB_PORT}" \
#        --db_name="${AIIDA_DB_NAME}" \
#        --db_user="${AIIDA_DB_USER}" \
#        --db_pass="${AIIDA_DB_PASS}" \
#        --repo="${AIIDA_REPO}" \
#        --first-name="${AIIDA_FIRST_NAME}" \
#        --last-name="${AIIDA_LAST_NAME}" \
#        --institution="${AIIDA_INSTITUTION}" \
#        --no-password
#
##
#verdi computer setup < "computer_setup.txt" 
#verdi computer configure slurmcontrol
##
## Configure kkr code
#verdi code setup < "code_setup_kkr.txt"
#verdi code setup < "code_setup_voronoi.txt"

#pytest -v --cov-report=html --cov=aiida_kkr --ignore=test_entrypoints.py --ignore=test_scf_wc_simple.py --ignore=test_common_workfunctions.py
#pytest -v --cov-report=html --cov=aiida_kkr --ignore=test_scf_wc_simple.py
#pytest -s --cov-report=html --cov=aiida_kkr --ignore=test_scf_wc_simple.py 
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_vorostart_workflow
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_kkrimp_full_workflow
#pytest -s --cov-report=html --cov=aiida_kkr -k Test_scf_workflow

pytest --cov-report=html --cov=aiida_kkr
