#!/usr/bin/env bash

export AIIDA_PATH='.';
mkdir -p '.aiida/access/test_profile'

# control test coverage by setting the envoronment variables
# SKIP_NOWORK   to run only workflow tests
# RUN_VORONOI   to run tests using voronoi
# RUN_KKRHOST   to run tests using KKRhost
# RUN_KKRIMP    to run tests using KKRimp
# RUN_ALL       to run all tests
# if the environment variable is unset or the empty string then the corresponding tests will be ignored

usage(){
  echo "$0 usage:" && grep " .)\ #" $0;
  echo
  echo "Default behavior is to run only tests that do not require running an actual calculation (i.e. skip workflow tests)."
  echo "Additional settings with the environment variables (set to something or unset to remove setting):";
  echo "  'RUN_ALL': run all tests";
  echo "  'SKIP_NOWORK': skip workflow tests";
  echo "  'RUN_VORONOI': run voronoi tests";
  echo "  'RUN_KKRHOST': run kkrhost tests";
  echo "  'RUN_KKRIMP': run kkrimp tests";
  echo "  'RUN_EOS': run eos workflow test";
  echo "  'GITHUB_SUITE': run tests as which should run on github action";
  echo
  echo "Setting that control the pytest run (also set as environment variables):";
  echo "  'HTML': create html coverage report";
  echo "  'CONTINUE_ON_ERROR': continue running tests even if we encounter errors";
  echo
  exit 0;
}

addopt=""
while getopts ":vhp:" option; do
  case $option in
    v) # add debug options for this script
       set -x # add debug output
       # Add verbosity flags '-sv' to pytest run
       addopt=" -sv " && echo "Found -v flag: adding option '$addopt' to pytest execution" ;;
    p) # run tests in parallel (specify number of cores with -p N)
       addopt=" -n $OPTARG " && echo "Found -p flag: adding option '$addopt' for parallel pytest execution using $OPTARG cores" ;;
    h) # Display help
       usage
  esac
done


# print settings before starting the tests

echo "Test settings"
echo "============="
if [[ ! -z "$HTML" ]]; then
  echo "Create html coverage report (unset 'HTML' to deactivate this)"
  repfmt="html"
else
  echo "Create XML coverage report (set 'HTML' to activate this)"
  repfmt="xml"
fi
if [[ ! -z "$CONTINUE_ON_ERROR" ]]; then
  echo "Continue running tests even if we encounter errors (unset 'CONTINUE_ON_ERROR' to deactivate this)"
else
  echo "Stop on first error (set 'CONTINUE_ON_ERROR' to deactivate this)"
  set -e # force stop on first error
fi
echo "== Test selection =="
if [[ ! -z "$GITHUB_SUITE" ]]; then
	echo "Running all test which should work on github action (unset 'GITHUB_SUITE' to prevent this)"
elif [[ ! -z "$RUN_ALL" ]]; then
  if [[ -z "$SKIP_NOWORK" ]]; then
    echo "run non-workflow tests (prevent this with 'SKIP_NOWORK' env)"
  else
    echo "skip non-workflow tests (unset 'SKIP_NOWORK' env to activate this)"
  fi
  if [[ ! -z "$RUN_VORONOI" ]]; then
    echo "run workflows using voronoi code (unset 'RUN_VORONOI' env to prevent this)"
  else
    echo "skip workflows tests using voronoi (set 'RUN_VORONOI' env to activate this)"
  fi
  if [[ ! -z "$RUN_KKRHOST" ]]; then
    echo "run workflows using KKRhost code (unset 'RUN_KKRHOST' env to prevent this)"
  else
    echo "skip workflows tests using KKRhost (set 'RUN_KKRHOST' env to activate this)"
  fi
  if [[ ! -z "$RUN_KKRIMP" ]]; then
    echo "run workflows using KKRimp code (unset 'RUN_KKRIMP' env to prevent this)"
  else
    echo "skip workflows tests using KKRimp (set 'RUN_KKRIMP' env to activate this)"
  fi
  if [[ ! -z "$RUN_EOS" ]]; then
    echo "run eos workflow (unset 'RUN_EOS' env to prevent this)"
  else
    echo "skip eos workflow tests (set 'RUN_EOS' env to activate this)"
  fi
  if [[ ! -z "$NO_RMQ" ]]; then
    echo "do not run workflows and workfuntions that need rabbitMQ (unset 'NO_RMQ' env to prevent this)"
  else
    echo "run workflows and workfunctions also with rabbitMQ (set 'NO_RMQ' env to deactivate this)"
  fi
fi
echo "============="
echo

# Now run tests

if [[ ! -z "$RUN_ALL" ]]; then
  echo "run all tests (first non-workflow tests, then workflow tests)"
  pytest --cov-report=$repfmt --cov=../ --ignore=jukkr --mpl -p no:warnings $addopt --ignore=workflows

  # now workflow tests
  pytest --cov-report=$repfmt --cov=../ --cov-append --ignore=jukkr workflows/ $addopt
elif [[ ! -z "$GITHUB_SUITE" ]]; then
  if [[ -z "$SKIP_NOWORK" ]]; then
    pytest --cov-report=$repfmt --cov=../ --cov-report xml:coverage.xml --ignore=workflows --ignore=jukkr --mpl -p no:warnings $addopt
  fi
  # only use list of currently working tests

  pytest --cov-report=$repfmt --cov-append --cov=../ -x \
      ./workflows/test_vorostart_wc.py \
	  ./workflows/test_scf_wc_simple.py \
	  ./workflows/test_dos_wc.py \
	  ./workflows/test_bs_wc.py \
	  ./workflows/test_gf_writeout_wc.py \
	  ./workflows/test_eos.py \
	  ./workflows/test_decimate.py \
	  ./workflows/test_kkrimp_dos_wc.py \
	  ./workflows/test_kkrimp_full_wc.py \
	  $addopt
  # full list, excludes:
  #  - currently tests for kkrimp_sub, jij and combine_imps are not working
  #  - tests for imp_BdG and STM are still missing
  # pytest --cov-report=$repfmt --cov-append --cov=../ -x \
     #  ./workflows/test_vorostart_wc.py \
	 #  ./workflows/test_scf_wc_simple.py \
	 #  ./workflows/test_dos_wc.py \
	 #  ./workflows/test_bs_wc.py \
	 #  ./workflows/test_gf_writeout_wc.py \
	 #  ./workflows/test_jij_wc.py \
	 #  ./workflows/test_eos.py \
	 #  ./workflows/test_decimate.py \
	 #  ./workflows/test_kkrimp_sub_wc.py \
	 #  ./workflows/test_kkrimp_dos_wc.py \
	 #  ./workflows/test_kkrimp_full_wc.py \
	 #  ./workflows/test_combine_imps.py \
	 #  $addopt
else
  # tests without running actual calculations
  if [[ -z "$SKIP_NOWORK" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run non-workflow tests"
    pytest --cov-report=$repfmt --cov=./.. --ignore=workflows --ignore=jukkr --mpl -p no:warnings $addopt
  else
    # skip things that need rabbitMQ
    if [[ -z "$SKIP_NOWORK" ]] && [[ ! -z "$NO_RMQ" ]]; then
      echo "run non-workflow tests"
      pytest --cov-report=$repfmt --cov=./.. --ignore=workflows --ignore=jukkr --ignore=calculations --ignore=test_common_workfunctions_with_rmq.py --ignore=test_plot_kkr.py --mpl -p no:warnings $addopt
    else
      echo "skipping tests that are not workflows"
    fi
  fi

  # test running full workflows, need compiled codes and execute them

  # tests using only voronoi

  if [[ ! -z "$RUN_VORONOI" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run vorostart workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_vorostart_wc.py -k test_kkr_startpot_parent_KKR $addopt

  else
    echo "skipping vorostart workflow test"
  fi

  # tests using kkrhost (and voronoi)

  if [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkr_dos workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_dos_wc.py $addopt
  else
    echo "skipping kkr_dos workflow test"
  fi
  if [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkr_bs workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_bs_wc.py $addopt
  else
    echo "skipping kkr_bs workflow test"
  fi
  if [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkr_gf_writeout workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_gf_writeout_wc.py $addopt
  else
    echo "skipping kkr_gf_writeout workflow test"
  fi
  if [[ ! -z "$RUN_VORONOI" ]] && [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run voro_start, kkr_scf workflow tests"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_vorostart_wc.py -k test_kkr_startpot_wc_Cu $addopt
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_scf_wc_simple.py $addopt
  else
    echo "skipping voro_start, kkr_scf workflow tests"
  fi
  if [[ ! -z "$RUN_VORONOI" ]] && [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]] && [[ ! -z "$RUN_EOS" ]]; then
    echo "run eos workflow tests"
    pytest --cov-report=$repfmt --cov-append --cov=./.. ./workflows/test_eos.py $addopt
  else
    echo "skipping eos workflow tests"
  fi

  # tests using kkrimp (and kkrhost/voronoi)

  if [[ ! -z "$RUN_KKRIMP" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkrimp_scf workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. --ignore=jukkr ./workflows/test_kkrimp_sub_wc.py $addopt
    #pytest --cov-report=$repfmt --cov-append --cov=./.. --ignore=jukkr ./workflows/test_kkrimp_{dos,sub}_wc.py $addopt
  else
    echo "skipping kkrimp_scf workflow test"
  fi
  if [[ ! -z "$RUN_KKRIMP" ]] && [[ ! -z "$RUN_KKRHOST" ]] && [[ ! -z "$RUN_VORONOI" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkrimp_full workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. --ignore=jukkr ./workflows/test_kkrimp_full_wc.py $addopt
  else
    echo "skipping kkrimp_full workflow test"
  fi
  if [[ ! -z "$RUN_KKRIMP" ]] && [[ ! -z "$RUN_KKRHOST" ]] && [[ -z "$NO_RMQ" ]]; then
    echo "run kkrimp_dos workflow test"
    pytest --cov-report=$repfmt --cov-append --cov=./.. --ignore=jukkr -k Test_kkrimp_dos_workflow $addopt
  else
    echo "skipping kkrimp_dos workflow test"
  fi

fi
