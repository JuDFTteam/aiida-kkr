# delete existing cache files and rerun creation
# requires working installation of JuKKR codes
# the list of tests that are relevant is taken from the GITHUB_SUITE list in run_all.sh

# first remove old files
rm data_dir/kkr_dos_wc-nodes-* \
   data_dir/kkr_flex_wc-nodes-* \
   data_dir/kkr_scf_wc-nodes-* \
   data_dir/kkr_eos_wc-nodes-* \
   data_dir/kkr_imp_sub_wc-nodes-* \
   data_dir/kkr_startpot_wc-nodes-*

# then rerun the tests
pytest -svx ./workflows/test_vorostart_wc.py \
            ./workflows/test_dos_wc.py \
            ./workflows/test_gf_writeout_wc.py \
            ./workflows/test_scf_wc_simple.py \
            ./workflows/test_eos.py \
            ./workflows/test_kkrimp_sub_wc.py
