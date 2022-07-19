'''
This small dirty scripts helps to migrate aiida exports files produced by tests
prob works only under linux.

Modified from from AiiDA-Fleur: https://github.com/JuDFTteam/aiida-fleur/blob/6ab4a45a579fc7ebc73ebc5b4b3ab1b5bcb26dc8/tests/migrate_exports.py
'''
import os
import shutil
import subprocess

data_dirs = ['data_dir/']

for folder in data_dirs:
    for infile in os.listdir(folder):
        print('migrating aiida export file: ' + folder + infile)
        subprocess.run(['verdi', 'archive', 'migrate', folder + infile, '--in-place'], check=True)
