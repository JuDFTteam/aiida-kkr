# -*- coding: utf-8 -*-
"""
A simple script that checks the consistency between the version number specified in
setup.json, and the version in the __init__.py file.
"""

import json
import os
import re
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.join(SCRIPT_DIR, os.path.pardir)

# Get the __init__.py version number
with open(os.path.join(ROOT_DIR, 'aiida_kkr/__init__.py')) as f:
    MATCH_EXPR = "__version__[^'\"]+(['\"])([^'\"]+)"
    VERSION_INIT = re.search(MATCH_EXPR, f.read()).group(2).strip()  # type: ignore

# Get the setup.py version number (from the setup.json)
with open(os.path.join(ROOT_DIR, 'setup.json')) as f:
    VERSION_SETUP = json.load(f)['version']

# Get the pyproject.toml version number
with open(os.path.join(ROOT_DIR, 'pyproject.toml')) as f:
    MATCH_EXPR = "version[^'\"]+(['\"])([^'\"]+)"
    VERSION_PYPROJECT = re.search(MATCH_EXPR, f.read()).group(2).strip()  # type: ignore

# Get the .bumpversion.cfg current_version number
with open(os.path.join(ROOT_DIR, '.bumpversion.cfg')) as f:
    MATCH_EXPR = r'current_version\s*=\s*([^\\S\r\n]+)'
    VERSION_BUMPVERSION = re.search(MATCH_EXPR, f.read()).group(1).strip()  # type: ignore

if VERSION_INIT != VERSION_SETUP:
    print("Version numbers don't match: init:'{}', setup:'{}' ".format(VERSION_INIT, VERSION_SETUP))
    sys.exit(1)

if VERSION_INIT != VERSION_PYPROJECT:
    print("Version numbers don't match: init:'{}', pyproject:'{}' ".format(VERSION_INIT, VERSION_PYPROJECT))
    sys.exit(1)

if VERSION_INIT != VERSION_BUMPVERSION:
    print("Version numbers don't match: init:'{}', bumpversion:'{}' ".format(VERSION_INIT, VERSION_BUMPVERSION))
    sys.exit(1)
