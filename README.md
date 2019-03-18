[![Documentation Status](https://readthedocs.org/projects/aiida-kkr/badge/?version=latest)](https://aiida-kkr.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/JuDFTteam/aiida-kkr.svg?branch=master)](https://travis-ci.org/JuDFTteam/aiida-kkr)
[![codecov](https://codecov.io/gh/JuDFTteam/aiida-kkr/branch/master/graph/badge.svg)](https://codecov.io/gh/JuDFTteam/aiida-kkr)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)
[![GitHub version](https://badge.fury.io/gh/JuDFTteam%2Faiida-kkr.svg)](https://badge.fury.io/gh/JuDFTteam%2Faiida-kkr)
[![PyPI version](https://badge.fury.io/py/aiida-kkr.svg)](https://badge.fury.io/py/aiida-kkr)


# aiida-kkr

AiiDA plugin for the KKR codes plus workflows and utility.

## Features

* KKR calculations for bulk and interfaces
* treatment of alloys using VCA or CPA
* self-consistency, DOS and bandstructure calculations
* extraction of magnetic exchange coupling parameters (*J_ij*, *D_ij*)
* impurity embedding solving the Dyson equation
* ~~import old calculations using the calculation importer~~ (only working with aiida-core<1.0, i.e. in aiida-kkr v0.1.2)


# Installation

```shell
$ pip install aiida-kkr  # install latest version of aiida-kkr (published on pypi.org)
$ reentry scan -r aiida  # update entry points, needed in order to find kkr.* entrypoints in aiida

# setupt aiida if this was not done already:
$ verdi quicksetup  # better to set up a new profile
$ verdi calculation plugins  # should now show kkr.* entrypoints
```

for developer version download the repository and install the downloaded version
```shell
$ git clone https://github.com/JuDFTteam/aiida-kkr.git
$ pip install -e aiida-kkr
$ reentry scan -r aiida
```

# Usage and Documentation

See http://aiida-kkr.readthedocs.io for user's guide and API reference.

# Contribting guide

* Under construction ...
* ...

