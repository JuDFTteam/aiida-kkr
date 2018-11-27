[![Documentation Status](https://readthedocs.org/projects/aiida-kkr/badge/?version=latest)](https://aiida-kkr.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/JuDFTteam/aiida-kkr.svg?branch=develop)](https://travis-ci.org/JuDFTteam/aiida-kkr)
[![codecov](https://codecov.io/gh/JuDFTteam/aiida-kkr/branch/develop/graph/badge.svg)](https://codecov.io/gh/JuDFTteam/aiida-kkr)
[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)


# aiida-kkr

AiiDA plugin for the KKR codes plus workflows and utility.


# Installation

```shell
$ git clone https://github.com/broeder-j/aiida-kkr

$ cd aiida-kkr
$ pip install -e .  # also installs aiida, if missing (but not postgres)
$ reentry scan -r aiida  
$ verdi quicksetup  # better to set up a new profile
$ verdi calculation plugins  # should now show kkr.* entrypoints
```

# Usage and Documentation

See http://aiida-kkr.readthedocs.io for user's guide and API reference.
