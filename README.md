[![Documentation Status](https://readthedocs.org/projects/aiida-kkr/badge/?version=latest)](https://aiida-kkr.readthedocs.io/en/latest/?badge=latest)

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
