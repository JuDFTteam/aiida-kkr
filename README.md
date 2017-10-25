# aiida-kkr

AiiDA plugin for the KKR codes plus workflows and utility.


# Documentation

# Installation

```shell
$ git clone https://github.com/broeder-j/aiida-kkr

$ cd aiida-kkr
$ pip install -e .  # also installs aiida, if missing (but not postgres)
$ reentry scan -r aiida  
$ verdi quicksetup  # better to set up a new profile
$ verdi calculation plugins  # should now show kkr.* entrypoints
```

# Usage
