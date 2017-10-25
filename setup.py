#!/usr/bin/env python

from setuptools import setup, find_packages
import json

if __name__ == '__main__':
    # Provide static information in setup.json
    # such that it can be discovered automatically
    with open('setup.json', 'r') as info:
        kwargs = json.load(info)
    setup(
        packages=find_packages(),
        **kwargs
    )
