#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'photometry_tools',
      description = 'Python photometry tools for WFC3 calibration',
      author = 'Space Telescope Science Institute',
      url = 'https://github.com/shannnnnyyy/photometry_tools',
      packages = find_packages(),
      install_requires = ['astropy', 'matplotlib', 'numpy', 'photutils'])