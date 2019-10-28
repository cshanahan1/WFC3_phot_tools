#!/usr/bin/env python

from setuptools import find_packages
from setuptools import setup

setup(name = 'standard_star_phot_tools',
      description = 'Tools for photometry of calibration stars',
      author = 'Space Telescope Science Institute',
      url = 'https://github.com/shannnnnyyy/standard_star_photometry',
      packages = find_packages(),
      install_requires = ['astropy', 'matplotlib', 'numpy', 'photutils'])