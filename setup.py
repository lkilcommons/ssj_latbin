# (C) 2021 University of Colorado AES-CCAR-SEDA (Space Environment Data Analysis) Group
# Written by Liam Kilcommons (Mar 2021)

import os
import glob

os.environ['DISTUTILS_DEBUG'] = "1"

from setuptools import setup, Extension
from setuptools.command import install as _install

setup(name='ssjlatbin',
      version = "0.1.1",
      description = "Bin SSJ data by latitude for ML purposes",
      author = "Liam Kilcommons",
      author_email = 'liam.kilcommons@colorado.edu',
      url = "https://github.com/lkilcommons/ssj_latbin",
      download_url = "https://github.com/lkilcommons/ssj_latbin",
      long_description = ("Code used to create latitude binned SSJ particle fluxes"),
      install_requires=['numpy','netCDF4','matplotlib','pandas','pyarrow','logbook','cdflib','toml','geospacepy'],
      packages=['ssjlatbin'],
      package_dir={'ssjlatbin' : 'ssjlatbin'},
      license='LICENSE.txt',
      zip_safe = False,
      classifiers = [
            "Development Status :: 4 - Beta",
            "Topic :: Scientific/Engineering",
            "Intended Audience :: Science/Research",
            "Natural Language :: English",
            "Programming Language :: Python"
            ]
      )