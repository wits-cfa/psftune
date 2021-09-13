#!/usr/bin/env python

import os
import sys
from setuptools import setup, find_packages


requirements = ["numpy",
                "click",
                "matplotlib",
                "astropy",
                ],

PACKAGE_NAME = "psftune"
__version__ = "0.1.0"

setup(name=PACKAGE_NAME,
      version=__version__,
      description="Compare PSFs from an array telescope based on sidelobe behaviour",
      author="Sphesihle Makhathini",
      author_email="sphemakh@gmail.com",
      url="https://github.com/wits-cfa/psftune",
      packages=find_packages(),
      include_package_data=True,
      python_requires='>=3.6',
      install_requires=requirements,
      entry_points="""
            [console_scripts]
            psftune = psftune.main:cli
      """,
      classifiers=[],
      )
