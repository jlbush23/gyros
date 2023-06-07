# -*- coding: utf-8 -*-
"""
Created on Wed Feb  3 21:44:18 2021

@author: jlbus
"""
from setuptools import setup

dependencies = ['numpy','pandas',
                'astroquery','astropy','matplotlib','scipy',
                'lightkurve',
                'tess_cpm']

setup(
      name = 'gyros',
      url = 'https://github.com/jlbush23/gyros',
      author = 'Jonathan Bush',
      author_email = 'jlbush23@gmail.com',
      packages = ['gyros'],
      install_requires = dependencies,
      version = '0.1',
      description = "Download light curves from NASA's Kepler, K2, and TESS space missions!",
      #long_description = open('README.md').read()      
     )