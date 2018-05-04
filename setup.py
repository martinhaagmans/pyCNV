#!/usr/bin/env python

from distutils.core import setup

setup(name='pyCNV',
      version='1.0',
      description='Python Copy Number Variant detector for AMC-NGS-pipeline',
      author='Martin Haagmans',
      author_email='mahaagmans@gmail.com',
      packages=['pycnv'],
      url='www.github.com/zaag/pyCNV',
      license='MIT',
      scripts=['pyCNV'],
      install_requires=['numpy', 'matplotlib', 'seaborn', 'pandas'],
      )
