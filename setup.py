#!/usr/bin/env python

from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pyCNV',
      version='1.0',
      description='Python Copy Number Variant detector for AMC-NGS-pipeline',
      long_description=readme(),
      long_description_content_type='text/markdown',      
      author='Martin Haagmans',
      author_email='mahaagmans@gmail.com',
      packages=['pycnv'],
      url='https://www.github.com/zaag/pyCNV',
      license='MIT',
      scripts=['CNV'],
      install_requires=['numpy>=1.14.3', 
                        'matplotlib>=2.2.2', 
                        'seaborn>=0.8.1', 
                        'pandas>=0.23.0'],
      )




