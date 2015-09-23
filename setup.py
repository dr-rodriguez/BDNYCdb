#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup, find_packages
    setup
except ImportError:
    from distutils.core import setup
    setup

from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'DESCRIPTION.rst'), encoding='utf-8') as f: long_description = f.read()

setup(
    name='BDNYCdb',
    version='0.1.0',
    description='A Python library to communicate with the BDNYC Data Archive SQL file.',
    long_description=long_description,
    url='https://github.com/BDNYC/BDNYCdb',
    author='Joe Filippazzo',
    author_email='bdnyc.group@gmail.com',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
    ],
    keywords='astrophysics',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy','scipy','astropy'],

)