# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ModelSEEDpy',
    version='0.1.0',
    description='Python package for building and analyzing models using ModelSEED',
    long_description=readme,
    author='Christopher Henry',
    author_email='chenry@anl.gov',
    url='https://github.com/ModelSEED/ModelSEEDpy',
    license=license,
    packages=find_packages(exclude=('tests', 'docs')),
    install_requires=[
        "networkx >= 2.4",
        "cobra >= 0.17.1"
    ]
)
