# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='ModelSEEDpy',
    version='0.2.2',
    description='Python package for building and analyzing models using ModelSEED',
    long_description=readme,
    author='Christopher Henry',
    author_email='chenry@anl.gov',
    url='https://github.com/ModelSEED/ModelSEEDpy',
    license=license,
    packages=find_packages(exclude=('docs')),
    package_data={
        'modelseedpy': ['config.cfg'],
    },
    install_requires=[
        "networkx >= 2.4",
        "cobra >= 0.17.1",
        "scikit-learn == 0.23.2",  # too support KBase pickle models
        "scipy >= 1.5.4",
        "chemicals >= 1.0.13",
        "chemw >= 0.3.2",
        "matplotlib >= 3.0.0",
        "pyeda"
    ],
    project_urls={
        'Documentation': 'https://modelseedpy.readthedocs.io/en/stable/',
        'Issues': 'https://github.com/ModelSEED/ModelSEEDpy/issues',
    }
)