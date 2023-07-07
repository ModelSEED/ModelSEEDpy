# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.rst") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

setup(
    name="ModelSEEDpy",
    version="0.3.1",
    description="Python package for building and analyzing models using ModelSEED",
    long_description_content_type="text/x-rst",
    long_description=readme,
    author="Christopher Henry",
    author_email="chenry@anl.gov",
    url="https://github.com/ModelSEED/ModelSEEDpy",
    license=license,
    packages=find_packages(exclude=("docs")),
    package_data={
        "modelseedpy": ["config.cfg", "data/*"],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Natural Language :: English",
    ],
    install_requires=[
        "networkx >= 2.4",
        "cobra >= 0.17.1",
        "scikit-learn == 1.2.0",  # version lock for pickle ML models
        "scipy >= 1.5.4",
        "chemicals >= 1.0.13",
        "chemw >= 0.3.2",
        "matplotlib >= 3.0.0",
        "pyeda",
    ],
    tests_require=[
        "pytest",
    ],
    project_urls={
        "Documentation": "https://modelseedpy.readthedocs.io/en/latest/",
        "Issues": "https://github.com/ModelSEED/ModelSEEDpy/issues",
    },
)
