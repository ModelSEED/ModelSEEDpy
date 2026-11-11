# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.rst") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

setup(
    name="ModelSEEDpy",
    version="0.4.3",
    description="Python package for building and analyzing models using ModelSEED",
    long_description_content_type="text/x-rst",
    long_description=readme,
    author="Christopher Henry",
    author_email="chenry@anl.gov",
    url="https://github.com/ModelSEED/ModelSEEDpy",
    license=license,
    packages=find_packages(exclude=("docs")),
    package_data={
        "modelseedpy": ["config.cfg", "data/*", "data/templates/*", "data/ml/*"],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Natural Language :: English",
    ],
    install_requires=[
        "networkx >= 2.4",
        "cobra >= 0.29.0",
        "pandas >= 2.2.2",
        "numpy >= 1.24",
        # the default genome classifier (data/ml/*.npz) needs only numpy; scikit-learn is used by
        # modelseedpy.ml. Legacy pickle classifiers (*_01_17_2023) still need scikit-learn 1.2.0.
        "scikit-learn >= 1.2.0",
        "scipy >= 1.5.4",
        "chemicals >= 1.0.13",
        "chemw >= 0.3.2",
        "matplotlib >= 3.0.0",
        "Jinja2 >= 3.1.4",
        "sympy >=1.12.0",
    ],
    tests_require=[
        "pytest",
    ],
    project_urls={
        "Documentation": "https://modelseedpy.readthedocs.io/en/latest/",
        "Issues": "https://github.com/ModelSEED/ModelSEEDpy/issues",
    },
)
