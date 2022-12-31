# -*- coding: utf-8 -*-
project = "ModelSEEDpy"
copyright = "2022, DOE KBase"
author = "Filipe Liu, Andrew P. Freiburger, Chris Henry"

master_doc = 'index'
release = '1'
version = '0.0.1'

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

extensions = [
#     "sphinx.ext.autodoc",
#     "sphinx.ext.intersphinx",
#     "sphinx.ext.mathjax",
#     "sphinx.ext.viewcode",
#     "sphinx.ext.napoleon",
#     "sphinx.ext.autosummary",
    "autoapi",
    "nbsphinx",
]


# import sys
# from os.path import dirname, join
# SRC_PATH = join(dirname(dirname(__file__)), "src")
# sys.path.insert(0, SRC_PATH)
# autoapi_dirs = [SRC_PATH]
