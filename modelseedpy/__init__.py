# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import sys
import logging
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname
from modelseedpy.helpers import config

logging_hash = {
    "debug":logging.DEBUG,
    "critical":logging.CRITICAL,
    "error":logging.ERROR,
    "warning":logging.WARNING,
    "info":logging.INFO
}

#Configuing modelseedpy logger
logger = logging.getLogger(__name__)
c_handler = logging.StreamHandler()
c_handler.setLevel(logging_hash[config.get("logging","console_level")])
c_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
c_handler.setFormatter(c_format)
logger.addHandler(c_handler)
if config.get("logging","log_file") == "yes":
    f_handler = logging.FileHandler(config.get("logging","filename"),mode="w")
    f_handler.setLevel(logging_hash[config.get("logging","file_level")])
    f_format = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    f_handler.setFormatter(f_format)
    logger.addHandler(f_handler)

if sys.version_info[0] == 2:
    logger.warning("Python 2 is reaching end of life (see "
        "https://www.python.org/dev/peps/pep-0373/) and many cobra "
        "dependencies have already dropped support. At the moment it *should* "
        "still work but we will no longer actively maintain Python 2 support.")

import modelseedpy
from modelseedpy.core import (
    RastClient, MSGenome, MSBuilder, MSMedia, MSGrowthPhenotypes,MSModelUtil,
    FBAHelper, MSEditorAPI, MSATPCorrection, MSGapfill, MSEquation
)
from modelseedpy.core.exceptions import *

from modelseedpy.community import (MSCommunity, MSCompatibility, CommKineticPkg)

from modelseedpy.fbapkg import (
    BaseFBAPkg, RevBinPkg, ReactionUsePkg, SimpleThermoPkg, TotalFluxPkg, BilevelPkg,
    KBaseMediaPkg, FluxFittingPkg, ProteomeFittingPkg, GapfillingPkg, MetaboFBAPkg, FlexibleBiomassPkg,
    ProblemReplicationPkg, FullThermoPkg, MSPackageManager, ObjConstPkg, ChangeOptPkg, ElementUptakePkg
)

from modelseedpy.multiomics import MSExpression

__version__ = "0.2.2"
