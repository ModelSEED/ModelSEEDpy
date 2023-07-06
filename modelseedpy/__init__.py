# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import sys
import logging
import cobra
import warnings as _warnings
from os import name as _name
from os.path import abspath as _abspath
from os.path import dirname as _dirname
from modelseedpy.helpers import config

__author__ = "Christopher Henry"
__email__ = "chenry@anl.gov"
__version__ = "0.2.2"

logger = logging.getLogger(__name__)

print("modelseedpy", __version__)

if sys.version_info[0] == 2:
    logger.warning(
        "Python 2 is reaching end of life (see "
        "https://www.python.org/dev/peps/pep-0373/) and many cobra "
        "dependencies have already dropped support. At the moment it *should* "
        "still work but we will no longer actively maintain Python 2 support."
    )

if "e0" not in cobra.medium.annotations.compartment_shortlist["e"]:
    cobra.medium.annotations.compartment_shortlist["e"].append("e0")

import modelseedpy
from modelseedpy.core import (
    RastClient,
    MSGenome,
    MSBuilder,
    MSMedia,
    MSGrowthPhenotypes,
    MSModelUtil,
    FBAHelper,
    MSEditorAPI,
    MSATPCorrection,
    MSGapfill,
    MSEquation,
)
from modelseedpy.core.exceptions import *

from modelseedpy.community import MSCommunity, MSCompatibility, CommKineticPkg

from modelseedpy.biochem import ModelSEEDBiochem

from modelseedpy.fbapkg import (
    BaseFBAPkg,
    RevBinPkg,
    ReactionUsePkg,
    SimpleThermoPkg,
    TotalFluxPkg,
    BilevelPkg,
    KBaseMediaPkg,
    FluxFittingPkg,
    ProteomeFittingPkg,
    GapfillingPkg,
    MetaboFBAPkg,
    FlexibleBiomassPkg,
    ProblemReplicationPkg,
    FullThermoPkg,
    MSPackageManager,
    ObjConstPkg,
    ChangeOptPkg,
    ElementUptakePkg,
)

from modelseedpy.multiomics import MSExpression
