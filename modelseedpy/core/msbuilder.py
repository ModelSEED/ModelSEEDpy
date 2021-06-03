import logging

import re
import copy
from optlang.symbolics import Zero, add
from cobra.core import Gene, Metabolite, Model, Reaction
from modelseedpy.core.msgenome import MSGenome

logger = logging.getLogger(__name__)

#Static functions
def build_model():
    