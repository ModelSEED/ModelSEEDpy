import logging

import re
import copy
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)

>b0001|TERMKEY1:wathahshf @ asdasdas @ asdasds @ asdsdsd @ sdsds|TERMKEY2:wathahshf

input(KBASE) (INPUT?)  ??? >>>>MSGenome >>> (KBASE) (FASTA ?)????

build_model_app(genome, template)
build_model_rast_app(genome, template)
build_model_ko_app(genome, template)
build_model_ec_app(genome, template)

class Sequence:

    def __init__(self, feature_id, sequence, description=None):
        self.id = feature_id
        self.seq = sequence
        self.description = description
        #self.ontology_terms = {'termkey' : []} :


class MSGenome:

    def __init__(self):
        self.features = DictList()
