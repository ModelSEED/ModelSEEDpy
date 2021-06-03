import logging

import re
import copy
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)

#Static factory functions:
            
#def build_from_kbase_gto:
        
#def build_from_kbase_id:

#def build_from_protein_sequence_hash:

class MSFeature:
    def __init__(self, feature_id, sequence, description=None):
        self.id = feature_id
        self.seq = sequence
        self.description = description
        self.ontology_terms = {}


class MSGenome:

    def __init__(self):
        self.features = DictList()

    @staticmethod
    def from_fasta(filename, contigs=0):
        pass

    @staticmethod
    def from_dna_fasta(filename):
        pass

    @staticmethod
    def from_protein_sequence_hash(sequence):
        pass
