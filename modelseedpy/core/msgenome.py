import logging

import re
import copy
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)


def normalize_role(s):
    #print(s)
    s = s.strip().lower()
    s = re.sub('[\W_]+', '', s)
    return s

#Static factory functions:
            
#def build_from_kbase_gto:


def read_fasta(f, split='|', h_func=None):
    with open(f, 'r') as fh:
        return parse_fasta_str(fh.read(), split, h_func)


def parse_fasta_str(faa_str, split='|', h_func=None):
    lines = faa_str.split('\n')
    features = []
    seq = None
    for line in lines:
        if line.startswith('>'):
            if seq:
                features.append(seq)
            seq_id = line[1:]
            desc = None
            if h_func:
                seq_id, desc = h_func(seq_id)
            elif split:
                header_data = line[1:].split(split, 1)
                seq_id = header_data[0]
                desc = header_data[1]

            seq = MSFeature(seq_id, "", desc)
        else:
            if seq:
                seq.seq += line.strip()
    if seq and seq.seq and len(seq.seq) > 0:
        features.append(seq)
    return features


class MSFeature:

    def __init__(self, feature_id, sequence, description=None):
        self.id = feature_id
        self.seq = sequence
        self.description = description  # temporary replace with proper parsing
        self.ontology_terms = {}

    def add_ontology_term(self, ontology, value):
        if ontology not in self.ontology_terms:
            self.ontology_terms[ontology] = []
        if value not in self.ontology_terms[ontology]:
            self.ontology_terms[ontology].append(value)


class MSGenome:

    def __init__(self):
        self.features = DictList()

    @staticmethod
    def from_fasta(filename, contigs=0, split='|', h_func=None):
        genome = MSGenome()
        genome.features += read_fasta(filename, split, h_func)
        return genome

    @staticmethod
    def from_dna_fasta(filename):
        pass

    @staticmethod
    def from_protein_sequences_hash(sequences):
        """

        """
        features = []
        for seq_id in sequences:
            features.append(MSFeature(seq_id, sequences[seq_id]))
        genome = MSGenome()
        genome.features += features
        return genome
