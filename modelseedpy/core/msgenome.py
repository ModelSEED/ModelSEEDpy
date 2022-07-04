import logging
logger = logging.getLogger(__name__)
from cobra.core.dictlist import DictList
import re


def normalize_role(s):
    # print(s)
    s = s.strip().lower()
    s = re.sub('[\W_]+', '', s)
    return s

def read_fasta(fasta, split='|', h_func=None):
    with open(fasta, 'r') as fh:
        return parse_fasta_str(fh.read(), split, h_func)

def parse_fasta_str(faa_str, split='|', h_func=None):
    features = []
    seq = None
    for line in faa_str.split('\n'):
        if line.startswith('>'):
            if seq:
                features.append(seq)
            desc = None
            seq_id = line[1:]
            if h_func:
                seq_id, desc = h_func(line[1:])
            elif split:
                header_data = line[1:].split(split, 1)
                seq_id = header_data[0]
                desc = header_data[1]  # The unit test throws an error when this is commented


            seq = MSFeature(seq_id, "", desc)
        else:
            if seq:
                seq.seq += line.strip()
    if seq and seq.seq and len(seq.seq) > 0:
        features.append(seq)
    return features


class MSFeature:
    def __init__(self, feature_id, sequence, description=None):
        self.id = feature_id; self.seq = sequence
        self.description = description  # temporary replace with proper parsing
        self.ontology_terms = {}
        self.aliases = []

    def add_ontology_term(self, ontology_term, value):
        if ontology_term not in self.ontology_terms:
            self.ontology_terms[ontology_term] = []
        if value not in self.ontology_terms[ontology_term]:
            self.ontology_terms[ontology_term].append(value)


class MSGenome:
    def __init__(self):
        self.features = DictList()

    @staticmethod
    def from_fasta(filename, split='|', h_func=None):
        genome = MSGenome()
        genome.features += read_fasta(filename, split, h_func)
        return genome

    # @staticmethod
    # def from_dna_fasta(filename):
    #     pass

    @staticmethod
    def from_protein_sequences_hash(sequences):
        features = [MSFeature(seq_id, sequences[seq_id]) for seq_id in sequences]
        genome = MSGenome()
        genome.features += features
        return genome
    
    def alias_hash(self):
        return {alias:gene for gene in self.features for alias in gene.aliases}
    
    def search_for_gene(self,query):
        if query in self.features:
            return self.features.get_by_id(query)
        aliases = self.alias_hash()
        return aliases[query] if query in aliases else None