# -*- coding: utf-8 -*-
import logging
import re
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)

DEFAULT_SPLIT = " "


def normalize_role(s):
    s = s.strip().lower()
    s = re.sub(r"[\W_]+", "", s)
    return s


def read_fasta(f, split=DEFAULT_SPLIT, h_func=None):
    with open(f, "r") as fh:
        return parse_fasta_str(fh.read(), split, h_func)


def parse_fasta_str(faa_str, split=DEFAULT_SPLIT, h_func=None):
    features = []
    seq = None
    for line in faa_str.split("\n"):
        if line.startswith(">"):
            if seq:
                features.append(seq)
            desc = None
            seq_id = line[1:]
            if h_func:
                seq_id, desc = h_func(line[1:])
            elif split:
                header_data = line[1:].split(split, 1)
                seq_id = header_data[0]
                if len(header_data) > 1:
                    desc = header_data[
                        1
                    ]  # The unit test throws an error when this is commented


            seq = MSFeature(seq_id, "", desc)
        else:
            if seq:
                seq.seq += line.strip()
    if seq and seq.seq and len(seq.seq) > 0:
        features.append(seq)
    return features


class MSFeature:
    def __init__(self, feature_id, sequence, description=None):
        """

        @param feature_id: identifier for the protein coding feature
        @param sequence: protein sequence
        @param description: description of the feature
        """

        self.id = feature_id
        self.seq = sequence
        self.description = description  # temporary replace with proper parsing
        self.ontology_terms = {}
        self.aliases = []

    def add_ontology_term(self, ontology_term, value):
        """
        Add functional term to the feature

        @param ontology_term: type of the ontology (e.g., RAST, EC)
        @param value: value for the ontology (e.g., pyruvate kinase)
        """
        if ontology_term not in self.ontology_terms:
            self.ontology_terms[ontology_term] = []
        if value not in self.ontology_terms[ontology_term]:
            self.ontology_terms[ontology_term].append(value)


class MSGenome:
    def __init__(self):
        self.features = DictList()

    def add_features(self, feature_list: list):
        """

        :param feature_list:
        :return:
        """
        duplicates = list(filter(lambda o: o.id in self.features, feature_list))
        if len(duplicates) > 0:
            raise ValueError(
                f"unable to add features {duplicates} already present in the genome"
            )

        for f in feature_list:
            f._genome = self

        self.features += feature_list

    @staticmethod
    def from_fasta(
        filename, contigs=0, split="|", h_func=None
    ):  # !!! the contigs argument is never used
        genome = MSGenome()
        genome.features += read_fasta(filename, split, h_func)
        return genome

    def to_fasta(self, filename, l=80, fn_header=None):
        with open(filename, "w") as fh:
            for feature in self.features:
                h = f">{feature.id}\n"
                if fn_header:
                    h = fn_header(feature)
                fh.write(h)
                lines = [
                    feature.seq[i : i + l] + "\n" for i in range(0, len(feature.seq), l)
                ]
                for line in lines:
                    fh.write(line)
        return filename

    @staticmethod
    def from_dna_fasta(filename):
        pass

    @staticmethod
    def from_protein_sequences_hash(sequences):
        features = [MSFeature(seq_id, sequences[seq_id]) for seq_id in sequences]
        genome = MSGenome()
        genome.features += features
        return genome

    def alias_hash(self):
        return {alias: gene for gene in self.features for alias in gene.aliases}

    def search_for_gene(self, query):
        if query in self.features:
            return self.features.get_by_id(query)
        aliases = self.alias_hash()
        return aliases[query] if query in aliases else None
