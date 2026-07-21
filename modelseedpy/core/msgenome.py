# -*- coding: utf-8 -*-
import logging
import re
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)

DEFAULT_SPLIT = " "


def to_fasta(features, filename, line_size=80, fn_header=None):
    with open(filename, "w") as fh:
        for feature in features:
            if feature.seq:
                h = f">{feature.id}\n"
                if fn_header:
                    h = fn_header(feature)
                fh.write(h)
                _seq = feature.seq
                lines = [
                    _seq[i : i + line_size] + "\n"
                    for i in range(0, len(_seq), line_size)
                ]
                for line in lines:
                    fh.write(line)
    return filename


def normalize_role(s):
    s = s.strip().lower()
    s = re.sub(r"[\W_]+", "", s)
    return s


def read_fasta(f, split=DEFAULT_SPLIT, h_func=None):
    if f.endswith(".gz"):
        import gzip

        with gzip.open(f, "rb") as fh:
            return parse_fasta_str(fh.read().decode("utf-8"), split, h_func)
    else:
        with open(f, "r") as fh:
            return parse_fasta_str(fh.read(), split, h_func)


def read_fasta2(f, split=DEFAULT_SPLIT, h_func=None):
    if f.endswith(".gz"):
        import gzip

        with gzip.open(f, "rb") as fh:
            return extract_features(fh.read().decode("utf-8"), split, h_func)
    else:
        with open(f, "r") as fh:
            return extract_features(fh.read(), split, h_func)


def parse_fasta_str(faa_str, split=DEFAULT_SPLIT, h_func=None):
    features = []
    seq = None
    for line in faa_str.split("\n"):
        if line.startswith(">"):
            if seq:
                features.append(seq)
            seq_id = line[1:]
            desc = None
            if h_func:
                seq_id, desc = h_func(seq_id)
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


def read_gbff_records_from_file(filename: str):
    if filename.endswith(".gbff"):
        with open(filename, "r") as fh:
            return read_gbff_records(fh)
    elif filename.endswith(".gz"):
        import gzip
        from io import StringIO

        with gzip.open(filename, "rb") as fh:
            return read_gbff_records(StringIO(fh.read().decode("utf-8")))


def read_gbff_records(handler):
    from Bio import SeqIO

    gbff_records = []
    for record in SeqIO.parse(handler, "gb"):
        gbff_records.append(record)
    return gbff_records


def extract_features(faa_str, split=DEFAULT_SPLIT, h_func=None):
    features = []
    active_seq = None
    seq_lines = []
    for line in faa_str.split("\n"):
        if line.startswith(">"):
            if active_seq is not None:
                active_seq.seq = "".join(seq_lines)
                features.append(active_seq)
                seq_lines = []
            seq_id = line[1:]
            desc = None
            if h_func:
                seq_id, desc = h_func(seq_id)
            elif split:
                header_data = line[1:].split(split, 1)
                seq_id = header_data[0]
                if len(header_data) > 1:
                    desc = header_data[1]
            active_seq = MSFeature(seq_id, "", desc)
        else:
            seq_lines.append(line.strip())

    # add last sequence
    if len(seq_lines) > 0:
        active_seq.seq = "".join(seq_lines)
        features.append(active_seq)

    return features


class MSFeature:
    def __init__(self, feature_id, sequence, description=None, aliases=[]):
        """

        @param feature_id: identifier for the protein coding feature
        @param sequence: protein sequence
        @param description: description of the feature
        """

        self.id = feature_id
        self.seq = sequence
        self.description = description  # temporary replace with proper parsing
        self.ontology_terms = {}
        self.aliases = aliases

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
        self.id = None
        self.annoont = None
        self.scientific_name = None

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

    def create_new_feature(self,id,sequence):
        newftr = MSFeature(id,sequence)
        self.add_features([newftr])
        return newftr

    @staticmethod
    def from_annotation_ontology(
        annoont, prioritized_event_list=None, ontologies=None, merge_all=False,feature_type=None, translate_to_rast=True
    ):
        gene_hash = annoont.get_gene_term_hash()
        genome = MSGenome()
        features = []
        for gene in gene_hash:
            feature = MSFeature(gene.id,"")
            features.append(feature)
            for term in gene_hash[gene]:
                feature.add_ontology_term(term.ontology.id, term.id)
                if term.ontology.id == "SSO":
                    feature.add_ontology_term("RAST",annoont.get_term_name(term))
        genome.add_features(features)
        return genome

    @staticmethod
    def from_fasta(filename, split=" ", h_func=None):
        genome = MSGenome()
        genome.features += read_fasta2(filename, split, h_func)
        return genome

    @staticmethod
    def from_gbff_sequence(filename):
        gbff_records = read_gbff_records_from_file(filename)
        genome = MSGenome()
        features = []
        for rec in gbff_records:
            feature = MSFeature(rec.id, str(rec.seq), description=rec.description)
            features.append(feature)
        genome.features += features
        return genome

    @staticmethod
    def from_gbff_features(
        filename, feature_id_qualifier="protein_id", description_qualifier="product"
    ):
        gbff_records = read_gbff_records_from_file(filename)
        genome = MSGenome()
        features = []
        for rec in gbff_records:
            for f in rec.features:
                if f.type == "CDS":
                    translations = f.qualifiers.get("translation", [])
                    if len(translations) == 1:
                        feature_id = f.qualifiers.get(feature_id_qualifier, [None])[0]
                        description = f.qualifiers.get(description_qualifier, [None])[0]
                        if feature_id:
                            feature = MSFeature(
                                feature_id, translations[0], description=description
                            )
                            features.append(feature)
                        else:
                            logger.warning(
                                f"skip feature: unable to fetch id from qualifier {feature_id_qualifier}"
                            )
                    elif len(translations) > 1:
                        logger.warning(f"skip feature: with multiple sequences {f}")
        genome.features += features
        return genome

    def to_fasta(self, filename, l=80, fn_header=None):
        to_fasta(self.features, filename, l, fn_header)
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
        output = {}
        for gene in self.features:
            for alias in gene.aliases:
                #Check if alias is a list
                if isinstance(alias,list):
                    if alias[1] not in output:
                        output[alias[1]] = gene
                else:
                    if alias not in output:
                        output[alias] = gene
        return output

    def search_for_gene(self, query):
        if query in self.features:
            return self.features.get_by_id(query)
        aliases = self.alias_hash()
        return aliases[query] if query in aliases else None

    def _repr_html_(self):
        return f"""
        <table>
            <tr>
                <td><strong>Memory address</strong></td>
                <td>{f"{id(self):x}"}</td>
            </tr><tr>
                <td><strong>Features</strong></td>
                <td>{len(self.features)}</td>
            </tr>
        </table>"""


class GenomeGff(MSGenome):
    def __init__(self, contigs):
        self.contigs = contigs
        super().__init__()

    @staticmethod
    def read_sequence(feature_id, gff_record, expected_sequence, contigs):
        from Bio.Seq import Seq
        from Bio import Align

        protein_seq_cds = expected_sequence
        feature_contig = contigs.features.get_by_id(gff_record.contig_id)
        seq = Seq(feature_contig.seq[gff_record.start - 1 : gff_record.end])
        if gff_record.strand == "-":
            seq = seq.reverse_complement()
        seq_from_dna = str(seq.translate())
        if len(seq_from_dna) > 0 and seq_from_dna[-1] == "*":
            seq_from_dna = seq_from_dna[:-1]
        if len(protein_seq_cds) > 0 and protein_seq_cds[-1] == "*":
            protein_seq_cds = protein_seq_cds[:-1]
        eq = protein_seq_cds == seq_from_dna

        score = None
        if not eq and len(seq_from_dna) > 0:
            try:
                aligner = Align.PairwiseAligner()
                res = aligner.align(protein_seq_cds, seq_from_dna)
                score = res.score
            except ValueError as ex:
                print("error", gff_record)
                raise ex

        feature = MSFeature(feature_id, protein_seq_cds)
        feature.description = f"score: {score}"
        feature.gff = gff_record
        return feature

    @staticmethod
    def from_fna_faa_gff(
        filename_fna, filename_faa, filename_gff, _fn_get_id, prodigal=False
    ):
        genome_gff_features = _read_gff_features(filename_gff)
        genome_faa = MSGenome.from_fasta(filename_faa)
        contigs = MSGenome.from_fasta(filename_fna)

        feature_lookup = {}
        if prodigal:
            for feature in genome_faa.features:
                attr = dict(
                    x.split("=")
                    for x in feature.description.split(" # ")[-1].split(";")
                )
                if attr["ID"] not in feature_lookup:
                    feature_lookup[attr["ID"]] = feature
                else:
                    raise ValueError("")
        else:
            feature_lookup = {feature.id: feature for feature in genome_faa.features}

        features = []
        for gff_record in genome_gff_features:
            if gff_record.feature_type == "CDS":
                feature_id = gff_record.attr.get("ID")
                if _fn_get_id:
                    feature_id = _fn_get_id(gff_record)

                feature_cds = feature_lookup.get(feature_id)

                if feature_cds:
                    protein_seq_cds = feature_cds.seq
                    f = GenomeGff.read_sequence(
                        feature_id, gff_record, protein_seq_cds, contigs
                    )
                    features.append(f)
                else:
                    print(f"not found {feature_id}")

        genome = GenomeGff(contigs)
        genome.features += features
        return genome
