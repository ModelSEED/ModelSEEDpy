# -*- coding: utf-8 -*-
from modelseedpy.core.rpcclient import RPCClient
from modelseedpy.core.msgenome import MSFeature  # !!! import is never used
import re
from modelseedpy.core.msgenome import (
    MSGenome,
    read_fasta,
    normalize_role,
)  # move this to this lib     # !!! read_fasta is never used

### delete this after ####
def aux_rast_result(res, g):
    search_name_to_genes = {}
    search_name_to_orginal = {}
    for f in res[0]["features"]:
        function = f["function"] if "function" in f else None
        if not function:
            function = g.features.get_by_id(f["id"]).description

        if function and len(function) > 0:
            functions = re.split("; | / | @ | => ", function)
            for f_ in functions:
                f_norm = normalize_role(f_)
                if f_norm not in search_name_to_genes:
                    search_name_to_genes[f_norm] = set()
                    search_name_to_orginal[f_norm] = set()
                search_name_to_orginal[f_norm].add(f_)
                search_name_to_genes[f_norm].add(f["id"])

    return search_name_to_genes, search_name_to_orginal


### delete this after ####
def aux_template(template):
    rxn_roles = {}
    roles = dict(map(lambda x: (x["id"], x), template.roles))
    for r in template.reactions:
        rxn_roles[r.id] = set()
        complex_roles = r.get_complex_roles()
        if len(complex_roles) > 0:
            for cpx_id in complex_roles:
                for role_id in complex_roles[cpx_id]:
                    rxn_roles[r.id].add(normalize_role(roles[role_id]["name"]))
                    # print(role_id, normalize_role(roles[role_id]['name']))
    return rxn_roles


class RastClient:
    def __init__(self):
        self.rpc_client = RPCClient(
            "https://tutorial.theseed.org/services/genome_annotation"
        )
        self.stages = [
            {"name": "annotate_proteins_kmer_v2", "kmer_v2_parameters": {}},
            {
                "name": "annotate_proteins_kmer_v1",
                "kmer_v1_parameters": {"annotate_hypothetical_only": 1},
            },
            {
                "name": "annotate_proteins_similarity",
                "similarity_parameters": {"annotate_hypothetical_only": 1},
            },
        ]

    def annotate_genome(self, genome):
        p_features = []
        for f in genome.features:
            if f.seq and len(f.seq) > 0:
                p_features.append({"id": f.id, "protein_translation": f.seq})
        res = self.f(p_features)

        for o in res[0]["features"]:
            feature = genome.features.get_by_id(o["id"])
            if "function" in o:
                functions = re.split("; | / | @ | => ", o["function"])
                for function in functions:
                    feature.add_ontology_term("RAST", function)

        return res[0]["analysis_events"]

    def annotate_genome_from_fasta(self, filepath, split="|"):
        genome = MSGenome.from_fasta(filepath, split)
        res = self.annotate_genome(genome)

        return genome, res

    def annotate_protein_sequence(self, protein_id: str, protein_seq: str):
        p_features = [{"id": protein_id, "protein_translation": protein_seq}]
        return self.f(p_features)

    def annotate_protein_sequences(self, protein_seqs: dict):
        p_features = [{"id": protein_id, "protein_translation": protein_seq}]
        return self.f(p_features)

    def f1(self, protein_id, protein_seq):
        p_features = [{"id": protein_id, "protein_translation": protein_seq}]
        return self.f(p_features)

    def f(self, p_features):
        params = [{"features": p_features}, {"stages": self.stages}]
        result = self.rpc_client.call("GenomeAnnotation.run_pipeline", params)
        return result
