# -*- coding: utf-8 -*-
from modelseedpy.core.rpcclient import RPCClient
import re
from modelseedpy.core.msgenome import MSGenome, normalize_role


def split_annotation(annotation) -> list[str]:
    if annotation:
        functions = re.split("; | / | @", annotation)
        return functions
    return []


class RastClient:

    def __init__(self):
        self.rpc_client = RPCClient(
            "https://tutorial.theseed.org/services/genome_annotation"
        )
        self.stages = [
            {"name": "annotate_proteins_kmer_v2", "kmer_v2_parameters": {}},
            # {"name": "annotate_proteins_kmer_v1",
            # "kmer_v1_parameters": {"annotate_hypothetical_only": 1},},
            {
                "name": "annotate_proteins_similarity",
                "similarity_parameters": {"annotate_hypothetical_only": 1},
            },
        ]

    def annotate_genome(self, genome, split_terms=True):
        p_features = []
        for f in genome.features:
            if f.seq and len(f.seq) > 0:
                p_features.append({"id": f.id, "protein_translation": f.seq})
        res = self.f(p_features)

        if res:
            for o in res[0]["features"]:
                feature = genome.features.get_by_id(o["id"])
                if "function" in o:
                    rast_function = o["function"]
                    if split_terms:
                        for function in split_annotation(rast_function):
                            feature.add_ontology_term("RAST", function)
                    else:
                        feature.add_ontology_term("RAST", rast_function)

        return res[0]["analysis_events"] if res else None

    def annotate_genome_from_fasta(self, filepath, split="|"):
        genome = MSGenome.from_fasta(filepath, split)
        res = self.annotate_genome(genome)

        return genome, res

    def annotate_protein_sequence(self, protein_id: str, protein_seq: str):
        p_features = [{"id": protein_id, "protein_translation": protein_seq}]
        return self.f(p_features)

    def annotate_protein_batch(self, protein_batch: dict):
        protein_annotations = {}
        features = [
            {"id": _id, "protein_translation": seq}
            for _id, seq in protein_batch.items()
        ]
        params = [{"features": features}, {"stages": self.stages}]
        result = self.rpc_client.call("GenomeAnnotation.run_pipeline", params)
        if result and len(result) > 0 and "features" in result[0]:
            for feature in result[0]["features"]:
                if "function" in feature:
                    protein_annotations[feature["id"]] = feature["function"]

        return protein_annotations

    def annotate_protein_sequences(self, protein_seqs: dict, chunk_size: int = 5000):
        protein_annotations = {}
        chunk = {}

        for protein_id, protein_seq in protein_seqs.items():
            chunk[protein_id] = protein_seq

            # Process chunk when it reaches the specified size
            if len(chunk) >= chunk_size:
                protein_annotations.update(self.annotate_protein_batch(chunk))
                chunk = {}

        # Process remaining proteins in the last chunk
        if chunk:
            protein_annotations.update(self.annotate_protein_batch(chunk))

        return protein_annotations

    def f1(self, protein_id, protein_seq):
        p_features = [{"id": protein_id, "protein_translation": protein_seq}]
        return self.f(p_features)

    def f(self, p_features):
        params = [{"features": p_features}, {"stages": self.stages}]
        result = self.rpc_client.call("GenomeAnnotation.run_pipeline", params)
        return result
