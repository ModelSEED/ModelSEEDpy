from modelseedpy.core.rpcclient import RPCClient


class RastClient:
    
    def __init__(self):
        self.rpc_client = RPCClient("https://tutorial.theseed.org/services/genome_annotation")
        self.stages = [
            {"name": 'annotate_proteins_kmer_v2', "kmer_v2_parameters": {}},
            {"name": 'annotate_proteins_kmer_v1', "kmer_v1_parameters": {"annotate_hypothetical_only": 1}},
            {"name": 'annotate_proteins_similarity', "similarity_parameters": {"annotate_hypothetical_only": 1}}
        ]

    def f1(self, protein_id, protein_seq):
        p_features = [
            {
                "id": protein_id,
                "protein_translation": protein_seq
            }
        ]
        return self.f(p_features)

    def f(self, p_features):
        params = [
            {"features": p_features},
            {"stages": self.stages}
        ]
        result = self.rpc_client.call("GenomeAnnotation.run_pipeline", params)
        return result

    def fasta(self, filepath):
        # process FASTA
        p_features = []
        return self.f(p_features)
