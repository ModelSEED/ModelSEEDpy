from modelseedpy.core.rpcclient import RPCClient
from modelseedpy.core.msgenome import Sequence
import re
from cobrakbase.core.kbasegenomesgenome import normalize_role #move this to this lib

def read_fasta(f, split='|'):
    features = []
    with open(f, 'r') as fh:
        lines = fh.read().split('\n')
        seq = None
        for line in lines:
            if line.startswith('>'):
                if seq:
                    features.append(seq)
                seq_id = line[1:]
                desc = ""
                if split:
                    header_data = line[1:].split(split)
                    seq_id = header_data[0]
                    desc = header_data[1]

                seq = Sequence(seq_id, "", desc)
            else:
                if seq:
                    seq.seq += line
    return features


def aux_rast_result(res, g):
    search_name_to_genes = {}
    search_name_to_orginal = {}
    for f in res[0]['features']:
        function = f['function'] if 'function' in f else None
        if not function:
            function = g.features.get_by_id(f['id']).description

        if function and len(function) > 0:
            functions = re.split('; | / | @ | => ', function)
            for f_ in functions:
                f_norm = normalize_role(f_)
                if f_norm not in search_name_to_genes:
                    search_name_to_genes[f_norm] = set()
                    search_name_to_orginal[f_norm] = set()
                search_name_to_orginal[f_norm].add(f_)
                search_name_to_genes[f_norm].add(f['id'])

    return search_name_to_genes, search_name_to_orginal


def aux_template(template):
    rxn_roles = {}
    roles = dict(map(lambda x: (x['id'], x), template.roles))
    for r in template.reactions:
        rxn_roles[r.id] = set()
        complex_roles = r.get_complex_roles()
        if len(complex_roles) > 0:
            for cpx_id in complex_roles:
                for role_id in complex_roles[cpx_id]:
                    rxn_roles[r.id].add(normalize_role(roles[role_id]['name']))
                    # print(role_id, normalize_role(roles[role_id]['name']))
    return rxn_roles


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

    def fasta(self, filepath, split='|'):
        seqs = read_fasta(filepath, split)
        print(len(seqs))
        seqs = list(filter(lambda x: len(x.seq) > 0, seqs))
        print(len(seqs))
        p_features = []
        for s in seqs:
            p_features.append({
                'id': s.id,
                'protein_translation': s.seq,
            })
        return self.f(p_features)
