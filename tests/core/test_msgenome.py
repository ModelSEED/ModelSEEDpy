from modelseedpy.core.msgenome import MSFeature, MSGenome, read_fasta


def test_read_fasta():
    features = read_fasta('../../examples/GCF_000005845.2.faa', ' ')
    assert len(features) == 3


def test_msfeature_add_ontology_term():
    feature = MSFeature('gene1', 'MKV')
    feature.add_ontology_term('ONTOLOGY', 'value1')
    assert 'ONTOLOGY' in feature.ontology_terms
    assert 'value1' in feature.ontology_terms['ONTOLOGY']


def test_msgenome_from_protein_sequences_hash1():
    genome = MSGenome.from_protein_sequences_hash({'gene1': 'MKV'})
    assert len(genome.features) == 1
    assert genome.features[0]
    assert genome.features[0].id == 'gene1'
    assert genome.features[0].seq == 'MKV'


def test_msgenome_from_protein_sequences_hash2():
    genome = MSGenome.from_protein_sequences_hash({'gene1': 'MKV', 'gene2': 'MKVLGD'})
    assert len(genome.features) == 2
