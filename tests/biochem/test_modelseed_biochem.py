# -*- coding: utf-8 -*-

"""
def test_compound_data_h2o(cpd):
    assert cpd.id == 'cpd00001'
    assert cpd.abbr == 'h2o'
    assert cpd.name == 'H2O'
    assert cpd.formula == 'H2O'
    assert cpd.source == 'Primary Database'
    assert cpd.inchi_key == 'XLYOFNOQVPJJNP-UHFFFAOYSA-N'
    assert cpd.mass == 18.0
    assert cpd.charge == 0
    assert cpd.is_core == True
    assert cpd.is_obsolete == False
    assert cpd.linked_compound == None
    assert cpd.is_cofactor == False
    assert cpd.delta_g == -37.54
    assert cpd.pka == '1:1:15.70'
    assert cpd.pkb == '1:1:-1.80'
    assert cpd.is_abstract == False
    assert cpd.inchi == 'InChI=1S/H2O/h1H2'
    assert cpd.annotation == {
        'bigg.metabolite': {'h2o', 'oh1'},
        'kegg.compound': {'C00001', 'C01328'},
        'metacyc.compound': {'OH', 'OXONIUM', 'WATER'},
        'metanetx.chemical': {'MNXM2'}
    }
    # assert cpd.comprised_of == False
    # assert cpd.aliases == False

    assert cpd.smiles == 'O'
    assert cpd.flags == {'GC', 'EQ', 'EQU'}

    #names !
    # {'H20',
    #  'H2O',
    #  'H3O+',
    #  'HO-',
    #  'Hydroxide ion',
    #  'OH',
    #  'OH-',
    #  'Water',
    #  'hydrogen oxide',
    #  'hydroxide',
    #  'hydroxide ion',
    #  'hydroxyl',
    #  'hydroxyl ion',
    #  'oxonium',
    #  'water'}

# The following was another doc string.

rxn = reactions[0]
assert rxn.id == 'rxn00001'
assert rxn.abbr == 'R00004'
assert rxn.name == 'diphosphate phosphohydrolase'
assert rxn.delta_g == -3.46
assert rxn.delta_g_error == 0.05
assert rxn.status == 'OK'
assert rxn.is_obsolete == False
assert rxn.is_abstract == False
assert rxn.source == 'Primary Database'
assert rxn.flags == {'EQC', 'EQU', 'GCC', 'HB'}
assert rxn.annotation == {
    'ec-code': {'3.6.1.1'},
    'bigg.reaction': {'IPP1', 'PPA', 'PPA_1', 'PPAm'},
    'kegg.reaction': {'R00004'},
    'metacyc.reaction': {'INORGPYROPHOSPHAT-RXN'},
    'metanetx.reaction': {'MNXR100808'},
    'rhea': {'24579'}
}
assert rxn.names == {
    'Diphosphate phosphohydrolase',
    'INORGPYROPHOSPHAT-RXN',
    'Inorganic diphosphatase',
    'Inorganic pyrophosphatase',
    'PPA',
    'Pyrophosphate phosphohydrolase',
    'diphosphate phosphohydrolase',
    'inorganic diphosphatase',
    'inorganic diphosphatase (one proton translocation)',
    'inorganicdiphosphatase',
    'pyrophosphate phosphohydrolase'
}

TODO: Fix to actually test something.

rxn = reactions[0]
assert rxn.code == '(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0]'
assert rxn.stoichiometry == '-1:cpd00001:0:0:"H2O";-1:cpd00012:0:0:"PPi";2:cpd00009:0:0:"Phosphate";1:cpd00067:0:0:"H+"'
assert rxn.is_transport == False
assert rxn.equation == '(1) cpd00001[0] + (1) cpd00012[0] <=> (2) cpd00009[0] + (1) cpd00067[0]'
assert rxn.definition == '(1) H2O[0] + (1) PPi[0] <=> (2) Phosphate[0] + (1) H+[0]'
assert rxn.reversibility == '>'
assert rxn.direction == '='
assert rxn.pathways == 'MetaCyc: Degradation (Degradation/Utilization/Assimilation); Glyphosate-Degradation (glyphosate degradation); Noncarbon-Nutrients (Inorganic Nutrient Metabolism); PWY-7805 ((aminomethyl)phosphonate degradation); PWY-7807 (glyphosate degradation III); Phosphorus-Compounds (Phosphorus Compound Metabolism)'
assert rxn.compound_ids == {'cpd00001', 'cpd00009', 'cpd00012', 'cpd00067'}
assert rxn.linked_reaction == 'rxn27946;rxn27947;rxn27948;rxn32487;rxn38157;rxn38158'
"""
