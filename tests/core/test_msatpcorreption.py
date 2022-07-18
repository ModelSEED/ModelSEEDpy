import pytest
import json
import cobra
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy import MSATPCorrection, MSMedia


@pytest.fixture
def template():
    with open('./tests/test_data/template_core_bigg.json', 'r') as fh:
        return MSTemplateBuilder.from_dict(json.load(fh)).build()


@pytest.fixture
def get_model():

    def _method(ko=None):
        if ko is None:
            ko = []
        with open('./tests/test_data/e_coli_core.json', 'r') as fh:
            model_json = json.load(fh)
            model_json['compartments'] = {k + '0': v for (k, v) in model_json['compartments'].items()}
            metabolites = {}
            for m in model_json['metabolites']:
                m['id'] += '0'
                m['compartment'] += '0'
                metabolites[m['id']] = m
            for r in model_json['reactions']:
                r['metabolites'] = {i + '0': v for (i, v) in r['metabolites'].items()}
                compartments = set([metabolites[k]['compartment'] for k in r['metabolites'].keys()])
                if r['id'].endswith('_e'):
                    r['id'] += '0'
                elif len(compartments) == 1:
                    r['id'] += '_' + list(compartments)[0]
                else:
                    r['id'] += '_' + 'c0'  # hack cause there is only combo between e0 and c0

            model_json['reactions'] = [x for x in model_json['reactions'] if x['id'] not in ko]
            model = cobra.io.from_json(json.dumps(model_json))
            model.reactions.ATPM_c0.lower_bound = 0
            model.reactions.ATPM_c0.upper_bound = 1000
            return model

    return _method


@pytest.fixture
def media_glucose_aerobic():
    media = MSMedia.from_dict({
        'glc__D': (-1, 1000),
        'o2': (-1000, 1000),
        'h': (-1000, 1000),
        'h2o': (-1000, 1000)
    })
    media.id = 'glc/o2'
    return media


@pytest.fixture
def media_acetate_aerobic():
    media = MSMedia.from_dict({
        'ac': (-1, 1000),
        'o2': (-1000, 1000),
        'h': (-1000, 1000),
        'h2o': (-1000, 1000)
    })
    media.id = 'glc/o2'
    return media


@pytest.fixture
def media_all_aerobic(media_glucose_aerobic, media_acetate_aerobic):
    return [media_glucose_aerobic, media_acetate_aerobic]


def test_ms_atp_correction1(get_model, template, media_all_aerobic):
    model = get_model(["GLCpts_c0", "NADH16_c0", "GLCpts_c0", "O2t_c0"])
    atp_correction = MSATPCorrection(model, template, media_all_aerobic, atp_hydrolysis_id='ATPM_c0')
    atp_correction.evaluate_growth_media()
    print('non_core_reactions', len(atp_correction.noncore_reactions),
          'other_compartments', len(atp_correction.other_compartments),
          'original_bounds', len(atp_correction.original_bounds))
    for media in atp_correction.media_gapfill_stats:
        print(media, atp_correction.media_gapfill_stats[media])
    atp_correction.determine_growth_media()
    print(atp_correction.selected_media)
    media_eval = atp_correction.evaluate_growth_media()
    print(media_eval)
    atp_correction.expand_model_to_genome_scale()
    tests = atp_correction.build_tests()
    print(tests)
    assert True
