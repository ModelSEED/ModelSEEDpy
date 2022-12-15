# -*- coding: utf-8 -*-
import os
import pytest
import json
import cobra
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy import MSATPCorrection, MSMedia


@pytest.fixture
def template():
    with open(
        os.path.join(
            os.path.dirname(__file__), "..", "test_data", "template_core_bigg.json"
        ),
        "r",
    ) as fh:
        return MSTemplateBuilder.from_dict(json.load(fh)).build()


@pytest.fixture
def template_genome_scale():
    with open(
        os.path.join(
            os.path.dirname(__file__),
            "..",
            "test_data",
            "template_genome_scale_bigg.json",
        ),
        "r",
    ) as fh:
        return MSTemplateBuilder.from_dict(json.load(fh)).build()


@pytest.fixture
def get_model():
    def _method(ko=None, added_compounds=None, added_reactions=None):
        if ko is None:
            ko = []
        with open(
            os.path.join(
                os.path.dirname(__file__), "..", "test_data", "e_coli_core.json"
            ),
            "r",
        ) as fh:
            model_json = json.load(fh)
            model_json["compartments"] = {
                k + "0": v for (k, v) in model_json["compartments"].items()
            }
            metabolites = {}
            for m in model_json["metabolites"]:
                m["id"] += "0"
                m["compartment"] += "0"
                metabolites[m["id"]] = m
            for r in model_json["reactions"]:
                r["metabolites"] = {i + "0": v for (i, v) in r["metabolites"].items()}
                compartments = set(
                    [metabolites[k]["compartment"] for k in r["metabolites"].keys()]
                )
                if r["id"].endswith("_e"):
                    r["id"] += "0"
                elif len(compartments) == 1:
                    r["id"] += "_" + list(compartments)[0]
                else:
                    r["id"] += (
                        "_" + "c0"
                    )  # hack cause there is only combo between e0 and c0

            model_json["reactions"] = [
                x for x in model_json["reactions"] if x["id"] not in ko
            ]

            if added_compounds:
                for o in added_compounds:
                    model_json["metabolites"].append(o)
            if added_reactions:
                for o in added_reactions:
                    model_json["reactions"].append(o)
            model = cobra.io.from_json(json.dumps(model_json))
            model.reactions.ATPM_c0.lower_bound = 0
            model.reactions.ATPM_c0.upper_bound = 1000
            return model

    return _method


@pytest.fixture
def media_glucose_aerobic():
    media = MSMedia.from_dict(
        {
            "glc__D": (-1, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
        }
    )
    media.id = "glc/o2"
    return media


@pytest.fixture
def media_acetate_aerobic():
    media = MSMedia.from_dict(
        {
            "ac": (-1, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
        }
    )
    media.id = "glc/o2"
    return media


@pytest.fixture
def media_genome_scale_glucose_aerobic():
    media = MSMedia.from_dict(
        {
            "glc__D": (-10, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
            "pi": (-1000, 1000),
            "co2": (-1000, 1000),
            "nh4": (-1000, 1000),
            "k": (-1000, 1000),
        }
    )
    return media


@pytest.fixture
def media_all_aerobic(media_glucose_aerobic, media_acetate_aerobic):
    return [media_glucose_aerobic, media_acetate_aerobic]


@pytest.fixture
def get_model_with_infinite_atp_loop(get_model, template_genome_scale):
    def _method(ko=None):
        added_compounds = [
            {
                "id": "k_c0",
                "name": "K [c]",
                "compartment": "c0",
                "charge": 1,
                "formula": "K",
                "notes": {},
                "annotation": {"sbo": "SBO:0000247"},
            },
            {
                "id": "k_e0",
                "name": "K [e]",
                "compartment": "e0",
                "charge": 1,
                "formula": "K",
                "notes": {},
                "annotation": {"sbo": "SBO:0000247"},
            },
        ]
        added_reactions = [
            {
                "id": "EX_k_e0",
                "name": "K exchange",
                "metabolites": {"k_e0": -1.0},
                "lower_bound": -1000,
                "upper_bound": 1000.0,
                "gene_reaction_rule": "",
                "subsystem": "Extracellular exchange",
                "notes": {},
                "annotation": {},
            }
        ]
        model = get_model(ko, added_compounds, added_reactions)
        model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM_c0").add_metabolites(
            {model.metabolites.get_by_id("k_c0"): 0.01}
        )
        model.add_reactions([template_genome_scale.reactions.Kt2r_c.to_reaction(model)])
        model.add_reactions([template_genome_scale.reactions.Ktex_c.to_reaction(model)])
        model.medium = {
            "EX_co2_e0": 1000.0,
            "EX_glc__D_e0": 10.0,
            "EX_h_e0": 1000.0,
            "EX_h2o_e0": 1000.0,
            "EX_nh4_e0": 1000.0,
            "EX_o2_e0": 1000.0,
            "EX_pi_e0": 1000.0,
            "EX_k_e0": 1000.0,
        }

        return model

    return _method


def test_infinite_atp_model_growth_boost(
    get_model_with_infinite_atp_loop,
    template,
    template_genome_scale,
    media_glucose_aerobic,
):
    model = get_model_with_infinite_atp_loop()
    solution = model.optimize()
    assert solution.objective_value > 1.2  # should be around 1.316
    assert solution.status == "optimal"


def test_ms_atp_correction1(get_model, template, media_all_aerobic):
    model = get_model(["GLCpts_c0", "NADH16_c0", "CYTBD_c0", "O2t_c0"])
    atp_correction = MSATPCorrection(
        model,
        template,
        media_all_aerobic,
        atp_hydrolysis_id="ATPM_c0",
        load_default_medias=False,
    )
    atp_correction.evaluate_growth_media()
    assert len(atp_correction.noncore_reactions) == 1  # the biomass
    assert len(atp_correction.other_compartments) == 0  # none
    assert len(atp_correction.original_bounds) == 70  # 70 reactions

    """
    glc/o2 {'reversed': {}, 'new': {'GLCpts_c0': '>'}}
    ac/o2 {'reversed': {}, 'new': {'CYTBD_c0': '>', 'NADH16_c0': '>', 'O2t_c0': '>'}}
    """
    atp_correction.determine_growth_media()

    assert len(atp_correction.selected_media) == 1  # selects glucose

    atp_correction.apply_growth_media_gapfilling()

    media_eval = atp_correction.evaluate_growth_media()

    atp_correction.expand_model_to_genome_scale()
    tests = atp_correction.build_tests()

    assert tests
    assert len(tests) == 1
    assert tests[0]["threshold"] > 0
    assert tests[0]["objective"] == "ATPM_c0"


def test_ms_atp_correction_and_gap_fill1(
    get_model_with_infinite_atp_loop,
    template,
    template_genome_scale,
    media_glucose_aerobic,
    media_genome_scale_glucose_aerobic,
):
    from modelseedpy import MSGapfill

    model = get_model_with_infinite_atp_loop(["GLCpts_c0", "GLUSy_c0", "GLUDy_c0"])
    model.reactions.ATPM_c0.lower_bound = 0
    model.reactions.ATPM_c0.upper_bound = 1000

    atp_correction = MSATPCorrection(
        model,
        template,
        [media_glucose_aerobic],
        atp_hydrolysis_id="ATPM_c0",
        load_default_medias=False,
    )
    tests = atp_correction.run_atp_correction()

    # expected tests = [{'media': MSMedia object, 'is_max_threshold': True, 'threshold': 21.0, 'objective': 'ATPM_c0'}]

    assert tests
    assert len(tests) == 1
    assert tests[0]["threshold"] > 0
    assert tests[0]["objective"] == "ATPM_c0"

    gap_fill = MSGapfill(model, [template_genome_scale], [], tests, {}, [])
    result = gap_fill.run_gapfilling(
        media_genome_scale_glucose_aerobic,
        "BIOMASS_Ecoli_core_w_GAM_c0",
        minimum_obj=0.1,
    )

    # either GLUSy_c0 or GLUDy_c0 should be gap filled for glutamate

    assert result
    assert len(result["new"]) == 1
    assert "GLUSy_c0" in result["new"] or "GLUDy_c0" in result["new"]

    gap_fill.integrate_gapfill_solution(result)

    # TODO: add some model testing assertion
