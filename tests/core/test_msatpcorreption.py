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
    media.id = "ac/o2"
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
    atp_hydrolysis_id = "ATPM_c0"
    model = get_model(["GLCpts_c0", "NADH16_c0", "CYTBD_c0", "O2t_c0"])
    atp_correction = MSATPCorrection(
        model,
        template,
        media_all_aerobic,
        atp_hydrolysis_id=atp_hydrolysis_id,
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
    assert len(tests) == 2  # glucose and empty
    for t in tests:
        if t["media"].id == "empty":
            assert t["threshold"] <= 1e-05
        else:
            assert t["threshold"] > 1e-05
        assert t["objective"] == atp_hydrolysis_id
        assert t["is_max_threshold"] is True


def test_ms_atp_correction_and_gap_fill1(
    get_model_with_infinite_atp_loop,
    template,
    template_genome_scale,
    media_glucose_aerobic,
    media_genome_scale_glucose_aerobic,
):
    from modelseedpy import MSGapfill

    atp_hydrolysis_id = "ATPM_c0"

    model = get_model_with_infinite_atp_loop(["GLCpts_c0", "GLUSy_c0", "GLUDy_c0"])
    model.reactions.ATPM_c0.lower_bound = 0
    model.reactions.ATPM_c0.upper_bound = 1000
    model.objective = atp_hydrolysis_id
    atp_correction = MSATPCorrection(
        model,
        template,
        [media_glucose_aerobic],
        atp_hydrolysis_id=atp_hydrolysis_id,
        load_default_medias=False,
    )
    tests = atp_correction.run_atp_correction()
    # expected tests = [{'media': MSMedia object, 'is_max_threshold': True, 'threshold': 21.0, 'objective': 'ATPM_c0'}]

    assert tests
    assert len(tests) == 2
    for t in tests:
        if t["media"].id == "empty":
            assert t["threshold"] <= 1e-05
        else:
            assert t["threshold"] > 1e-05
        assert t["objective"] == atp_hydrolysis_id
        assert t["is_max_threshold"] is True

    model.objective = "BIOMASS_Ecoli_core_w_GAM_c0"
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


def test_build_tests_preserves_expansion_filter_when_no_media_selected(
    get_model_with_infinite_atp_loop, template, media_glucose_aerobic
):
    """``build_tests`` must NOT restore non-core bounds once the ATP expansion
    filter has been applied.

    Regression for KBDLJobRunningPrototype#28, cause 2.
    ``run_atp_correction`` ends ``return self.build_tests()`` -- that is, AFTER
    ``expand_model_to_genome_scale`` has blocked each cycle-forming reaction by
    zeroing one direction.  ``build_tests`` used to open with a truthiness test
    on ``selected_media``, so a model where NO medium qualified re-entered the
    "not yet computed" branch, and that branch ends in
    ``restore_noncore_reactions()`` -- which writes ``original_bounds`` back
    over every non-core reaction and so silently undoes the filter, returning
    energy-generating cycles to the delivered model.  Measured on KBDL builds:
    ``rxn08941_c0`` went from blocked back to [-1000, 1000], and free ATP from
    0 to 1000.

    An empty ``selected_media`` is a real answer ("evaluated, nothing
    qualified"), not a missing one, so the guard must distinguish it from
    ``None`` ("never evaluated").  This pins the ``[]`` side; the ``None`` side
    is pinned by :func:`test_build_tests_computes_when_never_evaluated`.

    The filter is applied here the same way ``expand_model_to_genome_scale``
    applies it -- zero one direction of a reaction that ``disable_noncore_
    reactions`` has recorded in ``original_bounds`` -- rather than by running
    the genome-scale expansion, so the test needs no solver-dependent
    expansion result and no local biochemistry database.
    """
    model = get_model_with_infinite_atp_loop()
    atp_correction = MSATPCorrection(
        model,
        template,
        [media_glucose_aerobic],
        atp_hydrolysis_id="ATPM_c0",
        load_default_medias=False,
    )
    # Populates noncore_reactions and original_bounds, then selects media.
    atp_correction.evaluate_growth_media(no_gapfilling=True)
    atp_correction.determine_growth_media()
    atp_correction.restore_noncore_reactions(noncore=True, othercompartment=False)

    # Pick a victim shaped like the reactions that actually formed the cycles:
    # REVERSIBLE (original bounds negative-to-positive, as rxn08941_c0 was at
    # [-1000, 1000]), non-core, and still in the model, so that a restore is
    # unambiguously observable.  Biomass and exchange reactions are excluded --
    # ``disable_noncore_reactions`` means to skip them, and only misses the
    # biomass here because ``FBAHelper.is_biomass`` matches a lowercase "bio"
    # prefix while this fixture's id is uppercase.
    victim = None
    for reaction, _direction in atp_correction.noncore_reactions:
        bounds = atp_correction.original_bounds.get(reaction.id)
        if not bounds or reaction.id not in model.reactions:
            continue
        if reaction.id.lower().startswith(("bio", "ex_", "dm_", "sk_")):
            continue
        if bounds[0] < 0 and bounds[1] > 0:
            victim = reaction
            break
    if victim is None:
        pytest.skip("fixture produced no reversible non-core reaction to filter")
    original_upper = atp_correction.original_bounds[victim.id][1]
    assert original_upper > 0

    # Apply the filter exactly as expand_model_to_genome_scale does.
    atp_correction.filtered_noncore = [(victim, ">")]
    victim.upper_bound = 0

    # The exact state a KBDL build reaches: media evaluated, none qualified.
    atp_correction.selected_media = []
    atp_correction.build_tests()

    assert model.reactions.get_by_id(victim.id).upper_bound == 0, (
        f"build_tests restored {victim.id} to upper_bound "
        f"{model.reactions.get_by_id(victim.id).upper_bound} (original "
        f"{original_upper}), erasing the ATP expansion filter"
    )


def test_build_tests_computes_when_never_evaluated(
    get_model, template, media_glucose_aerobic
):
    """The other side of the same gate: a caller asking for tests before any
    evaluation must still get the media evaluated.

    ``selected_media`` is ``None`` until ``determine_growth_media`` runs, and
    that is the case the "not yet computed" branch exists for.  Narrowing the
    guard to ``is None`` must not disable it.
    """
    model = get_model()
    atp_correction = MSATPCorrection(
        model,
        template,
        [media_glucose_aerobic],
        atp_hydrolysis_id="ATPM_c0",
        load_default_medias=False,
    )
    assert atp_correction.selected_media is None
    atp_correction.build_tests()
    assert atp_correction.selected_media is not None
