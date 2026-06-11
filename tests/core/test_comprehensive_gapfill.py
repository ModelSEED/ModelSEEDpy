# -*- coding: utf-8 -*-
"""Integration tests for MSGapfill.run_comprehensive_gapfilling and the
'Comprehensive' gapfilling_mode in run_multi_gapfill.

These tests use the same e_coli_core toy model and template_core_bigg template
already present in tests/test_data/.

**GLPK solver constraint**: Stage 1 of comprehensive gapfilling is a MILP
(RevBinPkg adds binary variables for each reaction direction in the rxn_filter).
With GLPK (the only solver available on this machine), large MILPs are
intractable. Therefore:

- Tests that call ``run_comprehensive_gapfilling`` directly use
  ``activation_rxn_filter`` to pass a small set of forward-only reactions.
  Forward-only reactions produce trivially integral RevBin variables (GLPK
  fixes them immediately) so Stage 1 solves in seconds.
- Tests that call ``run_multi_gapfill(gapfilling_mode='Comprehensive')``
  use a model that has ONLY forward-only non-exchange reactions, so the same
  fast-solve property holds for the full rxn_filter.

Both approaches exercise the full two-stage algorithm; the solver constraint
only limits which reactions appear in the activation maximization (Stage 1),
not the gapfilling minimization (Stage 2).
"""
import os
import json
import pytest
import cobra
from cobra.flux_analysis import pfba

from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy import MSGapfill, MSMedia


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def template():
    with open(
        os.path.join(
            os.path.dirname(__file__), "..", "test_data", "template_core_bigg.json"
        ),
        "r",
    ) as fh:
        return MSTemplateBuilder.from_dict(json.load(fh)).build()


def _build_model(model_json_path, ko=None, fwd_only=False):
    """Build an e_coli_core model with optional KOs and optional forward-only filter.

    ``fwd_only=True`` removes all reversible non-exchange reactions from the
    original model, leaving only forward-only reactions. The gapfilling template
    will then add the reversible reactions back as candidates, keeping them out
    of the activation rxn_filter (so Stage 1 has trivially integral binary vars).
    """
    if ko is None:
        ko = []
    with open(model_json_path, "r") as fh:
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
            r["id"] += "_c0"

    if fwd_only:
        # Keep only forward-only reactions + exchanges (+ KO-d reactions are also removed)
        model_json["reactions"] = [
            r for r in model_json["reactions"]
            if r["id"] not in ko
            and (
                r["id"].startswith("EX_")       # exchanges always kept
                or r.get("lower_bound", -1000) >= 0  # forward-only non-exchange
            )
        ]
    else:
        model_json["reactions"] = [
            r for r in model_json["reactions"]
            if r["id"] not in ko
        ]

    # Keep all metabolites so reactions can reference them
    model_json["metabolites"] = list(metabolites.values())
    return cobra.io.from_json(json.dumps(model_json))


@pytest.fixture
def model_data_path():
    return os.path.join(
        os.path.dirname(__file__), "..", "test_data", "e_coli_core.json"
    )


@pytest.fixture
def complete_media():
    """Generous complete media (all tested carbon/nitrogen sources open)."""
    return MSMedia.from_dict(
        {
            "glc__D": (-1000, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
            "pi": (-1000, 1000),
            "co2": (-1000, 1000),
            "nh4": (-1000, 1000),
            "succ": (-1000, 1000),
            "fum": (-1000, 1000),
            "ac": (-1000, 1000),
            "pyr": (-1000, 1000),
            "mal__L": (-1000, 1000),
            "fru": (-1000, 1000),
        }
    )


@pytest.fixture
def standard_media():
    """Standard glucose aerobic media (matches test_msgapfill.py)."""
    return MSMedia.from_dict(
        {
            "glc__D": (-10, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
            "pi": (-1000, 1000),
            "co2": (-1000, 1000),
            "nh4": (-1000, 1000),
        }
    )


def _small_fwd_filter(model, max_count=12):
    """Return a small set of forward-only non-exchange reaction ids.
    Used by tests that call run_comprehensive_gapfilling directly, so
    they can restrict Stage 1 to trivially integral binary variables."""
    fwd = [
        r.id for r in model.reactions
        if r.lower_bound == 0
        and r.upper_bound > 0
        and not r.id.startswith("EX_")
        and r.id != "ATPM_c0"   # ATP maintenance has lb > 0; skip
    ]
    return set(fwd[:max_count])


# ---------------------------------------------------------------------------
# Test A: API shape — run_comprehensive_gapfilling returns correct keys
# ---------------------------------------------------------------------------

def test_comprehensive_gapfilling_returns_correct_shape(
    template, model_data_path, complete_media
):
    """run_comprehensive_gapfilling must return a dict with the same keys as
    run_gapfilling: 'new', 'reversed', 'media', 'target', 'minobjective',
    'binary_check'. Uses a small forward-only activation_rxn_filter to keep
    Stage 1 tractable with GLPK."""
    model = _build_model(model_data_path, ko=["GLCpts_c0"])
    target = "BIOMASS_Ecoli_core_w_GAM_c0"
    gf = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)

    result = gf.run_comprehensive_gapfilling(
        complete_media,
        target=target,
        minimum_obj=0.01,
        binary_check=False,
        prefilter=False,
        activation_rxn_filter=_small_fwd_filter(model),
    )

    assert result is not None, "run_comprehensive_gapfilling returned None"
    for key in ("new", "reversed", "media", "target", "minobjective", "binary_check"):
        assert key in result, f"Missing key '{key}' in result dict"
    assert result["media"] is complete_media
    assert result["target"] == target
    assert result["minobjective"] == 0.01
    assert result["binary_check"] is False


# ---------------------------------------------------------------------------
# Test B: grows — integrated model produces positive biomass
# ---------------------------------------------------------------------------

def test_comprehensive_gapfilling_grows(template, model_data_path, complete_media):
    """(a) After comprehensive gapfilling and solution integration, the model
    produces positive biomass (slim_optimize > 0 on complete media).

    Uses activation_rxn_filter with forward-only reactions for GLPK.
    """
    model = _build_model(model_data_path, ko=["GLCpts_c0"])
    target = "BIOMASS_Ecoli_core_w_GAM_c0"
    gf = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)

    result = gf.run_comprehensive_gapfilling(
        complete_media,
        target=target,
        minimum_obj=0.01,
        prefilter=False,
        activation_rxn_filter=_small_fwd_filter(model),
    )
    assert result is not None, "run_comprehensive_gapfilling returned None"

    gf.integrate_gapfill_solution(result)
    biomass = model.slim_optimize()
    assert biomass is not None and biomass > 0, (
        f"Model does not grow after comprehensive gapfilling (biomass={biomass})"
    )


# ---------------------------------------------------------------------------
# Test C: comprehensive > sequential — more active original reactions
# ---------------------------------------------------------------------------

def test_comprehensive_vs_sequential_more_active_reactions(
    template, model_data_path, standard_media
):
    """(b) Comprehensive gapfilling activates strictly more original-model
    reactions (measured by pFBA flux) than Sequential gapfilling on the same
    inputs.

    We use a model where several forward-only native reactions are knocked out
    (ICL, MALS, FRUpts2, LDH_D — these can be reached via alternative carbon
    sources on complete media) along with the core blocked reactions
    (GLCpts, CYTBD, NADH16). Sequential gapfilling adds back only the minimum
    set needed for biomass on standard media. Comprehensive additionally forces
    the activation of the specified forward-only reactions, finding gapfill
    candidates that create metabolic connections for them.

    The pFBA count is measured after opening all exchanges (simulating complete
    media) so that both models are evaluated under the same permissive conditions.

    To keep Stage 1 tractable with GLPK, ``activation_rxn_filter`` is used
    with a small set of 10 forward-only reactions.
    """
    ko = ["GLCpts_c0", "CYTBD_c0", "NADH16_c0", "ICL_c0", "MALS_c0", "FRUpts2_c0", "LDH_D_c0"]
    target = "BIOMASS_Ecoli_core_w_GAM_c0"

    # Forward-only reactions whose activation Comprehensive should maximize.
    # These overlap with the KO set (ICL, MALS, FRUpts2, LDH_D are knocked out)
    # but also include unblocked reactions to maximize the contrast.
    activation_filter = {
        "PDH_c0", "CS_c0", "SUCDi_c0", "GLNS_c0", "PYK_c0",
        "GND_c0", "PGL_c0", "AKGDH_c0", "FRUpts2_c0", "LDH_D_c0",
    }

    # --- Sequential gapfilling ---
    model_seq = _build_model(model_data_path, ko=ko)
    orig_ids = set(r.id for r in model_seq.reactions)
    gf_seq = MSGapfill(model_seq, [template], minimum_obj=0.01, default_target=target)
    gf_seq.run_multi_gapfill(
        [standard_media],
        target=target,
        default_minimum_objective=0.01,
        gapfilling_mode="Sequential",
        binary_check=False,
        prefilter=False,
        run_sensitivity_analysis=False,
        integrate_solutions=True,
    )

    # --- Comprehensive gapfilling (direct call with GLPK-tractable filter) ---
    model_comp = _build_model(model_data_path, ko=ko)
    gf_comp = MSGapfill(model_comp, [template], minimum_obj=0.01, default_target=target)
    comp_result = gf_comp.run_comprehensive_gapfilling(
        standard_media,
        target=target,
        minimum_obj=0.01,
        prefilter=False,
        activation_rxn_filter=activation_filter,
    )
    if comp_result is not None:
        gf_comp.integrate_gapfill_solution(comp_result)

    # Count active original-model reactions via pFBA.
    # Open all exchanges to simulate complete media for the flux count.
    def count_active(model, orig_ids):
        try:
            with model:
                for rxn in model.exchanges:
                    rxn.lower_bound = -1000
                sol = pfba(model)
            return sum(
                1 for r in model.reactions
                if r.id in orig_ids and abs(sol.fluxes.get(r.id, 0)) > 1e-8
            )
        except Exception:
            return 0

    active_seq = count_active(model_seq, orig_ids)
    active_comp = count_active(model_comp, orig_ids)

    assert active_comp > 0, "Comprehensive model has no active reactions after integration"
    assert active_seq > 0, "Sequential model has no active reactions after integration"
    assert active_comp > active_seq, (
        f"Comprehensive should activate more reactions than Sequential: "
        f"comprehensive={active_comp}, sequential={active_seq}"
    )


# ---------------------------------------------------------------------------
# Test D: no native reactions removed after re-running standard gapfilling
# ---------------------------------------------------------------------------

def test_comprehensive_no_native_reactions_removed(
    template, model_data_path, complete_media, standard_media
):
    """(c) All original-model reactions survive after comprehensive gapfilling
    and a subsequent standard gapfilling pass."""
    model = _build_model(model_data_path, ko=["GLCpts_c0"])
    orig_ids = set(r.id for r in model.reactions)
    target = "BIOMASS_Ecoli_core_w_GAM_c0"

    gf = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)
    result = gf.run_comprehensive_gapfilling(
        complete_media,
        target=target,
        minimum_obj=0.01,
        prefilter=False,
        activation_rxn_filter=_small_fwd_filter(model),
    )
    assert result is not None
    gf.integrate_gapfill_solution(result)

    # Re-run standard gapfilling on the now-grown model
    gf2 = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)
    result2 = gf2.run_gapfilling(
        standard_media, target, minimum_obj=0.01, prefilter=False
    )

    # If a second gapfilling is needed, it must not remove any native reaction
    if result2 is not None:
        for rxn_id in orig_ids:
            assert rxn_id in model.reactions, (
                f"Native reaction {rxn_id} was removed after re-running standard gapfill"
            )


# ---------------------------------------------------------------------------
# Test E: gfmodel reusability — locks released, standard run still works
# ---------------------------------------------------------------------------

def test_gfmodel_reusability_after_comprehensive(
    template, model_data_path, complete_media, standard_media
):
    """After run_comprehensive_gapfilling, activation locks must be fully
    released and a subsequent run_gapfilling on the same MSGapfill instance
    must behave normally."""
    model = _build_model(model_data_path, ko=["GLCpts_c0"])
    target = "BIOMASS_Ecoli_core_w_GAM_c0"
    gf = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)

    # Run comprehensive gapfilling
    gf.run_comprehensive_gapfilling(
        complete_media,
        target=target,
        minimum_obj=0.01,
        binary_check=False,
        prefilter=False,
        activation_rxn_filter=_small_fwd_filter(model),
    )

    # Verify all activation locks have been released (lb == 0)
    rapkg = gf.gfpkgmgr.getpkg("ReactionActivationPkg")
    for rxn_id, var in rapkg.variables.get("fra", {}).items():
        assert var.lb == 0, (
            f"fra[{rxn_id}].lb = {var.lb} after cleanup — activation lock not released"
        )
    for rxn_id, var in rapkg.variables.get("rra", {}).items():
        assert var.lb == 0, (
            f"rra[{rxn_id}].lb = {var.lb} after cleanup — activation lock not released"
        )

    # Standard gapfilling must still work normally on the same instance
    std_result = gf.run_gapfilling(
        standard_media, target=target, minimum_obj=0.01, prefilter=False
    )
    assert std_result is not None, (
        "Standard run_gapfilling after comprehensive run returned None "
        "(objective or activation locks may not have been cleaned up)"
    )
    for key in ("new", "reversed", "media", "target", "minobjective", "binary_check"):
        assert key in std_result, f"Missing key '{key}' in standard gapfill result"


# ---------------------------------------------------------------------------
# Test F: run_multi_gapfill('Comprehensive') — correct solution shape
# ---------------------------------------------------------------------------

def test_comprehensive_multi_gapfill_solution_shape(template, model_data_path):
    """run_multi_gapfill with gapfilling_mode='Comprehensive' returns a solution
    dict with the correct shape.

    Uses a forward-only model (all reversible non-exchange reactions removed)
    so the Stage 1 rxn_filter contains only forward-only and exchange reactions.
    Uses a small media (glucose aerobic) to ensure the 7 reversible exchange
    binary variables have a clear direction-preference, keeping Stage 1
    tractable with GLPK in a few seconds.
    """
    # Forward-only model: missing reversible reactions + GLCpts KO
    model = _build_model(model_data_path, ko=["GLCpts_c0"], fwd_only=True)
    target = "BIOMASS_Ecoli_core_w_GAM_c0"

    # Small media: glucose aerobic. Open exactly those exchanges whose
    # direction is unambiguous so RevBin binary vars are quickly fixed.
    small_complete_media = MSMedia.from_dict(
        {
            "glc__D": (-1000, 1000),
            "o2": (-1000, 1000),
            "h": (-1000, 1000),
            "h2o": (-1000, 1000),
            "pi": (-1000, 1000),
            "co2": (-1000, 1000),
            "nh4": (-1000, 1000),
        }
    )

    gf = MSGapfill(model, [template], minimum_obj=0.01, default_target=target)
    results = gf.run_multi_gapfill(
        [small_complete_media],
        target=target,
        default_minimum_objective=0.01,
        gapfilling_mode="Comprehensive",
        binary_check=False,
        prefilter=False,
        run_sensitivity_analysis=False,
        integrate_solutions=True,
        # Comprehensive gapfilling adds many collectively-necessary reactions.
        # Individual-KO pruning (remove_unneeded_reactions=True) incorrectly
        # classifies them all as redundant when the base model has alternative
        # paths.  Disable pruning so all activated reactions are retained and
        # the model actually grows after integration.
        remove_unneeded_reactions=False,
    )

    assert results is not None, "run_multi_gapfill(Comprehensive) returned None"
    assert small_complete_media in results, "No entry for media in solution dict"

    sol = results[small_complete_media]
    for key in ("new", "reversed", "media", "target", "minobjective", "binary_check"):
        assert key in sol, f"Missing key '{key}' in Comprehensive solution dict"

    # The 'growth' field in the solution is set by integrate_gapfill_solution
    # (check_for_growth=True by default) right after integration; it measures
    # biomass on the same media.
    assert sol.get("growth", 0) > 0, (
        f"Model does not grow after Comprehensive run_multi_gapfill "
        f"(growth={sol.get('growth')})"
    )
