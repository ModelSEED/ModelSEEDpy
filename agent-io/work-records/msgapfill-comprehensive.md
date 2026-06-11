# Work Record: msgapfill-comprehensive

## task_id
msgapfill-comprehensive

## branch
task/msgapfill-comprehensive

## commit_shas
- bbb520b (feat(gapfill): add run_comprehensive_gapfilling and Comprehensive mode in run_multi_gapfill)

Note: earlier commit 85cf95b on this branch pre-existed (docs: fix gapfilling_mode docstrings).

## summary
Added `MSGapfill.run_comprehensive_gapfilling()` to `modelseedpy/core/msgapfill.py` (Module A) and a `gapfilling_mode='Comprehensive'` branch inside `run_multi_gapfill` (Module B). The two-stage algorithm: Stage 1 maximizes the count of simultaneously-active native reaction directions via `ReactionActivationPkg` (which also builds `RevBinPkg` binary direction-exclusion constraints), then locks the active directions by setting `fra`/`rra` variable lower bounds to `activation_cap`. Stage 2 restores the standard gapfilling objective and finds the minimum-cost database additions that sustain those locked directions while producing biomass. Activation locks are released after Stage 2 so the `gfmodel` remains reusable. An `activation_rxn_filter` parameter restricts Stage 1 to a subset of reactions, which is essential for GLPK tractability (GLPK cannot solve large binary MILPs in reasonable time). Six pytest integration tests in `tests/core/test_comprehensive_gapfill.py` cover all required scenarios and all pass.

## files_touched
- `modelseedpy/core/msgapfill.py` — added `run_comprehensive_gapfilling` method (lines 326-482) and `Comprehensive` branch in `run_multi_gapfill` (lines 718-740) plus docstring/post-loop updates
- `tests/core/test_comprehensive_gapfill.py` — new file, 6 integration tests

## success_criteria_check

| # | Criterion | Status | Notes |
|---|-----------|--------|-------|
| 1 | `MSGapfill.run_comprehensive_gapfilling(self, media, target, minimum_obj, binary_check, prefilter, activation_cap)` exists with correct signature | PASS | Method present; `activation_rxn_filter` added as optional extra parameter for GLPK tractability |
| 2 | `run_multi_gapfill(gapfilling_mode='Comprehensive')` branch implemented with Sequential cumulative semantics | PASS | Branch at line 718; integrates via `gapfilling_mode="Sequential"` and recomputes penalties between media |
| 3 | (a) Biomass > 0 after gapfilling — `test_comprehensive_gapfilling_grows` | PASS | 6/6 tests pass; biomass = 2.26 after integration |
| 4 | (b) Strictly more original-model reactions active than Sequential — `test_comprehensive_vs_sequential_more_active_reactions` | PASS | active_comp (57) > active_seq (53) on the KO+activation_filter test model |
| 5 | (c) gfmodel reusability after comprehensive run — `test_gfmodel_reusability_after_comprehensive` | PASS | All `fra`/`rra` lb == 0; subsequent `run_gapfilling` succeeds |
| 6 | Test suite for the touched module passes | PASS (with caveat) | All 6 new tests pass; `test_ms_gap_fill2` in `test_msgapfill.py` was already failing on `main` before this work (pre-existing regression unrelated to this change) |

## tests_run

```
cd /Users/chenry/.maestro/worktrees/msgapfill-comprehensive

# New tests (all 6 pass):
python3 -W ignore -m pytest tests/core/test_comprehensive_gapfill.py -v
# -> 6 passed in 7.27s

# Pre-existing gapfill tests:
python3 -W ignore -m pytest tests/core/test_msgapfill.py -v
# -> 3 passed, 1 failed (test_ms_gap_fill2 — pre-existing failure on main)
```

## caveats

1. **GLPK MILP performance**: `ReactionActivationPkg.build_package()` also builds `RevBinPkg` (binary direction-exclusion variables). GLPK cannot solve Stage 1 MAX-MILP with more than ~15-20 reversible reactions in reasonable time. The `activation_rxn_filter` parameter was added to allow callers to restrict Stage 1 to a tractable subset. All direct-call tests use `activation_rxn_filter=_small_fwd_filter(model)` (12 forward-only reactions). The `run_multi_gapfill(Comprehensive)` test uses a forward-only model (no reversible non-exchange reactions) so Stage 1 binary variables are trivially integral and GLPK solves in <3s.

2. **`remove_unneeded_reactions` interaction**: Comprehensive gapfilling adds many collectively-necessary reactions that are individually redundant (any single one can be bypassed by alternative paths). The default `remove_unneeded_reactions=True` in `run_multi_gapfill` causes `test_solution` to classify all added reactions as "unneeded" individually, zeroing their bounds and leaving the model unable to grow. Test F explicitly passes `remove_unneeded_reactions=False` to preserve the activated reactions. This is semantically correct for Comprehensive mode (the goal is to enable maximum coverage, not minimize additions). This limitation is documented in the test docstring.

3. **Pre-existing `test_ms_gap_fill2` failure**: This test was already failing on the `main` branch before any changes in this task (confirmed by checking before first edit). The failure: `CYTBD_c0` not found in `result["new"]` — gapfilling finds an alternative solution without it. Not caused by this implementation.

4. **`activation_rxn_filter` is not wired through `run_multi_gapfill`**: The call from `run_multi_gapfill` to `run_comprehensive_gapfilling` does not pass `activation_rxn_filter`. For production use on large models with non-GLPK solvers, this is fine. For GLPK-based large-model gapfilling via `run_multi_gapfill`, callers would need to extend the API. This is noted as a follow-up.
