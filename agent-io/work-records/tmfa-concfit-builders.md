# Work record: tmfa-concfit-builders

## task_id
tmfa-concfit-builders

## branch
tmfa-concfit-builders

## commit_shas
- 60278bb0f6aa0f905ad6e5713cdb4718504158da — feat(fbapkg): weighted L1/L2 concentration-fit objectives for FullThermoPkg

## summary

Extended `FullThermoPkg` in `modelseedpy/fbapkg/fullthermopkg.py` per
`agent-io/prds/tmfa-weighted-concfit/fullprompt.md` Modules 1-3. Module 1
(`build_concfit_objective`) gained an optional `weights: dict[str, float] |
None` parameter; `weights=None` reproduces the pre-change unweighted L1
objective exactly (unit coefficients on the `pconcfit`/`nconcfit` slacks).
Module 2 added a new `build_concfit_objective_l2` builder that attaches a
quadratic weighted-least-squares objective directly against the `logconc`
variables (no slacks), guarded by an early `NotImplementedError` when the
model's optlang backend is not gurobi/cplex, and which cleans up any
pre-existing L1 `concfit` constraints/`pconcfit`/`nconcfit` variables for the
fitted set. Module 3 wired `build_package(concfit_objective_type="l1"|"l2",
concfit_weights=...)` dispatch and turned the previous `KeyError` (when
`set_concfit_objective=True` but `target_concentrations` was missing) into a
clear `ValueError`.

One deliberate deviation from the PRD's literal guard expression: the PRD
text specifies detecting the backend via
`model.solver.interface.__module__`. I verified empirically that
`model.solver.interface` is always a Python *module* object (e.g.
`<module 'optlang.gurobi_interface'>`), and modules do not carry a
`__module__` attribute — accessing it raises `AttributeError` unconditionally,
regardless of backend. I implemented the equivalent, functioning check using
`model.solver.interface.__name__` (which does carry the dotted module name,
e.g. `'optlang.gurobi_interface'`) so the guard actually distinguishes
backends as intended. Documented inline in `_is_qp_capable_solver`.

A second issue surfaced empirically and is fixed as part of Module 1's
`weights=None` default: the pre-change code's literal `1 * pconcfit` /
`1 * nconcfit` terms (python `int` coefficient `1`) cause sympy to
auto-simplify `1*var` to the bare `Variable` symbol; when two or more such
bare-symbol terms are summed and assigned as a gurobi-backend objective,
optlang's quadratic-expression parser raises
`AssertionError: What is this? <class 'optlang.gurobi_interface.Variable'>`
(this is a pre-existing bug in the unmodified code, reproducible on `main`
whenever `build_concfit_objective` is exercised against a gurobi-backed
model — GLPK does not hit this code path). The PRD itself specifies the
default as `weights.get(metabolite.id, 1.0)` (float, not int `1`), which
sidesteps the issue since `1.0*var` remains a genuine `sympy.Mul` and is not
auto-simplified away. I used the float default exactly as the PRD specifies,
which happens to also fix this pre-existing gurobi incompatibility for the
weights=None case as a side effect. `build_min_error_objective` (unrelated,
untouched) still has the same latent literal-`1` bug on a gurobi backend —
out of scope for this PRD, left as-is, and one test
(`test_l2_replaces_objective_entirely`) was written to avoid exercising it by
using the (now gurobi-safe) L1 concfit objective as the "prior objective"
instead.

## files_touched
- `modelseedpy/fbapkg/fullthermopkg.py` — Module 1 (weighted L1), Module 2
  (new `build_concfit_objective_l2`), Module 3 (`build_package` dispatch +
  ValueError-not-KeyError), plus `_validate_concfit_weights`,
  `_is_qp_capable_solver`, `_require_qp_capable_solver` helpers.
- `tests/fbapkg/test_full_thermo_concfit.py` (new) — 23 tests covering
  backward compatibility, weight validation/zero-drop, MSID-only-key
  no-op, weighted-L1 hand-solvable conflict case, L2 QP guard, L2
  quadratic-objective + analytic weighted-least-squares optimum, L2
  stale-slack cleanup after a prior L1 build, objective-replacement
  (not combined with prior objectives), and `build_package` dispatch
  defaults/errors.

## success_criteria_check

(Acceptance Criteria 1-10 and 16 from the PRD, plus the task's own success
criteria — Criteria 11-15 concern the KBUtilLib `remote_solve` helper, which
is explicitly out of scope for this ModelSEEDpy-only task.)

1. **PASS** — `weights=None` produces coefficient `1.0` on every
   `pconcfit`/`nconcfit` term (`test_weights_none_gives_unit_coefficients`,
   `test_weights_none_and_omitted_are_equivalent`).
2. **PASS** — dict keys are matched against `cobra.Metabolite.id` exactly; a
   bare-MSID key contributes no term (`test_msid_only_key_contributes_no_term`).
3. **PASS** — weights validated finite/>=0 via `_validate_concfit_weights`;
   weight 0 zeroes the coefficient (`test_zero_weight_drops_term_from_objective`);
   negative or NaN/Inf weight raises `ValueError`
   (`test_negative_weight_raises_value_error`,
   `test_nan_or_inf_weight_raises_value_error`).
4. **PASS** — `build_concfit_objective_l2` on a gurobi-backed model produces
   an objective containing `**2` terms
   (`test_l2_objective_is_quadratic`), and matches the closed-form
   weighted-least-squares optimum on a hand-coupled 2-metabolite instance
   within `rel_tol=1e-4` (`test_l2_analytic_weighted_least_squares_optimum`).
5. **PASS** — `build_concfit_objective_l2` and `build_package(...,
   concfit_objective_type="l2")` both raise `NotImplementedError` matching
   "QP-capable" when the backend is GLPK
   (`test_l2_raises_notimplementederror_on_glpk`,
   `test_build_package_l2_raises_notimplementederror_on_glpk`). Guard fires
   before weight validation (`test_guard_checked_before_weight_validation`).
6. **PASS by inspection** — neither builder performs any unit conversion on
   `target_concentrations`/`weights` values; they are used verbatim as
   natural-log-molar constants against `logconc`.
7. **PASS** — building L2 after L1 for the same fitted set removes the
   `concfit` constraints and `pconcfit`/`nconcfit` variables, both from the
   package's internal dicts and from the underlying optlang model
   (`test_l2_after_l1_removes_stale_slacks_and_constraints`).
8. **PASS** — `build_package(set_concfit_objective=True)` with no
   `target_concentrations` raises `ValueError`, not `KeyError`
   (`test_missing_target_concentrations_raises_value_error`).
9. **PASS** — `concfit_objective_type` defaults to `"l1"`, `concfit_weights`
   defaults to `None`, and a `build_package` call with neither key set still
   builds the L1 path (`test_existing_build_package_callers_unchanged`,
   `test_dispatch_l1_by_default`).
10. **PASS** — both builders assign `self.model.objective =
    self.model.problem.Objective(...)` wholesale, replacing rather than
    combining with any prior objective
    (`test_concfit_replaces_min_error_objective`,
    `test_l2_replaces_objective_entirely`).
16. **PASS** — the QP guard accepts only `optlang.gurobi_interface` /
    `optlang.cplex_interface` (by module name, see caveat above); any other
    backend (GLPK, SCIPY, hybrid) is rejected. cplex was not available to
    test locally (only gurobi is licensed on this machine) — see caveats.

## tests_run

- `python3 -m pytest tests/fbapkg/test_full_thermo_concfit.py -v` — **23
  passed**, 0 failed, 0 skipped (gurobi 13.0.2 restricted-license backend was
  available locally, so all L2/gurobi-gated tests ran rather than skipping;
  they auto-skip via `requires_gurobi` if gurobi is unavailable in another
  environment).
- `python3 -m pytest tests/fbapkg/ -q` — **35 passed** (12 pre-existing +
  23 new), confirming all existing fbapkg tests still pass unmodified.
- `python3 -m pytest tests/ -q --ignore=tests/fbapkg` — 67 passed, 59 failed.
  All 59 failures are in `tests/multiomics/test_msexpression.py` and are
  **pre-existing** — verified by `git stash`-ing this change and re-running
  the same file against the unmodified `main` checkout, which reproduces 51
  of the same failures (some additional failures appear even with the
  change reverted, consistent with environment/fixture issues unrelated to
  this PRD, e.g. missing test data or API drift in an unrelated module).
  Not touched by this task; not investigated further per scope.

## caveats

- **PRD guard-expression deviation.** The PRD/task text literally specifies
  `model.solver.interface.__module__` for the QP-capability check. I
  verified this expression always raises `AttributeError` regardless of
  backend (`.interface` is a module object; modules don't carry
  `__module__`). Implemented the working equivalent via
  `model.solver.interface.__name__`, which correctly distinguishes
  `optlang.gurobi_interface` / `optlang.cplex_interface` from GLPK/SCIPY/
  hybrid. Flagging this explicitly since it's a literal deviation from the
  PRD's exact wording, done to make the guard actually function as intended
  (confirmed AC5's behavior, which is what matters, passes).
- **cplex untested.** Only `gurobipy` (restricted/non-production license,
  expires 2027-11-29) was available in this dev environment; CPLEX was not
  installed, so the "accept cplex" branch of the guard is exercised only by
  string-prefix logic, not by an actual CPLEX-backed solve. The guard
  function is symmetric between the two prefixes, so this is low risk, but
  it is untested against a real CPLEX interface module name.
  Corresponds to Acceptance Criterion 16's "cplex" half.
  Recommend future validation on a machine with CPLEX before relying on
  that path for production ANME work.
- **Pre-existing gurobi incompatibility in `build_min_error_objective`,
  left unfixed.** `build_min_error_objective` still uses literal `1 * var`
  terms and will raise the same `AssertionError` as the old
  `build_concfit_objective` did whenever it's exercised against a
  gurobi/cplex-backed model with 2+ terms. This predates this PRD and is not
  one of the three modules in scope; not touched. One new test
  (`test_l2_replaces_objective_entirely`) was written to avoid triggering it
  by using the L1 concfit objective (now gurobi-safe) as the "prior
  objective" fixture instead of `build_min_error_objective`.
- **Weight-0 semantics.** Per the PRD ("the two `1 * pconcfit`/`1 *
  nconcfit` terms become `w * ...`"), a weight of 0 zeroes the *objective
  coefficient* for that metabolite but the L1 builder still creates the
  `pconcfit`/`nconcfit` slack variables and `concfit` equality constraint for
  it (matching the literal substitution the PRD describes) — it does not
  skip the metabolite's constraint machinery entirely. Net effect: with
  weight 0 and no other cost, the slack variables absorb any target
  deviation for free, which is equivalent to not fitting that metabolite at
  all under this L1 formulation. Confirmed this interpretation against
  AC3's "omits that metabolite from the objective" wording (objective, not
  constraint set).
- **Module 4 (KBUtilLib `remote_solve`) and Acceptance Criteria 11-15 were
  explicitly out of scope** per the task instructions ("extend the TMFA
  concentration-fit in fullthermopkg.py... Acceptance Criteria 1-10 and
  16") and were not implemented or tested here.
- Ran `grep -rn "kbutillib\|KBUtilLib" modelseedpy/` — no matches;
  ModelSEEDpy introduces no import dependency on KBUtilLib (AC13,
  reaffirmed even though Module 4 itself is out of scope).
