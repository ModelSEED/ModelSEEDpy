# -*- coding: utf-8 -*-
import pytest
from cobra.core import Metabolite, Model, Reaction

from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager


@pytest.fixture
def model():
    """One internal reaction plus one exchange, objective on the internal reaction."""
    m = Model("boundtest")
    cpd = Metabolite("cpd00027_e0", compartment="e0")
    ex = Reaction("EX_cpd00027_e0", lower_bound=-10, upper_bound=10)
    ex.add_metabolites({cpd: -1})
    rxn = Reaction("rxn00001_c0", lower_bound=0, upper_bound=1000)
    rxn.add_metabolites({cpd: 1})
    m.add_reactions([ex, rxn])
    m.objective = "rxn00001_c0"
    return m


def _objc_constraints(model):
    return [c for c in model.solver.constraints if c.name.endswith("_objc")]


def test_objconst_reuses_constraint_with_identical_coefficients(model):
    """
    Building the same objective constraint twice must update the existing
    constraint's bounds in place rather than adding a second one.
    """
    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("ObjConstPkg", True)

    first = pkg.build_package(0.1, 1000)
    assert len(_objc_constraints(model)) == 1

    second = pkg.build_package(0.5, 900)
    assert len(_objc_constraints(model)) == 1
    assert second is first
    assert second.lb == 0.5
    assert second.ub == 900


def test_objconst_build_package_returns_the_constraint(model):
    """build_package must return the constraint so the caller can adjust it."""
    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("ObjConstPkg", True)
    constraint = pkg.build_package(0.0, 1000)

    assert constraint is not None
    assert constraint.name.endswith("_objc")


def test_objconst_new_constraint_when_objective_changes(model):
    """A different objective has different coefficients, so it gets its own constraint."""
    mgr = MSPackageManager.get_pkg_mgr(model)
    pkg = mgr.getpkg("ObjConstPkg", True)
    pkg.build_package(0.1, 1000)

    model.objective = "EX_cpd00027_e0"
    pkg.build_package(0.2, 500)

    assert len(_objc_constraints(model)) == 2


def test_media_bounds_applied_when_uptake_exceeds_current_upper_bound(model):
    """
    Setting a lower bound more negative than the current upper bound must not leave
    the reaction transiently inverted: the upper bound is raised first.
    """
    # Upper bound sits BELOW the lower bound the media will apply (-100), which is the
    # case the ordering fix exists for: assigning the lower bound first would leave
    # lb > ub.
    model.reactions.EX_cpd00027_e0.lower_bound = -300
    model.reactions.EX_cpd00027_e0.upper_bound = -200

    mdlutl = MSModelUtil.get(model)
    pkg = mdlutl.pkgmgr.getpkg("KBaseMediaPkg", True)
    pkg.build_package(None, default_uptake=100, default_excretion=100)

    exchange = model.reactions.EX_cpd00027_e0
    assert exchange.lower_bound == -100
    assert exchange.upper_bound == 100
    assert exchange.lower_bound <= exchange.upper_bound


def test_media_package_records_current_media(model):
    """current_media is exposed so later packages can see what media is applied."""
    mdlutl = MSModelUtil.get(model)
    pkg = mdlutl.pkgmgr.getpkg("KBaseMediaPkg", True)

    assert pkg.current_media is None
    pkg.build_package(None, default_uptake=100, default_excretion=100)
    assert pkg.current_media is None  # no media object was supplied
