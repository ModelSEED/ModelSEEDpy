# -*- coding: utf-8 -*-
import pytest
from cobra.core import Metabolite, Model, Reaction
from optlang.exceptions import ContainerAlreadyContains

from modelseedpy.fbapkg.basefbapkg import BaseFBAPkg
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager


def _model(name, extra_reaction_ids=()):
    m = Model(name)
    cpd = Metabolite("cpd00027_e0", compartment="e0")
    ex = Reaction("EX_cpd00027_e0", lower_bound=-10, upper_bound=10)
    ex.add_metabolites({cpd: -1})
    rxns = [ex]
    for rid in ("rxn00001_c0",) + tuple(extra_reaction_ids):
        r = Reaction(rid, lower_bound=0, upper_bound=1000)
        r.add_metabolites({cpd: 1})
        rxns.append(r)
    m.add_reactions(rxns)
    m.solver.update()
    return m


def _with_totalflux(name):
    """A model whose solver already carries a built TotalFluxPkg."""
    m = _model(name)
    pkg = MSPackageManager.get_pkg_mgr(m).getpkg("TotalFluxPkg", True)
    pkg.build_package()
    m.solver.update()
    return m, pkg


def test_adoption_is_off_by_default():
    """
    Default behaviour is unchanged: a newly constructed package claims nothing
    that is already on the solver. (Characterises existing behaviour -- it holds
    before and after this change.)
    """
    model = _model("default_off", extra_reaction_ids=("rxn00042_tf",))

    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("TotalFluxPkg", True)

    assert pkg.variables["tf"] == {}
    assert pkg.constraints["tf"] == {}


def test_rebuilding_a_package_on_a_copied_model_fails_without_adoption():
    """
    The failure this capability exists for, pinned so the motivation is visible.

    ``Model.copy()`` carries a package's variables and constraints across, but
    the copy gets its own package manager and therefore a FRESH, empty package
    object. Rebuilding on the copy re-adds names the solver already holds::

        ContainerAlreadyContains: Container ... already contains an object with
        name 'EX_cpd00027_e0_tf'

    (Characterises existing behaviour -- it holds before and after this change.)
    """
    model, _ = _with_totalflux("copy_fails")
    copied = model.copy()
    copied.solver.update()

    fresh = MSPackageManager.get_pkg_mgr(copied).getpkg("TotalFluxPkg", True)
    assert fresh.variables["tf"] == {}

    with pytest.raises(ContainerAlreadyContains):
        fresh.build_package()


def test_adopt_existing_objects_rebinds_a_package_on_a_copied_model():
    """
    Adoption rebinds the fresh package to what the copy already carries, so the
    package is usable again and nothing is rebuilt.
    """
    model, original = _with_totalflux("copy_adopt")
    expected = set(original.variables["tf"])

    copied = model.copy()
    copied.solver.update()
    fresh = MSPackageManager.get_pkg_mgr(copied).getpkg("TotalFluxPkg", True)

    fresh.adopt_existing_objects()

    assert set(fresh.variables["tf"]) == expected
    for key, var in fresh.variables["tf"].items():
        assert var.name == key + "_tf"
        # the copy's own objects, not the original model's
        assert var is copied.variables.get(var.name)
    # and building again is now a no-op rather than a collision
    fresh.build_package()
    copied.solver.update()
    assert set(fresh.variables["tf"]) == expected


def test_adopt_existing_can_be_requested_from_the_constructor():
    """``adopt_existing=True`` runs the scan during __init__."""
    model, original = _with_totalflux("ctor")
    copied = model.copy()
    copied.solver.update()

    pkg = BaseFBAPkg(
        copied,
        "totalflux probe",
        {"tf": "reaction"},
        {"tf": "reaction"},
        adopt_existing=True,
    )

    assert set(pkg.variables["tf"]) == set(original.variables["tf"])


def test_adoption_matches_on_name_only_and_can_claim_a_foreign_variable():
    """
    The hazard that makes this opt-in, pinned so it cannot be forgotten.

    The match is the final underscore-delimited segment of the name, so an
    ORDINARY reaction whose id ends in ``_tf`` is adopted as the total-flux
    variable for the reaction named by the remaining prefix. ``build_package``
    then finds that key taken and never builds the real variable or its
    coupling constraint -- silently, with no error.
    """
    model = _model("hazard", extra_reaction_ids=("rxn00042", "rxn00042_tf"))
    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("TotalFluxPkg", True)

    pkg.adopt_existing_objects()

    foreign = model.reactions.get_by_id("rxn00042_tf").forward_variable
    assert pkg.variables["tf"]["rxn00042"] is foreign

    pkg.build_package()
    model.solver.update()

    assert pkg.variables["tf"]["rxn00042"] is foreign
    assert pkg.variables["tf"]["rxn00042"].ub == 1000  # not the tf bound of 100
    assert "rxn00042" not in pkg.constraints["tf"]
