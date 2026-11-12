# -*- coding: utf-8 -*-
import pytest
from cobra.core import Metabolite, Model, Reaction

from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager


def _two_reaction_model(name):
    """One internal reaction plus one exchange on a shared compound."""
    m = Model(name)
    cpd = Metabolite("cpd00027_e0", compartment="e0")
    ex = Reaction("EX_cpd00027_e0", lower_bound=-10, upper_bound=10)
    ex.add_metabolites({cpd: -1})
    rxn = Reaction("rxn00001_c0", lower_bound=-7, upper_bound=13)
    rxn.add_metabolites({cpd: 1})
    m.add_reactions([ex, rxn])
    return m


@pytest.fixture
def host():
    return _two_reaction_model("host")


@pytest.fixture
def other():
    return _two_reaction_model("other")


def test_replicated_variables_belong_to_the_host_solver_interface(host, other):
    """
    The replicated variables must be built in the HOST model's solver interface.

    ``optlang.Variable`` is bound to whichever interface optlang selected as its
    default at import time, so ``Variable.clone`` produced a variable of that
    interface regardless of the host model's solver. Adding it to a host using
    any other solver raised::

        TypeError: Cannot add Variable ... of interface type
        optlang.glpk_interface to model of type optlang.glpk_exact_interface

    which made ``build_package`` unusable for every model whose solver is not
    optlang's default.
    """
    host.solver = "glpk_exact"
    pkg = MSPackageManager.get_pkg_mgr(host).getpkg("ProblemReplicationPkg", True)

    pkg.build_package({"models": [other]})
    host.solver.update()

    expected = host.problem.Variable
    replicated = pkg.variables["0"]
    assert replicated, "no variables were replicated"
    for var in replicated.values():
        assert isinstance(var, expected)
        assert host.variables.get(var.name, None) is not None


def test_replicated_variables_keep_their_bounds_and_type(host, other):
    """
    Bounds, type and the suffixed name must carry across to the replica.

    Run against a non-default solver so this fails on the unpatched code for the
    same interface reason as the test above: on optlang's default solver the old
    ``Variable.clone`` call preserved bounds too, so a test run there would pass
    either way and pin nothing.
    """
    host.solver = "glpk_exact"
    pkg = MSPackageManager.get_pkg_mgr(host).getpkg("ProblemReplicationPkg", True)

    pkg.build_package({"models": [other]})
    host.solver.update()

    source = {v.name: v for v in other.variables}
    for original_name, replica in pkg.variables["0"].items():
        original = source[original_name]
        assert replica.name == original_name + ".0"
        assert replica.lb == original.lb
        assert replica.ub == original.ub
        assert replica.type == original.type


def test_shared_variable_packages_are_resolved_by_package_name(host, other):
    """
    ``shared_variable_packages`` is keyed by package NAME, the way every other
    package lookup in the codebase is (``pkgmgr.getpkg("KBaseMediaPkg")``).

    The keys were used directly as package objects (``pkg.variables``), which on
    a ``str`` is a guaranteed ``AttributeError``, so this path could never run.
    Resolving the name through the package manager is what makes it work, and a
    shared variable must then NOT be replicated: the replicas share the host's.
    """
    host_tf = MSPackageManager.get_pkg_mgr(host).getpkg("TotalFluxPkg", True)
    host_tf.build_package()
    MSPackageManager.get_pkg_mgr(other).getpkg("TotalFluxPkg", True).build_package()
    host.solver.update()

    shared_name = host_tf.variables["tf"]["rxn00001_c0"].name
    assert shared_name == "rxn00001_c0_tf"

    pkg = MSPackageManager.get_pkg_mgr(host).getpkg("ProblemReplicationPkg", True)
    pkg.build_package(
        {"models": [other], "shared_variable_packages": {"TotalFluxPkg": ["tf"]}}
    )
    host.solver.update()

    # The shared variable is NOT copied ...
    assert shared_name not in pkg.variables["0"]
    assert host.variables.get(shared_name + ".0", None) is None
    # ... while an unshared one is.
    assert "rxn00001_c0" in pkg.variables["0"]
    assert host.variables.get("rxn00001_c0.0", None) is not None


def test_shared_variables_are_substituted_into_replicated_constraints(host, other):
    """
    A replicated constraint that referenced a shared variable must reference the
    HOST's variable object, not a copy -- that sharing is the point of the
    package.
    """
    host_tf = MSPackageManager.get_pkg_mgr(host).getpkg("TotalFluxPkg", True)
    host_tf.build_package()
    MSPackageManager.get_pkg_mgr(other).getpkg("TotalFluxPkg", True).build_package()
    host.solver.update()

    shared_var = host_tf.variables["tf"]["rxn00001_c0"]

    pkg = MSPackageManager.get_pkg_mgr(host).getpkg("ProblemReplicationPkg", True)
    pkg.build_package(
        {"models": [other], "shared_variable_packages": {"TotalFluxPkg": ["tf"]}}
    )
    host.solver.update()

    replicated = pkg.constraints["0"]["rxn00001_c0_tf"]
    assert any(v is shared_var for v in replicated.variables)
