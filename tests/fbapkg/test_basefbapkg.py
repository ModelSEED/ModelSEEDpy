# -*- coding: utf-8 -*-
import pytest
from cobra.core import Metabolite, Model, Reaction

from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.fbapkg.mspackagemanager import MSPackageManager


@pytest.fixture
def model():
    """One internal reaction plus one exchange on a shared compound."""
    m = Model("cleartest")
    cpd = Metabolite("cpd00027_e0", compartment="e0")
    ex = Reaction("EX_cpd00027_e0", lower_bound=-10, upper_bound=10)
    ex.add_metabolites({cpd: -1})
    rxn = Reaction("rxn00001_c0", lower_bound=0, upper_bound=1000)
    rxn.add_metabolites({cpd: 1})
    m.add_reactions([ex, rxn])
    return m


def test_clear_keeps_the_registered_type_keys(model):
    """
    ``clear`` replaced self.variables / self.constraints with ``{}``, dropping
    the per-type keys that ``__init__`` registered. Every build method indexes
    those keys directly -- ``build_variable`` does
    ``len(self.variables[obj_type])`` -- so rebuilding after a clear raised
    ``KeyError: 'tf'`` and a package could not be reused.
    """
    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("TotalFluxPkg", True)
    pkg.build_package()
    model.solver.update()
    assert len(pkg.variables["tf"]) == 2

    pkg.clear()

    # The type keys survive, now empty ...
    assert set(pkg.variables) == {"tf"}
    assert set(pkg.constraints) == {"tf"}
    assert pkg.variables["tf"] == {}
    assert pkg.constraints["tf"] == {}

    # ... so the package is reusable.
    pkg.build_package()
    model.solver.update()
    assert len(pkg.variables["tf"]) == 2


def test_clear_removes_the_objects_from_the_solver(model):
    """clear must still detach what it built from the solver."""
    pkg = MSPackageManager.get_pkg_mgr(model).getpkg("TotalFluxPkg", True)
    pkg.build_package()
    model.solver.update()
    built = [v.name for v in pkg.variables["tf"].values()]
    assert all(model.variables.get(n, None) is not None for n in built)

    pkg.clear()
    model.solver.update()

    for name in built:
        assert model.variables.get(name, None) is None


def test_get_current_media_reports_the_applied_media(model):
    """
    Any package can ask what media is in force without re-deriving it, instead
    of reaching into the package manager itself.
    """
    mdlutl = MSModelUtil.get(model)
    pkg = mdlutl.pkgmgr.getpkg("TotalFluxPkg", True)

    assert pkg.get_current_media() is None

    media_pkg = mdlutl.pkgmgr.getpkg("KBaseMediaPkg", True)
    media_pkg.build_package(None, default_uptake=100, default_excretion=100)

    # No media object was supplied, so there is still none to report, but the
    # accessor resolves to KBaseMediaPkg's attribute rather than raising.
    assert pkg.get_current_media() is media_pkg.current_media


def test_get_current_media_is_not_shadowed_on_the_media_package(model):
    """
    ``KBaseMediaPkg`` stores the media in an instance attribute called
    ``current_media``. Naming this accessor ``current_media`` would have been
    shadowed by that attribute on that subclass, so ``pkg.current_media()``
    would raise ``TypeError: 'NoneType' object is not callable``. The accessor
    must stay callable on every package, the media package included.
    """
    mdlutl = MSModelUtil.get(model)
    media_pkg = mdlutl.pkgmgr.getpkg("KBaseMediaPkg", True)

    assert callable(media_pkg.get_current_media)
    assert media_pkg.get_current_media() is media_pkg.current_media
