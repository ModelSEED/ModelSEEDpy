import inspect

import modelseedpy
from modelseedpy import fbapkg
from modelseedpy.fbapkg import BaseFBAPkg


def _fbapkg_package_classes():
    """Every concrete FBA package class exported by modelseedpy.fbapkg."""
    return {
        name: obj
        for name, obj in inspect.getmembers(fbapkg, inspect.isclass)
        if name.endswith("Pkg") and name != "BaseFBAPkg" and issubclass(obj, BaseFBAPkg)
    }


def test_every_fbapkg_package_is_exported_at_top_level():
    # MSPackageManager discovers packages from the TOP-LEVEL modelseedpy
    # namespace, so a package exported only from modelseedpy.fbapkg cannot be
    # found by getpkg() (e.g. getpkg("ObjectivePkg") raised PackageError).
    missing = sorted(
        name for name in _fbapkg_package_classes() if not hasattr(modelseedpy, name)
    )
    assert missing == []


def test_objectivepkg_is_available_to_package_manager():
    from modelseedpy.fbapkg.mspackagemanager import MSPackageManager

    class _Model:  # MSPackageManager.__init__ only stores the model
        pass

    manager = MSPackageManager(_Model())
    assert "ObjectivePkg" in manager.list_available_packages()
