import importlib


def test_import_combaero_package() -> None:
    """Basic smoke test that the top-level combaero package is importable."""
    mod = importlib.import_module("combaero")
    assert mod is not None


def test_import_core_module() -> None:
    """Basic smoke test that the internal _core extension can be imported."""
    core = importlib.import_module("combaero._core")
    assert core is not None


def test_no_module_imports_a_dev_only_dependency() -> None:
    """Every module of the package imports without sympy (a dev-only
    dependency; the MSVC smoke job installs numpy and scipy only). Run in a fresh interpreter so another test's import of sympy
    cannot hide a leak."""
    import subprocess
    import sys

    probe = (
        "import importlib, pkgutil, sys\n"
        "import combaero\n"
        "for m in pkgutil.walk_packages(combaero.__path__, 'combaero.'):\n"
        "    importlib.import_module(m.name)\n"
        "print(sorted(n for n in ('sympy', 'mpmath') if n in sys.modules))\n"
    )
    out = subprocess.run(
        [sys.executable, "-c", probe], capture_output=True, text=True, check=True
    ).stdout.strip()
    assert out == "[]", f"dev-only modules imported by the package: {out}"


def test_the_input_structs_pickle() -> None:
    """The member-count guard of the pickle helper (#1164) holds on this
    compiler: the MSVC smoke job runs this file only."""
    import copy
    import pickle

    import combaero as cb

    names = [
        "OrificeGeometry", "RibGeometry", "PinFinGeometry", "RibCorrelationSet",
        "RibRatioSet", "JetArrayCorrelationSet", "SingleJetImpingementSet",
        "PinFinNuSet", "PinFinFrictionSet", "PinFinRatioModifier", "RibTerm",
        "RibRange", "StatedAccuracy", "RatioBaseline", "ImpingementRange",
        "JetArrayGeometricFit",
    ]  # fmt: skip
    objs = [getattr(cb, n)() for n in names] + [cb.JetRowGeometry(1.0, 2.0, 3.0, 4.0, 5.0)]
    for o in objs:
        for twin in (pickle.loads(pickle.dumps(o)), copy.deepcopy(o)):
            assert type(twin) is type(o)
            assert repr(sorted(twin.__getstate__())) == repr(sorted(o.__getstate__()))
