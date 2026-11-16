"""A resolved or solved network deep-copies (#1164).

Three things stopped it:
- write-only debug stores of C++ results on nodes (a node's wall couplings,
  a combustor's last mixer result), now gone;
- input structs held by elements (OrificeGeometry, JetRowGeometry, the
  correlation sets) without pickle support, now pickled from their
  read-write fields;
- the heat-transfer wrappers' forwarding __getattr__, which recursed on the
  blank instance copy builds, now refusing _r and dunder names.

A copy must be the same model: it re-solves to the same root, and its F and
J at the solution are bit-identical to the original's at matched evaluation
history. An evaluation leaves the Fanno kernels' warm starts (_last_G,
_last_M_exit) on the elements, so two back-to-back evaluations at one x can
differ in the last bits; here both start from the solve's own warm starts.
"""

from __future__ import annotations

import copy
import math
import pickle
import random
import struct
import sys
import warnings
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
    RibbedModel,
)
from combaero.network.components import (
    _ChamberHtc,
    _DuctHtc,
    _FilmedHtc,
    _FilmSideHtc,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
from networks import FAMILIES  # noqa: E402

# --- field-for-field comparison ---------------------------------------------


def _fields(o: Any) -> Any:
    """A comparable image of an object: pybind structs by their pickled
    state (recursively), enums by type and value, floats by their bits (NaN
    equals NaN, -0.0 differs from 0.0)."""
    if hasattr(o, "__getstate__") and type(type(o)).__name__ == "pybind11_type":
        if hasattr(type(o), "__members__"):  # a pybind enum
            return (type(o).__name__, int(o))
        state = o.__getstate__()
        assert isinstance(state, dict), type(o).__name__
        return (type(o).__name__, {k: _fields(v) for k, v in state.items()})
    if isinstance(o, float):
        return ("f", struct.pack("<d", o))
    if isinstance(o, (list, tuple)):
        return [_fields(v) for v in o]
    return o


# --- the network families -----------------------------------------------------


def _ribbed_channel() -> FlowNetwork:
    """A mass flow through a two-wall ribbed duct (Han 1988 orthogonal
    ribs) into a pressure boundary: the network holds a RibCorrelationSet."""
    T, P, D = 300.0, 2.0e5, 0.05
    Y = cb.mole_to_mass(cb.species.dry_air())
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("in", m_dot=0.05, Tt=T, Y=Y))
    out = PressureBoundary("out")
    out.Pt, out.Tt, out.Y = P, T, Y
    g.add_node(out)
    g.add_element(
        ChannelElement(
            id="ch",
            from_node="in",
            to_node="out",
            length=1.0,
            diameter=D,
            roughness=0.0,
            surface=ConvectiveSurface(
                area=0.1,
                model=RibbedModel(
                    correlation_set=cb.han_1988_orthogonal(),
                    e_D=0.063,
                    p_e=10.0,
                    W_H=1.0,
                    n_ribbed_walls=2,
                ),
            ),
        )
    )
    return g


# Each seed failed to deep-copy on main after a solve, through the type named.
CASES = {
    "hotcold": (lambda: FAMILIES["hotcold"](0)[0], "WallCouplingResult"),
    "wall": (lambda: FAMILIES["pin_fin"](1)[0], "WallCouplingResult"),
    "combustor": (lambda: FAMILIES["fuel_controller"](0)[0], "MixerResult"),
    "orifice": (lambda: FAMILIES["effusion"](0)[0], "OrificeGeometry"),
    "impingement": (lambda: FAMILIES["impingement"](0)[0], "JetArrayCorrelationSet"),
    "rib": (_ribbed_channel, "RibCorrelationSet"),
}


def _solve(net: FlowNetwork) -> dict[str, Any]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(net).solve(timeout=60)
    assert r["__success__"], r.get("__message__")
    return r


def _f_and_j(net: FlowNetwork, x: np.ndarray) -> tuple[bytes, bytes]:
    """F and dense J at x from a fresh solver, set up as solve() sets one up."""
    s = NetworkSolver(net)
    s._build_x0()
    F, J = s._residuals_and_jacobian(x)
    Jd = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
    return np.asarray(F, dtype=float).tobytes(), np.asarray(Jd, dtype=float).tobytes()


@pytest.mark.parametrize("case", sorted(CASES))
def test_a_network_deep_copies_before_and_after_solve(case: str) -> None:
    build, _ = CASES[case]
    net = build()
    copy.deepcopy(net)  # as built
    net.resolve_all_topology()
    NetworkSolver(net)._build_x0()
    copy.deepcopy(net)  # resolved
    _solve(net)
    copy.deepcopy(net)  # solved


@pytest.mark.parametrize("case", sorted(CASES))
def test_a_copied_solved_network_is_the_same_model(case: str) -> None:
    """Same F and J, bit for bit, at the solution (each the first evaluation
    after the solve, so from the same warm starts); the same root re-solved."""
    net = CASES[case][0]()
    r = _solve(net)
    x = np.asarray(r["__x_solution__"], dtype=float)
    twin = copy.deepcopy(net)

    F0, J0 = _f_and_j(net, x)
    F1, J1 = _f_and_j(twin, x)
    assert F1 == F0
    assert J1 == J0

    x_net = np.asarray(_solve(net)["__x_solution__"], dtype=float)
    x_twin = np.asarray(_solve(twin)["__x_solution__"], dtype=float)
    np.testing.assert_allclose(x_twin, x_net, rtol=1e-10, atol=0.0)


def test_the_combustor_keeps_no_mixer_result() -> None:
    """The stores are gone, not just made copyable: nothing read them."""
    net = CASES["combustor"][0]()
    _solve(net)
    for node in net.nodes.values():
        assert not hasattr(node, "_last_mix_res")
        assert not hasattr(node, "_wall_couplings")
    net = CASES["hotcold"][0]()
    _solve(net)
    for node in net.nodes.values():
        assert not hasattr(node, "_wall_couplings")


# --- the input structs pickle field for field ----------------------------------


def _orifice_geometry() -> Any:
    g = cb.OrificeGeometry()
    g.d, g.D, g.t, g.r, g.bevel = 0.02, 0.05, 0.003, 4e-4, 0.1
    g.tap = cb.OrificeTap.Corner
    return g


def _pin_fin_geometry() -> Any:
    return cb.PinFinGeometry(2.5, 2.5, 1.0, 8, cb.PinArrangement.Staggered)


def _rib_geometry() -> Any:
    return cb.RibGeometry(e_D=0.063, p_e=10.0, W_H=1.0, alpha_deg=90.0)


def _jet_row_geometry() -> Any:
    return cb.JetRowGeometry(d=2.54e-3, xn_d=5.0, yn_d=4.0, z_d=2.0, n_holes=12.0)


def _stated_accuracy() -> Any:
    return cb.StatedAccuracy.stated(0.06)


# (builder, how results are compared: a callable of the object -> numbers)
PICKLED: dict[str, tuple[Any, Any]] = {
    "OrificeGeometry": (
        _orifice_geometry,
        lambda g: cb.Cd_sharp_thin_plate_and_dRe(g, _orifice_state()),
    ),
    "JetArrayCorrelationSet": (
        cb.florschuetz_1981_staggered,
        lambda s: _jet(s),
    ),
    "JetRowGeometry": (
        _jet_row_geometry,
        lambda g: _row(g),
    ),
    "SingleJetImpingementSet": (
        cb.goldstein_1986_single_jet,
        lambda s: cb.single_jet_impingement_nu(
            s, cb.ImpingementThermalBC.ConstantHeatFlux, 2.3e4, 7.75, 5.0
        ),
    ),
    "RibCorrelationSet": (
        cb.han_park_1988_angled,
        lambda s: _rib(cb.evaluate_rib(s, cb.RibGeometry(0.0625, 10.0, 2.0, 60.0), 3e4)),
    ),
    "RibRatioSet": (
        lambda: cb.taslim_spring_1987(1.0, 0.083),
        lambda s: _rib(cb.evaluate_rib_ratio(s, cb.RibGeometry(0.083, 10.0, 1.0, 90.0), 3e4, 0.7)),
    ),
    "RibGeometry": (
        _rib_geometry,
        lambda g: _rib(cb.evaluate_rib(cb.han_1988_orthogonal(), g, 3e4)),
    ),
    "PinFinGeometry": (
        _pin_fin_geometry,
        lambda g: cb.evaluate_pin_fin_nu(cb.metzger_1986_staggered_nu(), g, 2e4, 0.7).Nu,
    ),
    "PinFinNuSet": (
        cb.metzger_1986_staggered_nu,
        lambda s: cb.evaluate_pin_fin_nu(s, _pin_fin_geometry(), 2e4, 0.7).Nu,
    ),
    "PinFinFrictionSet": (
        cb.metzger_1982_staggered_friction,
        lambda s: cb.evaluate_pin_fin_friction(s, _pin_fin_geometry(), 2e4).f,
    ),
    "PinFinRatioModifier": (
        cb.chyu_1998_inline_over_staggered,
        lambda s: _attrs(cb.evaluate_pin_fin_modifier(s, _inline_pins(), 2e4)),
    ),
    "StatedAccuracy": (_stated_accuracy, lambda a: (a.value, a.usable_as_band())),
    "RatioBaseline": (
        lambda: cb.RatioBaseline(0.023, 0.8, 0.4),
        lambda b: (b.coeff, b.re_exponent, b.pr_exponent),
    ),
    "RibTerm": (lambda: cb.RibTerm(0.35, 10.0), lambda t: (t.exponent, t.reference)),
    "RibRange": (lambda: cb.RibRange(1e4, 6e4), lambda r: (r.lo, r.hi)),
    "ImpingementRange": (lambda: cb.ImpingementRange(5.0, 15.0), lambda r: (r.lo, r.hi)),
    "JetArrayGeometricFit": (
        lambda: cb.JetArrayGeometricFit(0.1, 0.2, 0.3, 0.4),
        lambda f: (f.C, f.nx, f.ny, f.nz),
    ),
}


def _orifice_state() -> Any:
    st = cb.OrificeState()
    st.Re_D = 1e5
    return st


def _inline_pins() -> Any:
    return cb.PinFinGeometry(2.5, 2.5, 1.0, 8, cb.PinArrangement.Inline)


def _jet(s: Any) -> tuple[float, ...]:
    r = cb.jet_array_impingement_nu(s, 2e4, 0.2, 0.7, 6.0, 5.0, 2.0)
    return (r.Nu, r.dNu_dRe_j, r.dNu_dGc_Gj)


def _row(g: Any) -> tuple[float, ...]:
    return _attrs(
        cb.jet_row_heat_transfer(cb.florschuetz_1981_inline(), g, 0.01, 0.002, 2e-5, 0.03, 0.7)
    )


def _rib(r: Any) -> tuple[float, ...]:
    return _attrs(r)


def _attrs(r: Any) -> tuple[Any, ...]:
    names = sorted(n for n in dir(r) if not n.startswith("_"))
    vals = [getattr(r, n) for n in names]
    return tuple(v for v in vals if isinstance(v, (int, float, bool)))


def _same_numbers(a: Any, b: Any) -> bool:
    a, b = np.atleast_1d(np.asarray(a, dtype=float)), np.atleast_1d(np.asarray(b, dtype=float))
    return a.shape == b.shape and a.tobytes() == b.tobytes()


@pytest.mark.parametrize("name", sorted(PICKLED))
def test_a_copied_input_struct_gives_identical_results(name: str) -> None:
    """A realistic instance (the shipped factory, or a working geometry):
    the copy compares field for field and evaluates to the same bits."""
    build, evaluate = PICKLED[name]
    obj = build()
    assert type(obj).__name__ == name
    for twin in (pickle.loads(pickle.dumps(obj)), copy.deepcopy(obj), copy.copy(obj)):
        assert type(twin) is type(obj)
        assert twin is not obj
        assert _fields(twin) == _fields(obj)
        assert _same_numbers(evaluate(twin), evaluate(obj))


# --- every field away from its default ------------------------------------------
#
# A factory leaves many fields at their T{} default (83 of 269 leaves in the
# fixtures above), and a restore that skipped one of those would still
# compare equal. Here every leaf, nested ones included, is first set to a
# value its default does not have.


def _default(name: str) -> Any:
    """The T{} instance: every bound constructor default is the C++ one."""
    cls = getattr(cb, name)
    if name == "JetRowGeometry":  # no default constructor bound
        return cls(0.0, 0.0, 0.0, 0.0, 0.0)
    return cls()


def _rw_fields(cls: type) -> list[str]:
    out: list[str] = []
    for klass in cls.__mro__:
        for k, v in vars(klass).items():
            if isinstance(v, property) and v.fset is not None and k not in out:
                out.append(k)
    return out


def _is_struct(v: Any) -> bool:
    return type(type(v)).__name__ == "pybind11_type" and not hasattr(type(v), "__members__")


# Awkward floats a lossy copy would get wrong: NaN with a payload, -0.0,
# inf, a subnormal.
_ODD_FLOATS = (
    struct.unpack("<d", struct.pack("<Q", 0x7FF8000000000123))[0],
    -0.0,
    math.inf,
    5e-324,
)


def _different(v: Any, rng: random.Random) -> Any:
    if isinstance(v, bool):
        return not v
    if hasattr(type(v), "__members__"):  # a pybind enum
        others = [m for m in type(v).__members__.values() if m != v]
        assert others, f"{type(v).__name__} has one member"
        return rng.choice(others)
    if isinstance(v, int):
        return v + rng.randint(1, 99)
    if isinstance(v, float):
        while True:
            w = rng.choice((*_ODD_FLOATS, rng.uniform(-1e3, 1e3)))
            if struct.pack("<d", w) != struct.pack("<d", v):
                return w
    if isinstance(v, str):
        return v + "x\u00e9" + str(rng.random())
    raise AssertionError(f"unhandled field type {type(v).__name__}")


def _perturb(obj: Any, rng: random.Random) -> Any:
    """Set every read-write field, recursively, away from its value."""
    for k in _rw_fields(type(obj)):
        v = getattr(obj, k)
        if _is_struct(v):
            setattr(obj, k, _perturb(copy.deepcopy(v), rng))
        else:
            setattr(obj, k, _different(v, rng))
    return obj


def _leaves(obj: Any, pre: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    for k in _rw_fields(type(obj)):
        v = getattr(obj, k)
        if _is_struct(v):
            out.update(_leaves(v, f"{pre}{k}."))
        elif isinstance(v, float):
            out[pre + k] = ("f", struct.pack("<d", v))
        else:
            out[pre + k] = (type(v).__name__, repr(v))
    return out


def _set_leaf(obj: Any, path: str, value: Any) -> None:
    head, _, rest = path.partition(".")
    if not rest:
        setattr(obj, head, value)
        return
    inner = getattr(obj, head)
    _set_leaf(inner, rest, value)
    setattr(obj, head, inner)


def _get_leaf(obj: Any, path: str) -> Any:
    for part in path.split("."):
        obj = getattr(obj, part)
    return obj


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("name", sorted(PICKLED))
def test_every_field_round_trips(name: str, seed: int) -> None:
    """Every leaf away from its default, then pickle (protocols 2-5),
    deepcopy and copy: every leaf comes back, bit for bit."""
    rng = random.Random(f"{name}:{seed}")
    obj = _perturb(_default(name), rng)
    ref = _leaves(obj)
    base = _leaves(_default(name))
    assert ref.keys() == base.keys()
    assert not [k for k in ref if ref[k] == base[k]], "a leaf is still at its default"
    twins = [pickle.loads(pickle.dumps(obj, protocol=p)) for p in range(2, 6)]
    twins += [copy.deepcopy(obj), copy.copy(obj)]
    for twin in twins:
        assert type(twin) is type(obj)
        assert _leaves(twin) == ref


@pytest.mark.parametrize("name", sorted(PICKLED))
def test_a_restore_that_skips_any_one_field_is_caught(name: str) -> None:
    """The comparison above sees every field: a copy with any single leaf
    left at its default (a __setstate__ that skipped it) compares unequal."""
    rng = random.Random(name)
    obj = _perturb(_default(name), rng)
    ref = _leaves(obj)
    default = _default(name)
    for path in ref:
        lossy = copy.deepcopy(obj)
        _set_leaf(lossy, path, _get_leaf(default, path))
        assert _leaves(lossy) != ref, path


@pytest.mark.parametrize("name", sorted(PICKLED))
def test_a_state_of_the_wrong_size_is_refused(name: str) -> None:
    """The pickle helper knows each struct's C++ member count (asserted at
    compile time against the struct itself, not its bindings). A state with
    a field missing, or one too many, raises instead of leaving a member at
    its default."""
    obj = PICKLED[name][0]()
    state = obj.__getstate__()
    missing = dict(state)
    missing.pop(next(iter(missing)))
    blank = type(obj).__new__(type(obj))
    with pytest.raises(RuntimeError, match="C\\+\\+ members"):
        blank.__setstate__(missing)
    extra = dict(state, not_a_member=1.0)
    blank = type(obj).__new__(type(obj))
    with pytest.raises(RuntimeError, match="C\\+\\+ members"):
        blank.__setstate__(extra)


def test_a_changed_field_survives_the_round_trip() -> None:
    """Not just the factory defaults: a field the caller set is kept."""
    s = cb.florschuetz_1981_inline()
    s.standard_error = 0.123
    s.A_fit = cb.JetArrayGeometricFit(9.0, 8.0, 7.0, 6.0)
    s.valid_z_d = cb.ImpingementRange(0.5, 4.5)
    s.pattern = cb.JetHolePattern.Staggered
    t = pickle.loads(pickle.dumps(s))
    assert t.standard_error == 0.123
    assert (t.A_fit.C, t.A_fit.nx, t.A_fit.ny, t.A_fit.nz) == (9.0, 8.0, 7.0, 6.0)
    assert (t.valid_z_d.lo, t.valid_z_d.hi) == (0.5, 4.5)
    assert t.pattern == cb.JetHolePattern.Staggered
    g = _orifice_geometry()
    u = copy.deepcopy(g)
    assert u.tap == cb.OrificeTap.Corner
    assert math.isclose(u.beta, 0.4)


# --- the forwarding wrappers ----------------------------------------------------


@pytest.mark.parametrize("cls", [_DuctHtc, _ChamberHtc, _FilmedHtc, _FilmSideHtc])
def test_a_blank_wrapper_raises_attribute_error_not_recursion(cls: type) -> None:
    """copy and pickle build the instance without __init__ and look up hooks
    on it; with _r unset, forwarding recursed without end."""
    blank = cls.__new__(cls)
    with pytest.raises(AttributeError):
        _ = blank._r
    with pytest.raises(AttributeError):
        _ = blank.anything
    with pytest.raises(AttributeError):
        _ = blank.__deepcopy__
    copy.copy(blank)
    copy.deepcopy(blank)


class _Picklable:
    """A stand-in surface result that copies cleanly."""

    def __init__(self) -> None:
        self.h = 123.0
        self.extra = [1.0, 2.0]


def test_a_populated_wrapper_deep_copies_and_still_forwards() -> None:
    w = _DuctHtc.__new__(_DuctHtc)
    w._r = _Picklable()
    w.dh_dT = 4.0
    c = copy.deepcopy(w)
    assert c.h == 123.0
    assert c.extra == [1.0, 2.0]
    assert c.extra is not w.extra
    assert c.dh_dT == 4.0
