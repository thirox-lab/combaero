"""EffusionPlateElement's opt-in flat-plate gas-side baseline.

`gas_baseline="flat_plate"` swaps the UNBLOWN coefficient that
`gas_augmentation` multiplies: the chamber's own surface correlation (the
default, fully developed duct) for a turbulent boundary layer restarted
`gas_start_length` upstream of the panel -- CR-2648's St_o, the unblown
reference of the Stanford full-coverage data, averaged over the panel. The
chamber still supplies the stream (T_gas, U_gas, rho); only the coefficient
changes. The default stays the chamber and is bit-identical to before.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    EffusionPlateElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
GEOM = {
    "hole_diameter": 0.8e-3,
    "wall_thickness": 2e-3,
    "pitch": 6e-3,
    "panel_length": 0.1,
    "panel_width": 0.1,
    "angle_deg": 30.0,
}


def _pb(name: str, Pt: float, Tt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _solve(chamber: bool = True, **kw) -> tuple[NetworkSolver, dict]:
    g = FlowNetwork()
    g.add_node(_pb("cool", 1.03e5, 600.0))
    if chamber:
        g.add_node(_pb("gas_in", 1.003e5, 1500.0))
        g.add_node(_pb("exit", 1.0e5, 1500.0))
        g.add_node(MomentumChamberNode("ch", main_inlet="duct"))
        g.add_element(
            ChannelElement("duct", "gas_in", "ch", length=0.2, diameter=0.08, roughness=0.0)
        )
        g.add_element(
            ChannelElement("tail", "ch", "exit", length=0.2, diameter=0.08, roughness=0.0)
        )
    else:
        g.add_node(_pb("ch", 1.0e5, 1500.0))
    g.add_element(EffusionPlateElement("eff", "cool", "ch", **{**GEOM, **kw}))
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


def _diag(**kw) -> dict:
    return _solve(**kw)[1]["__element_diag__"]["eff"]


# ---------------------------------------------------------------------------
# The correlation (C++), against its source and an independent integral
# ---------------------------------------------------------------------------


def test_the_local_value_is_cr2648s_printed_st0() -> None:
    re, pr = 4.0e5, 0.71
    st0 = 0.0295 * re**-0.2 * pr**-0.4  # "ST*PR**0.4=0.0295*REX**(-.2)"
    nu, dnu = cb.flat_plate_turbulent_nusselt_and_derivative(re, pr)
    assert nu == pytest.approx(st0 * re * pr, rel=1e-13)
    h = 1e-6 * re
    fd = (
        cb.flat_plate_turbulent_nusselt_and_derivative(re + h, pr)[0]
        - cb.flat_plate_turbulent_nusselt_and_derivative(re - h, pr)[0]
    ) / (2 * h)
    assert dnu == pytest.approx(fd, rel=1e-8)


def test_the_panel_average_is_the_integral_of_the_local_value() -> None:
    re_l, pr, s = 5.0e4, 0.7, 2.0
    n = 100000
    integral = (
        sum(
            cb.flat_plate_turbulent_nusselt_and_derivative(re_l * (s + (i + 0.5) / n), pr)[0]
            / (s + (i + 0.5) / n)
            for i in range(n)
        )
        / n
    )
    fp = cb.flat_plate_turbulent_panel_nusselt(re_l, pr, s)
    assert fp.Nu == pytest.approx(integral, rel=1e-8)
    nu, dnu = cb.flat_plate_turbulent_panel_nusselt_and_derivative(re_l, pr, s)
    h = 1e-6 * re_l
    fd = (
        cb.flat_plate_turbulent_panel_nusselt(re_l + h, pr, s).Nu
        - cb.flat_plate_turbulent_panel_nusselt(re_l - h, pr, s).Nu
    ) / (2 * h)
    assert nu == fp.Nu and fp.dNu_dRe_L == dnu
    assert dnu == pytest.approx(fd, rel=1e-8)
    # Re_x 1e5 to 1.5e5: the leading edge is below the Stanford span.
    assert fp.Re_x_start == pytest.approx(1.0e5) and fp.Re_x_end == pytest.approx(1.5e5)
    assert fp.extrapolated


# ---------------------------------------------------------------------------
# The option
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kw", "match"),
    [
        ({"gas_baseline": "flat_plate"}, "gas_start_length >= 0"),
        ({"gas_baseline": "flat_plate", "gas_start_length": -0.01}, "gas_start_length >= 0"),
        ({"gas_start_length": 0.01}, "flat-plate baseline's input"),
        ({"gas_baseline": "duct"}, "'chamber' or 'flat_plate'"),
    ],
)
def test_the_option_is_validated(kw: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        EffusionPlateElement("eff", "cool", "ch", **{**GEOM, **kw})


def test_the_flat_plate_needs_a_panel_length() -> None:
    geom = {k: v for k, v in GEOM.items() if k not in ("panel_length", "panel_width")}
    with pytest.raises(ValueError, match="needs panel_length"):
        EffusionPlateElement(
            "eff",
            "cool",
            "ch",
            panel_area=0.01,
            gas_baseline="flat_plate",
            gas_start_length=0.0,
            **geom,
        )


def test_the_flat_plate_needs_a_flowing_gas_side() -> None:
    """On a plenum there is no gas velocity, hence no boundary layer: refused
    rather than silently ignored."""
    with pytest.raises(ValueError, match="needs a flowing gas side"):
        _solve(chamber=False, gas_baseline="flat_plate", gas_start_length=0.01)


# ---------------------------------------------------------------------------
# The default is unchanged
# ---------------------------------------------------------------------------


def test_the_default_is_the_chamber_and_adds_nothing() -> None:
    """Bit-identical to before: the explicit chamber baseline gives the same
    diagnostics, key for key and bit for bit, and no flat-plate key."""
    d_default = _diag()
    d_chamber = _diag(gas_baseline="chamber")
    assert d_default.keys() == d_chamber.keys()
    for k, v in d_default.items():
        w = d_chamber[k]
        assert (v == w) or (isinstance(v, float) and math.isnan(v) and math.isnan(w)), k
    assert not any(k.startswith("h_gas_unblown_x") for k in d_default)
    assert "h_gas_unblown_extrapolated" not in d_default
    assert d_default["h_gas_unblown_model"].startswith("SmoothModel")


# ---------------------------------------------------------------------------
# The flat-plate gas side
# ---------------------------------------------------------------------------


def test_the_flat_plate_coefficient_is_the_correlation_at_the_chambers_stream() -> None:
    x0 = 0.05
    s, r = _solve(gas_baseline="flat_plate", gas_start_length=x0)
    d = r["__element_diag__"]["eff"]
    dc = _diag()
    ch = s.network.nodes["ch"]
    # Recomputed here from the chamber's own static state.
    state = s._derived_states["ch"]
    from combaero.network import NetworkMixtureState

    st = NetworkMixtureState(P=r["ch.P"], Pt=r["ch.Pt"], Tt=state[0], m_dot=0.0, Y=Y_AIR)
    T_static = ch._closure(st, ch._total_m_dot).T_static
    ts = cb.transport_state(T_static, r["ch.P"], st.X)
    L = GEOM["panel_length"]
    re_l = cb.density(T_static, r["ch.P"], st.X) * d["U_gas"] * L / ts.mu
    expected = cb.flat_plate_turbulent_panel_nusselt(re_l, ts.Pr, x0 / L).Nu * ts.k / L
    assert d["h_gas_unblown"] == pytest.approx(expected, rel=1e-6)
    assert d["h_gas_unblown_model"] == "flat_plate_turbulent(CR-2648 St_o, panel average)"
    assert d["h_gas_unblown_Re"] == pytest.approx(re_l, rel=1e-6)
    assert d["h_gas_unblown_x_start"] == x0
    assert d["h_gas_unblown_x_end"] == pytest.approx(x0 + L)
    assert d["h_gas_unblown_Re_x_start"] == pytest.approx(re_l * x0 / L, rel=1e-12)
    # Only the coefficient changes: the stream is the chamber's.
    for k in ("T_gas", "U_gas", "blowing_ratio", "density_ratio", "T_gas_from_main_inlet"):
        assert d[k] == pytest.approx(dc[k], rel=1e-6), k
    assert d["h_gas_unblown"] != pytest.approx(dc["h_gas_unblown"], rel=0.05)
    assert d["eta_overall"] != pytest.approx(dc["eta_overall"], rel=1e-3)


def test_augmentation_multiplies_the_flat_plate() -> None:
    d1 = _diag(gas_baseline="flat_plate", gas_start_length=0.0)
    d2 = _diag(gas_baseline="flat_plate", gas_start_length=0.0, gas_augmentation=1.5)
    assert d2["h_gas_unblown"] == pytest.approx(d1["h_gas_unblown"], rel=1e-9)
    assert d2["h_gas"] == pytest.approx(1.5 * d1["h_gas_unblown"], rel=1e-9)
    assert d2["T_wall_hot"] > d1["T_wall_hot"]


def test_a_longer_run_thickens_the_boundary_layer() -> None:
    """h ~ x^-0.2: the farther the boundary layer has run, the lower the
    coefficient and the cooler the wall. The start length moves the result."""
    h = [
        _diag(gas_baseline="flat_plate", gas_start_length=x0)["h_gas_unblown"]
        for x0 in (0.0, 0.05, 0.5)
    ]
    assert h[0] > h[1] > h[2]
    # A fully turbulent Re_x from zero is outside the data span, and says so.
    assert (
        _diag(gas_baseline="flat_plate", gas_start_length=0.0)["h_gas_unblown_extrapolated"] == 1.0
    )


def test_no_gas_flow_keeps_the_flat_plate_keys_as_nan() -> None:
    """A chamber below the no-through-flow floor has no gas side: the keys
    stay, NaN, the flat plate's with them (#998's keys-stay rule)."""
    plate = EffusionPlateElement(
        "eff", "cool", "ch", gas_baseline="flat_plate", gas_start_length=0.01, **GEOM
    )
    s, _r = _solve(gas_baseline="flat_plate", gas_start_length=0.01)
    ch = s.network.nodes["ch"]
    plate._discharge_node = ch
    saved = ch._total_m_dot
    try:
        ch._total_m_dot = 0.0
        from combaero.network import NetworkMixtureState

        cool = NetworkMixtureState(P=1.03e5, Pt=1.03e5, Tt=600.0, m_dot=1e-3, Y=Y_AIR)
        gas = NetworkMixtureState(P=1.0e5, Pt=1.0e5, Tt=1500.0, m_dot=0.0, Y=Y_AIR)
        out = plate.wall_heat_transfer(cool, gas)
    finally:
        ch._total_m_dot = saved
    assert out["gas_side_chamber"] == 1.0
    for k in ("h_gas_unblown_x_start", "h_gas_unblown_extrapolated", "h_gas_unblown"):
        assert math.isnan(out[k]), k
