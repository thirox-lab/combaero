"""Regression and Jacobian tests for #1128: PressureLossElement, AreaChangeElement,
and LosslessConnectionElement honour PressureBoundary.exit_head_lost at outflow
boundaries, with analytic Jacobians (C1 across sonic choke).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    AreaChangeElement,
    ChannelElement,
    CombustorNode,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
)
from combaero.network.components import FLOW_REVERSAL_M_SCALE
from combaero.network.pressure_loss import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    LinearThetaHeadLoss,
)

Y_AIR = list(cb.species.dry_air_mass())
X_AIR = list(cb.species.dry_air())
D = 0.05
A = math.pi * D * D / 4.0
KEY = "exit_head_recovered"


def _fuel_mix_mass() -> list[float]:
    X = [0.0] * cb.num_species()
    for name, x in (("N2", 0.72), ("O2", 0.18), ("CH4", 0.07), ("CO2", 0.03)):
        X[cb.species_index_from_name(name)] = x
    return list(cb.mole_to_mass(X))


def _solve(net: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(net).solve(timeout=30)
    assert r["__success__"], r.get("__message__")
    return r


def _s4(kind: str, coupling: str = "auto") -> dict:
    """Probe S4 from #1128: PB 1.05e5 -> channel -> chamber(A) -> X -> PB 1.0e5."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("a", area=A))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR, coupling=coupling))
    g.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
    if kind == "short_channel":
        g.add_element(ChannelElement("x", "a", "out", length=1e-4, diameter=D))
    elif kind == "pressure_loss_zero":
        g.add_element(
            PressureLossElement("x", "a", "out", ConstantHeadLoss(zeta=0.0, area=A), area=A)
        )
    elif kind == "lossless_chamber":
        g.add_element(LosslessConnectionElement("x", "a", "out"))
    elif kind == "area_change_4":
        g.add_element(AreaChangeElement("x", "a", "out", F0=A, F1=4.0 * A))
    elif kind == "area_change_passthrough":
        g.add_element(AreaChangeElement("x", "a", "out", F0=A, F1=A))
    elif kind == "orifice_passthrough":
        g.add_element(OrificeElement("x", "a", "out", Cd=0.6, diameter=D, correlation="fixed"))
    else:
        raise ValueError(kind)
    return _solve(g)


@pytest.mark.parametrize("coupling", ["auto", "static"])
@pytest.mark.parametrize(
    "kind",
    [
        "pressure_loss_zero",
        "lossless_chamber",
        "area_change_passthrough",
        "orifice_passthrough",
    ],
)
def test_s4_zero_loss_elements_match_short_channel(kind: str, coupling: str) -> None:
    """Zero-loss PressureLoss, Lossless out of a chamber, equal-face AreaChange,
    and pass-through Orifice match a same-area short channel into an
    exit_head_lost boundary (#1128)."""
    ref = _s4("short_channel", coupling=coupling)
    r = _s4(kind, coupling=coupling)
    # L = 1e-4 m friction is ~5e-5 of the exit head; flows match to 1e-3.
    assert r["x.m_dot"] == pytest.approx(ref["x.m_dot"], rel=1e-3)
    if KEY in r["__element_diag__"]["x"]:
        assert r["__element_diag__"]["x"][KEY] == 0.0
    assert r["a.P"] == pytest.approx(1.0e5, rel=1e-6)


@pytest.mark.parametrize("coupling", ["auto", "static"])
def test_s4_area_change_loses_exit_head(coupling: str) -> None:
    """AreaChange (F0=A -> F1=4A) loses q(F1) at a static outflow boundary and
    reports exit_head_recovered == 0."""
    r_lost = _s4("area_change_4", coupling=coupling)
    r_tot = _s4("area_change_4", coupling="total")
    assert r_lost["__element_diag__"]["x"][KEY] == 0.0
    assert r_tot["__element_diag__"]["x"][KEY] > 300.0
    assert r_lost["x.m_dot"] < r_tot["x.m_dot"]
    # Losing ~394 Pa out of 5 kPa reduces flow by ~3-5%.
    assert 0.90 < r_lost["x.m_dot"] / r_tot["x.m_dot"] < 0.98


def test_s4_total_coupling_still_recovers_head() -> None:
    """coupling='total' keeps exact total-to-total recovery on all three elements."""
    for kind in (
        "pressure_loss_zero",
        "lossless_chamber",
        "area_change_4",
        "area_change_passthrough",
    ):
        r = _s4(kind, coupling="total")
        assert r["__element_diag__"]["x"][KEY] > 300.0


def test_lossless_with_explicit_area_from_plenum_loses_exit_head() -> None:
    """A LosslessConnectionElement with an explicit area loses exit head from a plenum."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c", "in", "p", length=0.5, diameter=D))
    lc = LosslessConnectionElement("lc", "p", "out")
    lc.area = A
    g.add_element(lc)
    r = _solve(g)
    ref = _s4("short_channel")
    assert r["lc.m_dot"] == pytest.approx(ref["x.m_dot"], rel=1e-3)
    assert r["__element_diag__"]["lc"][KEY] == 0.0


def test_pressure_loss_chokes_smoothly_into_static_boundary() -> None:
    """A zero-loss PressureLossElement of area A between 3 bar and 1 bar chokes
    at the exact sonic mass flow of area A, independent of back pressure."""
    flows = []
    for p_back in (1.0e5, 0.7e5):
        g = FlowNetwork()
        g.add_node(PressureBoundary("hi", Pt=3.0e5, Tt=300.0, Y=Y_AIR))
        g.add_node(PressureBoundary("lo", Pt=p_back, Tt=300.0, Y=Y_AIR))
        g.add_element(
            PressureLossElement("pl", "hi", "lo", ConstantHeadLoss(zeta=0.0, area=A), area=A)
        )
        r = _solve(g)
        flows.append(r["pl.m_dot"])
        assert r["__element_diag__"]["pl"][KEY] == 0.0
    assert flows[0] == pytest.approx(flows[1], rel=1e-10)
    p_star = cb.critical_pressure_ratio(300.0, 3.0e5, X_AIR) * 3.0e5
    m_star = A * cb.mass_flux_isentropic(300.0, 3.0e5, p_star, X_AIR)
    assert flows[0] == pytest.approx(m_star, rel=1e-9)


@pytest.mark.parametrize("incompressible", [False, True])
@pytest.mark.parametrize("G", [0.0, 80.0, 220.0, 350.0])
def test_discharge_total_pressure_jacobian_vs_fd(incompressible: bool, G: float) -> None:
    """Analytic partials of discharge_total_pressure vs central FD in subsonic,
    near-sonic, choked, and incompressible regimes, including multi-species Y."""
    P, Tt = 1.0e5, 850.0
    Y = _fuel_mix_mass()
    res = _solver_tools.discharge_total_pressure(P, Tt, G, Y, incompressible)

    def eval_pt(p: float = P, t: float = Tt, g_val: float = G, y_val: list[float] = Y) -> float:
        return float(_solver_tools.discharge_total_pressure(p, t, g_val, y_val, incompressible).Pt)

    hP = 1.0
    fd_P = (eval_pt(p=P + hP) - eval_pt(p=P - hP)) / (2.0 * hP)
    assert res.dPt_dP == pytest.approx(fd_P, rel=1e-6, abs=1e-9)

    hT = 1e-2
    fd_T = (eval_pt(t=Tt + hT) - eval_pt(t=Tt - hT)) / (2.0 * hT)
    assert res.dPt_dTt == pytest.approx(fd_T, rel=1e-6, abs=1e-9)

    hG = 1e-3
    fd_G = (eval_pt(g_val=G + hG) - eval_pt(g_val=G - hG)) / (2.0 * hG)
    assert res.dPt_dG == pytest.approx(fd_G, rel=1e-6, abs=1e-6)

    for k, yk in enumerate(Y):
        if yk <= 0.0:
            continue
        hY = 1e-6
        yp = list(Y)
        ym = list(Y)
        yp[k] += hY
        ym[k] -= hY
        fd_Y = (eval_pt(y_val=yp) - eval_pt(y_val=ym)) / (2.0 * hY)
        assert res.dPt_dY[k] == pytest.approx(fd_Y, rel=1e-5, abs=1e-6 * res.Pt)


def test_discharge_total_pressure_is_c1_across_sonic_choke() -> None:
    """Subsonic and choked branches of discharge_total_pressure meet C1 at P = P*."""
    Tt, G = 600.0, 200.0
    Y = _fuel_mix_mass()
    P_star = float(_solver_tools.sonic_static_pressure(Tt, G, Y).P)
    eps = 1e-6 * P_star
    sub = _solver_tools.discharge_total_pressure(P_star + eps, Tt, G, Y, False)
    chk = _solver_tools.discharge_total_pressure(P_star - eps, Tt, G, Y, False)
    assert not sub.choked and chk.choked
    assert sub.Pt == pytest.approx(chk.Pt, rel=1e-10)
    assert abs(sub.dPt_dP) < 1e-5 and chk.dPt_dP == 0.0
    assert sub.dPt_dG == pytest.approx(chk.dPt_dG, rel=1e-5)
    assert sub.dPt_dTt == pytest.approx(chk.dPt_dTt, rel=1e-5)


@pytest.mark.parametrize("model", ["sharp", "conical"])
@pytest.mark.parametrize("incompressible", [False, True])
@pytest.mark.parametrize("drive", [-50.0, -0.5, 0.0, 0.5, 50.0, 2.0e4, 2.5e5])
def test_area_change_element_flow_exit_head_lost_jacobian_vs_fd(
    model: str, incompressible: bool, drive: float
) -> None:
    """Analytic partials of area_change_element_flow with exit_head_lost_b=True
    (and exit_head_lost_a=False) vs central FD across zero drive, sub-floor,
    subsonic, and choked states."""
    if incompressible and abs(drive) > 1.0e5:
        pytest.skip("incompressible test at moderate drive")
    Y_a = _fuel_mix_mass()
    Y_b = Y_AIR
    Pt_a = 3.0e5
    Pt_b = Pt_a - drive
    args = {
        "Pt_a": Pt_a,
        "Tt_a": 420.0,
        "Y_a": Y_a,
        "Pt_b": Pt_b,
        "Tt_b": 600.0,
        "Y_b": Y_b,
    }
    geo = {
        "F_a": 1.0e-3,
        "F_b": 4.0e-3,
        "model": model,
        "length": 0.08,
        "D_h": 0.0,
        "incompressible": incompressible,
        "exit_head_lost_b": True,
        "exit_head_lost_a": False,
    }
    f = cb.area_change_element_flow(**args, **geo, with_derivatives=True)

    def m(**kw: object) -> float:
        return float(
            cb.area_change_element_flow(**{**args, **kw}, **geo, with_derivatives=False).m_dot
        )

    scale = abs(f.m_dot) + abs(f.dm_dPt_a) * 1.0 + 1e-12
    h_p = 1e-4 if abs(drive) < 1.0 else 1e-2
    for key, attr, h in (
        ("Pt_a", "dm_dPt_a", h_p),
        ("Pt_b", "dm_dPt_b", h_p),
        ("Tt_a", "dm_dTt_a", 1e-2),
        ("Tt_b", "dm_dTt_b", 1e-2),
    ):
        fd = (m(**{key: args[key] + h}) - m(**{key: args[key] - h})) / (2.0 * h)
        assert getattr(f, attr) == pytest.approx(fd, rel=1e-4, abs=1e-8 * scale), (
            model,
            incompressible,
            drive,
            key,
        )
    for node, grad in (("Y_a", f.dm_dY_a), ("Y_b", f.dm_dY_b)):
        Y = args[node]
        for k, yk in enumerate(Y):
            if yk <= 0.0:
                continue
            h = 1e-6
            yp = list(Y)
            ym = list(Y)
            yp[k] += h
            ym[k] -= h
            fd = (m(**{node: yp}) - m(**{node: ym})) / (2.0 * h)
            assert grad[k] == pytest.approx(fd, rel=1e-4, abs=1e-6 * scale), (node, k)


def test_area_change_low_re_exit_face_choke_jacobian_vs_fd() -> None:
    """At very low Re where r(Re) > 1 causes the exit face to choke (Me = 1)
    before M1 = 1, area_change_flow still chokes smoothly with exact partials."""
    # Tiny orifice-scale expansion at low pressure: Re < 50 at high drive.
    Pt0, Tt0 = 1.0e4, 1000.0
    F_in, F_out = 1.0e-8, 4.0e-8
    f = cb.area_change_flow(
        Pt0, Tt0, Y_AIR, 1.0e2, F_in, F_out, "sharp", 0.0, 0.0, False, True, True
    )
    assert f.choked
    assert f.mach < 1.0  # exit face choked before inlet reached M = 1
    assert f.dm_dPt_target == 0.0

    def m(pt: float = Pt0, tt: float = Tt0) -> float:
        return float(
            cb.area_change_flow(
                pt, tt, Y_AIR, 1.0e2, F_in, F_out, "sharp", 0.0, 0.0, False, False, True
            ).m_dot
        )

    fd_pt = (m(pt=Pt0 + 1.0) - m(pt=Pt0 - 1.0)) / 2.0
    fd_tt = (m(tt=Tt0 + 0.1) - m(tt=Tt0 - 0.1)) / 0.2
    assert f.dm_dPt0 == pytest.approx(fd_pt, rel=1e-5)
    assert f.dm_dTt0 == pytest.approx(fd_tt, rel=1e-5)


def _check_network_jacobian(g: FlowNetwork, x: np.ndarray, tol: float = 1e-5) -> None:
    s = NetworkSolver(g)
    s._residuals(x)
    s._build_x0()
    _, jac = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    scale = np.abs(fd * x).max(axis=1, keepdims=True)
    scale[scale == 0.0] = 1.0
    err = float((np.abs(jac.toarray() - fd) * np.abs(x) / scale).max())
    assert err < tol, f"max scaled Jacobian error = {err:.3e}"


@pytest.mark.parametrize(
    "kind",
    [
        "pressure_loss_zero",
        "lossless_chamber",
        "area_change_4",
        "area_change_passthrough",
    ],
)
def test_s4_global_network_jacobian_vs_fd(kind: str) -> None:
    """Full network Jacobian vs central FD on all S4 elements with exit_head_lost."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("a", area=A))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
    if kind == "pressure_loss_zero":
        g.add_element(
            PressureLossElement("x", "a", "out", ConstantHeadLoss(zeta=0.0, area=A), area=A)
        )
    elif kind == "lossless_chamber":
        g.add_element(LosslessConnectionElement("x", "a", "out"))
    elif kind == "area_change_4":
        g.add_element(AreaChangeElement("x", "a", "out", F0=A, F1=4.0 * A))
    else:
        g.add_element(AreaChangeElement("x", "a", "out", F0=A, F1=A))
    r = _solve(g)
    _check_network_jacobian(g, np.array(r["__x_solution__"], dtype=float))


def test_combustor_theta_pressure_loss_exit_head_lost_jacobian_vs_fd() -> None:
    """LinearThetaHeadLoss after a CombustorNode into a static outflow boundary:
    both burned-T theta term and discharge_total_pressure Tt/Y terms match FD."""
    Y_CH4 = [0.0] * cb.num_species()
    Y_CH4[cb.species_index_from_name("CH4")] = 1.0
    g = FlowNetwork()
    g.add_node(PressureBoundary("air", Pt=2.0e5, Tt=500.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("fuel", m_dot=0.015, Tt=300.0, Y=Y_CH4))
    g.add_node(CombustorNode("cmb", method="complete", area=0.01))
    g.add_node(PressureBoundary("out", Pt=1.5e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("inlet", "air", "cmb", length=0.2, diameter=0.06))
    g.add_element(LosslessConnectionElement("fline", "fuel", "cmb"))
    g.add_element(
        PressureLossElement(
            "exit",
            "cmb",
            "out",
            correlation=LinearThetaHeadLoss(k=0.5, zeta0=1.5, area=0.005),
            area=0.005,
        )
    )
    r = _solve(g)
    assert r["__element_diag__"]["exit"][KEY] == 0.0
    _check_network_jacobian(g, np.array(r["__x_solution__"], dtype=float), tol=1e-4)


@pytest.mark.parametrize("kind", ["passthrough", "lossless_area"])
@pytest.mark.parametrize("m_over_scale", [-20.0, -1.0, 0.5, 20.0])
def test_passthrough_and_lossless_reversal_jacobian_vs_fd(kind: str, m_over_scale: float) -> None:
    """Pass-through AreaChange and area-bearing Lossless into static boundaries
    have exact Jacobians across the zero-flow reversal blend."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=1.05e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PlenumNode("p"))
    g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    g.add_element(ChannelElement("ch", "in", "p", length=0.5, diameter=D))
    if kind == "passthrough":
        g.add_element(AreaChangeElement("e", "p", "out", F0=A, F1=A))
    else:
        lc = LosslessConnectionElement("e", "p", "out")
        lc.area = A
        g.add_element(lc)
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30)
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    x[s.unknown_names.index("e.m_dot")] = m_over_scale * FLOW_REVERSAL_M_SCALE
    _check_network_jacobian(g, x)


def test_fraction_loss_with_area_honours_exit_head_lost_and_reverses_symmetrically() -> None:
    """A ConstantFractionLoss with m_ref and an explicit area honours exit_head_lost
    symmetrically in forward and reversed orientation."""

    def build(mirrored: bool) -> FlowNetwork:
        g = FlowNetwork()
        g.add_node(PressureBoundary("a", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
        g.add_node(PressureBoundary("b", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
        g.add_node(PlenumNode("p"))
        e1, e2 = (("p", "a"), ("b", "p")) if mirrored else (("a", "p"), ("p", "b"))
        g.add_element(
            PressureLossElement(
                "loss", *e1, correlation=ConstantFractionLoss(xi=0.2, m_ref=0.05), area=1e-3
            )
        )
        g.add_element(
            PressureLossElement(
                "part", *e2, correlation=ConstantHeadLoss(zeta=2.0, area=1e-3), area=1e-3
            )
        )
        return g

    against = _solve(build(False))
    along = _solve(build(True))
    assert against["loss.m_dot"] == pytest.approx(-along["loss.m_dot"], rel=1e-8)
    assert against["p.Pt"] == pytest.approx(along["p.Pt"], rel=1e-8)
    assert against["__element_diag__"]["loss"][KEY] == 0.0
    assert along["__element_diag__"]["loss"][KEY] == 0.0
    # Because exit dynamic head is lost at boundary 'a', p.Pt > 1e5 / 0.8.
    assert along["p.Pt"] > 1.0e5 / 0.8 + 100.0


def test_passthrough_orifice_exit_head_lost_jacobian_and_choke() -> None:
    """Pass-through OrificeElement (bore == D_ch) honours exit_head_lost with
    exact analytic Jacobian (subsonic and choked) matching LosslessConnectionElement."""
    for p_in in (1.05e5, 3.0e5):
        g = FlowNetwork()
        g.add_node(PressureBoundary("in", Pt=p_in, Tt=300.0, Y=Y_AIR))
        g.add_node(MomentumChamberNode("a", area=A))
        g.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
        g.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
        g.add_element(OrificeElement("x", "a", "out", Cd=0.6, diameter=D, correlation="fixed"))
        r = _solve(g)
        x = np.array(r["__x_solution__"], dtype=float)
        if p_in == 1.05e5:
            _check_network_jacobian(g, x)
        else:
            # At p_in = 3 bar the upstream same-area channel sits on its own
            # choke cusp (M_out = 1); verify the pass-through orifice row's
            # choked Jacobian directly against central FD.
            s = NetworkSolver(g)
            s._build_x0()
            s._residuals(x)
            _, jac = s._residuals_and_jacobian(x)
            fd = fd_jacobian(s, x)
            row = s.unknown_names.index("x.m_dot")
            j_row = jac.toarray()[row]
            fd_row = fd[row]
            scale = float(np.abs(fd_row * x).max()) or 1.0
            err = float((np.abs(j_row - fd_row) * np.abs(x) / scale).max())
            assert err < 1e-6, f"choked pass-through orifice row err = {err:.3e}"
        # Compare against LosslessConnectionElement in the same position.
        g_ref = FlowNetwork()
        g_ref.add_node(PressureBoundary("in", Pt=p_in, Tt=300.0, Y=Y_AIR))
        g_ref.add_node(MomentumChamberNode("a", area=A))
        g_ref.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
        g_ref.add_element(ChannelElement("c", "in", "a", length=0.5, diameter=D))
        g_ref.add_element(LosslessConnectionElement("x", "a", "out"))
        ref = _solve(g_ref)
        assert r["x.m_dot"] == pytest.approx(ref["x.m_dot"], rel=1e-12)
        assert r["a.Pt"] == pytest.approx(ref["a.Pt"], rel=1e-12)


@pytest.mark.parametrize("pt_seed", [95000.0, 101678.0, 102000.0, 102322.0, 105000.0])
@pytest.mark.parametrize("coupling", ["auto", "total"])
def test_blocked_fraction_loss_converges_from_all_seeds(pt_seed: float, coupling: str) -> None:
    """Regression for the #1004/#1128 blocked-fraction NO_PROGRESS stall:
    ConstantFractionLoss without m_ref has dF/dm = 0 outside the narrow
    FLOW_REVERSAL_M_SCALE blend, so hybr can exit with NO_PROGRESS from seeds
    such as p.Pt = 102000, 102322, 105000 Pa (and on main at zeta=1, Pa=1.1e5).
    LM fallback rescues every seed to the exact blend root."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("a", Pt=1.1e5, Tt=300.0, Y=Y_AIR, coupling=coupling))
    g.add_node(PressureBoundary("b", Pt=1.0e5, Tt=300.0, Y=Y_AIR, coupling=coupling))
    g.add_node(PlenumNode("p"))
    g.add_element(
        PressureLossElement("main", "a", "b", ConstantHeadLoss(zeta=2.0, area=1e-3), area=1e-3)
    )
    g.add_element(PressureLossElement("loss", "a", "p", ConstantFractionLoss(xi=0.3), area=1e-3))
    g.add_element(
        PressureLossElement("part", "p", "b", ConstantHeadLoss(zeta=2.0, area=1e-3), area=1e-3)
    )
    s = NetworkSolver(g)
    orig_x0 = s._build_x0

    def _seeded_x0() -> np.ndarray:
        x0 = orig_x0()
        x0[s.unknown_names.index("p.P")] = pt_seed
        x0[s.unknown_names.index("p.Pt")] = pt_seed
        return x0

    s._build_x0 = _seeded_x0  # type: ignore[method-assign]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30)
    assert r["__success__"], r.get("__message__")
    assert abs(r["loss.m_dot"]) < FLOW_REVERSAL_M_SCALE
    assert r["__fraction_unanchored__"] == ["loss"]
