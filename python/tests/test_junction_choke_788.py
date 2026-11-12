"""A channel into a junction: branch selection and the choked-bounded exit (#788).

The compressible momentum-chamber closure Pt = Pt_cl(P, T0, |m|/A) falls to a
minimum Pt* at the sonic static pressure P* and rises again below it, so one
Pt had two static roots and Newton could land on the supersonic one. Below P*
the closure is now continued by its mirror image about Pt*, 2 Pt* - Pt_cl(P):
one root per Pt, the subsonic one whenever the chamber can carry its flux at
that Pt.

A Pt below Pt* then roots on the continuation: the chamber cannot carry its
flux at that stagnation pressure. Fed by a choked channel of its own area,
with nothing reading its static pressure, that is the channel's
under-expanded exit -- the "choked bounded" exit a pressure boundary takes
(#636) -- and the root is the one a PlenumNode in the junction's place gives.
Such a root is accepted as a choked-bounded junction; any other is refused.

The network is the issue's: 32 kg/s air at 897 K and 1 kg/s CH4 through an
orifice into a plenum, a combustor, a 1 m x 0.1 m channel, an exit orifice
and a 20 bar outlet.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero import _solver_tools
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    ChannelElement,
    CombustorNode,
    ConvectiveSurface,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    RibbedModel,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
_X_CH4 = [0.0] * cb.num_species()
_X_CH4[cb.species_index_from_name("CH4")] = 1.0
Y_CH4 = list(cb.mole_to_mass(_X_CH4))
JCT = "__junction__ch__ox"
P_OUT = 20.0e5


def _net(
    bore: float | None,
    *,
    correlation: str = "fixed",
    in_pipe: bool = True,
    junction: str = "momentum",
    exit_Pt: float = P_OUT,
    junction_area: float | None = None,
) -> FlowNetwork:
    """bore None: the channel discharges straight into the outlet at exit_Pt."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=32.0, Tt=897.0, Y=Y_AIR))
    n.add_node(MassFlowBoundary("fuel", m_dot=1.0, Tt=300.0, Y=Y_CH4))
    n.add_node(PlenumNode("pl"))
    n.add_node(CombustorNode("cmb"))
    n.add_node(PressureBoundary("out", Pt=exit_Pt, Tt=300.0, Y=Y_AIR))
    n.add_element(LosslessConnectionElement("ea", "air", "pl"))
    n.add_element(OrificeElement("of", "fuel", "pl", Cd=0.6, diameter=0.02, correlation="fixed"))
    n.add_element(LosslessConnectionElement("ep", "pl", "cmb"))
    if bore is None:
        n.add_element(ChannelElement("ch", "cmb", "out", length=1.0, diameter=0.1))
        return n
    if junction == "plenum":
        n.add_node(PlenumNode(JCT))
    else:
        n.add_node(MomentumChamberNode(JCT, area=junction_area))
    n.add_element(ChannelElement("ch", "cmb", JCT, length=1.0, diameter=0.1))
    kw: dict = {"diameter": bore, "correlation": correlation, "in_pipe": in_pipe}
    if correlation == "fixed":
        kw["Cd"] = 0.6
    n.add_element(OrificeElement("ox", JCT, "out", **kw))
    return n


def _solve(net: FlowNetwork, **kw: object) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    # As the GUI and the stress harness solve: with the Newton-homotopy retry.
    s.homotopy_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30.0, **kw)
    return s, r


# --- the closure: one root per Pt --------------------------------------------

_T0, _A, _M_DOT = 2000.0, 7.854e-3, 33.0


def _closure(P: float, Pt: float, m: float = _M_DOT, T0: float = _T0, Y=Y_AIR):
    return _solver_tools.momentum_chamber_residual_and_jacobian(P, Pt, m, T0, list(Y), _A)


@pytest.mark.parametrize("pt_over_ptstar", [0.4, 0.9, 1.1, 2.0])
def test_the_closure_has_one_root_per_total_pressure(pt_over_ptstar: float) -> None:
    """At fixed (T0, flux) the residual is monotone in P across the sonic
    point, so every Pt has exactly one root: subsonic for Pt >= Pt*, past the
    sonic point (bounded) below it. Before #788 a Pt above Pt* had two roots
    and one below none."""
    Pt_star = _closure(1.0e6, 5.0e6).Pt_sonic  # bounded at 10 bar: P* ~ 26 bar
    P_star = _closure(1.0e6, 5.0e6).P_sonic
    assert Pt_star > 0.0 and P_star > 0.0
    Pt = pt_over_ptstar * Pt_star
    Ps = np.geomspace(0.02 * P_star, 3.0 * Pt_star, 3000)
    r = np.array([_closure(P, Pt).residual for P in Ps])
    assert np.all(np.diff(r) < 0.0), "not monotone in P"
    roots = np.flatnonzero(np.sign(r[:-1]) != np.sign(r[1:]))
    assert roots.size == 1
    P_root = Ps[roots[0]]
    assert bool(P_root > P_star) is (pt_over_ptstar > 1.0)
    assert _closure(P_root, Pt).bounded is bool(P_root < P_star)


@pytest.mark.parametrize("P_over_Pstar", [0.95, 0.6, 0.2])
@pytest.mark.parametrize("sign", [1.0, -1.0])
def test_the_continuation_jacobian_matches_central_differences(
    P_over_Pstar: float, sign: float
) -> None:
    """Every partial of the continued row, on a CH4-burnt mix: P, Pt, m_dot
    (either sign), T0 and each species present."""
    Y = list(0.93 * np.array(Y_AIR) + 0.07 * np.array(Y_CH4))
    P_star = _closure(1.0e6, 5.0e6, Y=Y).P_sonic
    P, Pt, m = P_over_Pstar * P_star, 3.0e6, sign * _M_DOT
    r = _closure(P, Pt, m, Y=Y)
    assert r.bounded

    def cd(f, x: float, h: float) -> float:
        return (f(x + h) - f(x - h)) / (2.0 * h)

    checks = {
        "P": (r.d_res_dP, lambda v: _closure(v, Pt, m, Y=Y).residual, P, 1e-6 * P),
        "Pt": (r.d_res_dP_total, lambda v: _closure(P, v, m, Y=Y).residual, Pt, 1e-6 * Pt),
        "m": (r.d_res_dmdot, lambda v: _closure(P, Pt, v, Y=Y).residual, m, 1e-6 * abs(m)),
        "T0": (r.d_res_dT, lambda v: _closure(P, Pt, m, T0=v, Y=Y).residual, _T0, 1e-6 * _T0),
    }
    for k, (an, f, x, h) in checks.items():
        assert an == pytest.approx(cd(f, x, h), rel=1e-6), k
    for k, y in enumerate(Y):
        if y < 1e-3:
            continue

        def fy(v: float, k: int = k) -> float:
            Yk = list(Y)
            Yk[k] = v
            return _closure(P, Pt, m, Y=Yk).residual

        assert r.d_res_dY[k] == pytest.approx(cd(fy, y, 1e-6), rel=1e-5), f"Y[{k}]"


def test_the_subsonic_closure_is_unchanged() -> None:
    """Below the sonic point nothing moves: same residual and partials as the
    closure's own Pt - P exp(ln(Pt/P)), bounded False."""
    r = _closure(4.0e6, 5.0e6)
    assert not r.bounded and r.M < 1.0
    assert r.residual == pytest.approx(5.0e6 - 4.0e6 - r.q_dyn, rel=1e-12)
    assert r.P_sonic == 0.0 and r.Pt_sonic == 0.0


# --- part 1: Newton cannot land on the supersonic root -----------------------


@pytest.mark.parametrize(("P0", "Pt0"), [(1.7e6, 5.67e6), (2.0e6, 5.2e6), (1.0e6, 6.0e6)])
def test_newton_seeded_on_the_supersonic_root_finds_the_subsonic_one(P0: float, Pt0: float) -> None:
    """The issue's 95 mm root: seeded where it once "converged" (junction at
    Mach 1.46, P 17 bar, Pt 57 bar), the solve lands on the subsonic root,
    Mach 0.74 at 52 bar (fixed Cd 0.6). Before #788 these seeds ended on the
    supersonic root (Mach 1.15-1.5), refused INCONSISTENT or stalled."""
    net = _net(0.095)
    net.nodes[JCT].initial_guess = {f"{JCT}.P": P0, f"{JCT}.Pt": Pt0}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(net).solve(auto_retry=False)
    assert r["__success__"], r["__message__"]
    assert r["__supersonic_chambers__"] == [] and r["__choked_bounded_junctions__"] == []
    assert r[f"{JCT}.mach"] == pytest.approx(0.739, abs=0.005)
    assert r[f"{JCT}.Pt"] == pytest.approx(5.221e6, rel=2e-3)
    assert r["ch.choked"] == 0.0
    assert r[f"{JCT}.P"] > r["ox.P_out"]


# --- part 2: the choked-bounded junction -------------------------------------

_BOUNDED = [
    pytest.param(0.097, "ReaderHarrisGallagher", True, id="RHG-97mm"),
    pytest.param(0.099, "fixed", True, id="fixed-99mm"),
    pytest.param(0.15, "fixed", False, id="lumped-150mm"),
    pytest.param(0.1, "fixed", True, id="bore-equals-pipe"),
]


@pytest.fixture(scope="module")
def choked_reference() -> dict:
    """The channel straight into the outlet: choked, the outlet taking it
    choked bounded (#636)."""
    _, r = _solve(_net(None))
    assert r["__success__"] and r["ch.choked"] == 1.0
    return r


@pytest.mark.parametrize(("bore", "correlation", "in_pipe"), _BOUNDED)
def test_a_choked_channel_into_a_junction_is_choked_bounded(
    bore: float, correlation: str, in_pipe: bool, choked_reference: dict
) -> None:
    """Before #788 these had no root (no_progress). Now: converged and
    accepted, the junction reported as the plenum its root is; the channel
    choked as into the outlet directly, so the combustor's pressure is that
    network's; and the junction's Pt is the PlenumNode network's."""
    s, r = _solve(_net(bore, correlation=correlation, in_pipe=in_pipe))
    assert r["__success__"], r["__message__"]
    assert r["__choked_bounded_junctions__"] == [JCT]
    assert r["__supersonic_chambers__"] == []
    assert r["ch.m_dot"] == pytest.approx(33.0, rel=1e-9)
    assert r["ch.choked"] == 1.0 and r["ch.mach_out"] == pytest.approx(1.0, abs=1e-6)
    assert r["cmb.Pt"] == pytest.approx(choked_reference["cmb.Pt"], rel=1e-7)
    # The plenum in its place: the same root, to the solves' tolerance.
    _, rp = _solve(_net(bore, correlation=correlation, in_pipe=in_pipe, junction="plenum"))
    assert rp["__success__"]
    assert r[f"{JCT}.Pt"] == pytest.approx(rp[f"{JCT}.Pt"], rel=1e-6)
    assert r["ox.m_dot"] == pytest.approx(rp["ox.m_dot"], rel=1e-9)
    # Reported as a plenum: P = Pt, at rest; the sonic state of the flux and
    # the expansion loss beside it.
    Pt_j = r[f"{JCT}.Pt"]
    assert r[f"{JCT}.P"] == Pt_j and r[f"{JCT}.mach"] == 0.0
    assert r[f"{JCT}.T"] == pytest.approx(r[f"{JCT}.Tt"], rel=1e-12)
    assert r[f"{JCT}.choked_bounded"] == 1.0
    assert r[f"{JCT}.Pt_sonic"] > Pt_j
    assert r[f"{JCT}.expansion_loss"] == pytest.approx(r[f"{JCT}.Pt_sonic"] - Pt_j, rel=1e-12)
    # Static pressure falls along the flow: channel exit -> junction -> outlet.
    assert r["ch.P_out"] == Pt_j
    assert Pt_j >= P_OUT * (1.0 - 1e-9)
    # Mass and energy close; no element creates flow work.
    assert r["__second_law_failures__"] == []
    assert r["__second_law__"]["ch"]["verdict"] == "pass"
    assert r["__second_law__"]["ch"]["sigma_norm"] > 0.0
    assert abs(r["__ledger__"]["mass_closure_rel"]) < 1e-8
    assert abs(r["__ledger__"]["energy_closure_rel"]) < 1e-8
    # The unknown itself sits on the continuation, which no other row reads.
    x = np.array(r["__x_solution__"])
    assert x[s.unknown_names.index(f"{JCT}.P")] < r[f"{JCT}.P_sonic"]


def test_the_network_jacobian_at_a_choked_bounded_root_matches_central_differences() -> None:
    """FD of every row at the root and off it, the junction's continued row
    included. On the lumped orifice: an orifice IN the pipe reads its approach
    from the junction's Pt at the pipe's flux, which past the pipe's own
    sonic limit carries a pre-existing Jacobian gap in its flow columns (the
    same with a PlenumNode junction, on main too; not this row)."""
    s, r = _solve(_net(0.15, in_pipe=False))
    assert r["__choked_bounded_junctions__"] == [JCT]
    x_root = np.array(r["__x_solution__"])
    rng = np.random.default_rng(788)
    i_P = s.unknown_names.index(f"{JCT}.P")
    for x in (x_root, x_root * (1.0 + 0.003 * rng.standard_normal(x_root.size))):
        F, J = s._residuals_and_jacobian(x)
        J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
        assert x[i_P] == s._get_node_state(s.network.nodes[JCT], x).P
        assert s.network.nodes[JCT]._closure(s._get_node_state(s.network.nodes[JCT], x)).bounded
        Jfd = fd_jacobian(s, x)
        for i in range(J.shape[0]):
            scale = max(np.max(np.abs(J[i])), 1e-12)
            assert np.max(np.abs(J[i] - Jfd[i])) / scale < 1e-5, s.unknown_names[i]


# --- refused: past the sonic limit, but not a choked duct's exit -------------


def _lossless_fed_chamber() -> FlowNetwork:
    """A chamber too small for the flow a wide duct draws through it, fed
    by a lossless connection: nothing chokes upstream of it. (A wide ORIFICE
    there, the fixture before #1016, is now refused at resolve: its bore
    was 3.5x the chamber's area.)"""
    n = FlowNetwork()
    n.add_node(PressureBoundary("in", Pt=10.0e5, Tt=300.0, Y=Y_AIR))
    n.add_node(MomentumChamberNode("mc", area=2.0e-4))
    n.add_node(PressureBoundary("out", Pt=1.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(LosslessConnectionElement("ll", "in", "mc"))
    n.add_element(ChannelElement("o", "mc", "out", length=0.05, diameter=0.03))
    return n


def test_a_bounded_chamber_fed_by_a_lossless_connection_is_refused() -> None:
    _, r = _solve(_lossless_fed_chamber())
    assert not r["__success__"]
    assert r["__supersonic_chambers__"] == ["mc"]
    assert r["__choked_bounded_junctions__"] == []
    assert "'ll' (LosslessConnectionElement) feeds it, not a choked channel" in r["__message__"]


def test_a_junction_narrower_than_its_choked_channel_is_refused() -> None:
    """The chamber is the throat: it would choke before the channel. The
    outlet is no restriction: a bore equal to the junction's area is no
    orifice (#1016; the fixture's 0.15 m bore, 3.75x the junction's area,
    is now refused at resolve)."""
    A_jct = 0.6 * np.pi * 0.05**2
    _, r = _solve(_net(np.sqrt(4.0 * A_jct / np.pi), in_pipe=False, junction_area=A_jct))
    assert not r["__success__"]
    assert r["__choked_bounded_junctions__"] == []
    assert JCT in r["__supersonic_chambers__"]


def test_a_junction_whose_static_pressure_is_read_is_refused() -> None:
    """A ribbed channel downstream reads the junction's static P (its density):
    past the sonic limit that static has no meaning, so the root is refused."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=40.0, Tt=1500.0, Y=Y_AIR))
    n.add_node(PlenumNode("pl"))
    n.add_node(MomentumChamberNode(JCT))
    n.add_node(PressureBoundary("out", Pt=5.0e5, Tt=300.0, Y=Y_AIR))
    n.add_element(LosslessConnectionElement("ea", "air", "pl"))
    n.add_element(ChannelElement("ch", "pl", JCT, length=1.0, diameter=0.1))
    rib = ChannelElement(
        "rb",
        JCT,
        "out",
        length=0.05,
        diameter=0.15,
        surface=ConvectiveSurface(area=0.06, model=RibbedModel(e_D=0.05, p_e=10.0)),
    )
    n.add_element(rib)
    _, r = _solve(n)
    assert r["ch.choked"] == 1.0
    assert not r["__success__"]
    assert r["__choked_bounded_junctions__"] == []
    assert r["__supersonic_chambers__"] == [JCT]
    assert "['rb.m_dot'] read its static pressure" in r["__message__"]


@pytest.mark.parametrize("correlation", ["fixed", "ReaderHarrisGallagher"])
def test_a_bore_equal_to_the_pipe_converges_cold(correlation: str) -> None:
    """No orifice at the junction (bore == pipe passes the flow unchanged), so
    the outlet's 20 bar fixes the junction's Pt far below the sonic Pt* of the
    pipe's flux. The start's P = Pt is then past the sonic point, and the
    first Newton step crossed the closure's zero-slope fold to a negative P:
    no_progress without retries. The start is seeded on the closure's root
    (_seed_chambers_past_sonic), and the plain cold solve converges."""
    s = NetworkSolver(_net(0.1, correlation=correlation))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(auto_retry=False)
    assert r["__success__"], r["__message__"]
    assert r["__choked_bounded_junctions__"] == [JCT]
    assert r[f"{JCT}.Pt"] == pytest.approx(P_OUT, rel=1e-9)


def test_a_start_on_the_subsonic_side_is_not_moved() -> None:
    """The seed moves only a chamber that starts past its sonic point; the
    bore = pipe start is moved onto its root (its Pt from the outlet)."""
    for bore, moved in ((0.08, False), (0.1, True)):
        net = _net(bore)
        net.resolve_all_topology()
        s = NetworkSolver(net)
        x0 = s._build_x0()
        s._propagate_states(x0, with_relay=False)
        st = s._get_node_state(net.nodes[JCT], x0)
        assert net.nodes[JCT].area > 0.0
        c = net.nodes[JCT]._closure(st)
        if moved:
            assert c.bounded and abs(c.residual) < 1e-6 * st.Pt
        else:
            assert st.Pt == st.P and not c.bounded
