"""#903: the flow splitter -- one inlet to N outlets in given fractions.

Each outlet is a mass-flow controller: its flow is its fraction of the
splitter's throughput (a linear row), its total-pressure drop is free and
reported as the valve's metering, and one reference outlet is open
(Pt_out = Pt_in). These tests pin the fractions to solver precision, the
reported drops and metering against a hand calculation, the exact Jacobian,
the inlet state carried to every outlet, the refusals, and the C3X coolant
supply through one splitter against its per-hole flow boundaries. The
collector side (a many-inlet plenum or merge chamber) is checked against a
hand energy and species balance at very different inlet temperatures.
"""

from __future__ import annotations

import importlib.util
import math
import sys
import warnings
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from scipy.optimize import brentq

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    FlowSplitter,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
    SplitterOutletElement,
)
from combaero.network.components import SPLITTER_FRACTION_SUM_TOL
from combaero.network.pressure_loss import ConstantHeadLoss

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
N_SPECIES = len(cb.species.names)
P_SINK = 1.0e5
T_SUPPLY = 500.0
# Outlet node -> (fraction, head-loss zeta, loss area [m^2]) downstream of it.
PATHS = {
    "o1": (0.2, 6.0, 1.0e-4),
    "o2": (0.3, 1.5, 1.6e-4),
    "o3": (None, 1.0, 3.0e-4),
}
#: The path needing the most pressure at these fractions: the reference.
REFERENCE = "o1"
VALVE_AREA = 2.0e-4
M_SUPPLY = 0.05  # [kg/s] a mass-flow supply


def _solve(net: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return s, s.solve()


def _supply(net: FlowNetwork, kind: str, m_dot: float = M_SUPPLY) -> str:
    """The splitter's inlet node: a mass-flow boundary supplying it, or a
    plenum fed through an orifice from a pressure boundary."""
    if kind == "mass":
        net.add_node(MassFlowBoundary("sup", m_dot=m_dot, Tt=T_SUPPLY, Y=AIR))
        return "sup"
    net.add_node(PressureBoundary("res", Pt=3.0e5, Tt=T_SUPPLY, Y=AIR))
    net.add_node(PlenumNode("sup"))
    net.add_element(
        OrificeElement("feed", "res", "sup", diameter=0.03, Cd=0.8, correlation="fixed")
    )
    return "sup"


def _network(
    kind: str = "mass",
    outlets: dict[str, float | None] | None = None,
    reference: str | None = REFERENCE,
    m_supply: float = M_SUPPLY,
    auto_reference: bool = True,
    sink_coupling: str = "auto",
) -> tuple[FlowNetwork, FlowSplitter]:
    """supply -> splitter -> o_i -> head loss -> sink_i (one per outlet)."""
    net = FlowNetwork()
    inlet = _supply(net, kind, m_supply)
    for o, (_, zeta, area) in PATHS.items():
        net.add_node(PlenumNode(o))
        net.add_node(PressureBoundary(f"{o}_x", Pt=P_SINK, Tt=300.0, Y=AIR, coupling=sink_coupling))
        net.add_element(
            PressureLossElement(
                f"{o}_loss", o, f"{o}_x", correlation=ConstantHeadLoss(zeta, area), area=area
            )
        )
    sp = FlowSplitter(
        "sp",
        inlet,
        outlets if outlets is not None else {o: f for o, (f, _, _) in PATHS.items()},
        reference=reference,
        areas=dict.fromkeys(PATHS, VALVE_AREA),
        auto_reference=auto_reference,
    )
    sp.add_to(net)
    return net, sp


# --- the fractions -----------------------------------------------------------


@pytest.mark.parametrize("kind", ["mass", "pressure"])
@pytest.mark.parametrize(
    "outlets",
    [
        {"o1": 0.2, "o2": 0.3, "o3": None},  # a remainder
        {"o1": 0.2, "o2": 0.3, "o3": 0.5},  # all given, summing to 1
        {"o1": 0.1, "o2": 0.2, "o3": 0.7},  # 1 + 2.2e-16 in floating point
    ],
)
def test_the_fractions_hold_to_solver_precision(kind: str, outlets: dict) -> None:
    _, r = _solve(_network(kind, outlets)[0])
    assert r["__success__"], r["__message__"]
    m = {o: r[f"sp__{o}.m_dot"] for o in PATHS}
    M = sum(m.values())
    f = {
        o: (1.0 - sum(v for v in outlets.values() if v is not None) if v is None else v)
        for o, v in outlets.items()
    }
    for o in PATHS:
        assert m[o] / M == pytest.approx(f[o], abs=1e-13), o
    if kind == "mass":
        assert pytest.approx(M_SUPPLY, rel=1e-13) == M
    # The network's control volume closes and every outlet passes the second law.
    assert abs(r["__ledger__"]["mass_closure_rel"]) < 1e-12
    assert abs(r["__ledger__"]["energy_closure_rel"]) < 1e-12
    for o in PATHS:
        assert r["__second_law__"][f"sp__{o}"]["verdict"] == "pass"


def test_each_outlet_carries_the_inlet_state() -> None:
    """Tt and composition of the supply at every outlet; Pt per its path."""
    _, r = _solve(_network("mass")[0])
    assert r["__success__"], r["__message__"]
    for o in PATHS:
        assert r[f"{o}.Tt"] == pytest.approx(T_SUPPLY, rel=1e-14)
        for k in range(N_SPECIES):
            assert r[f"{o}.Y[{k}]"] == pytest.approx(AIR[k], abs=1e-15)
    pts = {r[f"{o}.Pt"] for o in PATHS}
    assert len(pts) == len(PATHS)  # each outlet at its own path's pressure


# --- the valve drops against a hand calculation ------------------------------


def _hand_Pt(m: float, zeta: float, area: float) -> float:
    """Pt at an outlet plenum whose head loss zeta m^2 / (2 rho A^2), rho at
    the plenum (Pt, T_SUPPLY), takes it to P_SINK: the root of
    Pt (Pt - P_SINK) = zeta m^2 R T / (2 A^2)."""
    R = float(cb.specific_gas_constant(cb.mass_to_mole(AIR)))
    c = zeta * m * m * R * T_SUPPLY / (2.0 * area * area)
    return 0.5 * (P_SINK + math.sqrt(P_SINK * P_SINK + 4.0 * c))


def _hand_flux(Pt: float, P: float) -> float:
    """Isentropic jet mass flux [kg/m^2/s] from (T_SUPPLY, Pt) to P, ideal gas
    at the supply's gamma. The C++ flux follows cp(T) along the jet; on these
    jets (a few percent of Pt) the two differ by 1.1e-4."""
    X = cb.mass_to_mole(AIR)
    g = float(cb.cp_mass(T_SUPPLY, X)) / float(cb.cv_mass(T_SUPPLY, X))
    R = float(cb.specific_gas_constant(X))
    pr = P / Pt
    return (
        Pt
        / math.sqrt(R * T_SUPPLY)
        * math.sqrt(2.0 * g / (g - 1.0) * (pr ** (2.0 / g) - pr ** ((g + 1.0) / g)))
    )


@pytest.mark.parametrize("m_supply", [M_SUPPLY, 0.3 * M_SUPPLY])
def test_the_valve_drops_and_metering_match_a_hand_calculation(m_supply: float) -> None:
    """Supplied at a fixed flow, every outlet flow is known; each outlet's
    Pt is the closed-form head-loss root at that flow; the inlet sits at the
    reference's; each valve takes the difference, metered by an isentropic
    jet. At 0.3 of the flow the drops are ~3% of Pt, and the valve's zeta is
    the incompressible K with the jet's first-order compressibility. (Lower, the head loss's flow
    reversal blend, exp(-m / 1e-4 kg/s), would show in the hand comparison.)"""
    _, r = _solve(_network("mass", m_supply=m_supply, sink_coupling="total")[0])
    assert r["__success__"], r["__message__"]
    fr = {"o1": 0.2, "o2": 0.3, "o3": 0.5}
    Pt = {o: _hand_Pt(fr[o] * m_supply, z, a) for o, (_, z, a) in PATHS.items()}
    Pt_in = Pt[REFERENCE]
    assert r["sup.Pt"] == pytest.approx(Pt_in, rel=1e-10)
    d = r["__element_diag__"]
    checked_low_mach = False
    for o in PATHS:
        dg = d[f"sp__{o}"]
        assert r[f"{o}.Pt"] == pytest.approx(Pt[o], rel=1e-10), o
        assert dg["valve_dPt"] == pytest.approx(Pt_in - Pt[o], rel=1e-8, abs=1e-6), o
        assert dg["fraction_achieved"] == pytest.approx(fr[o], abs=1e-13)
        assert dg["valve_needs_pump"] == 0.0
        if o == REFERENCE:
            assert dg["is_reference"] == 1.0
            assert "valve_CdA" not in dg
            continue
        m = fr[o] * m_supply
        CdA = m / _hand_flux(Pt_in, Pt[o])
        assert dg["valve_CdA"] == pytest.approx(CdA, rel=1e-3), o
        assert dg["valve_Cd"] == pytest.approx(dg["valve_CdA"] / VALVE_AREA, rel=1e-14)
        assert dg["valve_zeta"] == pytest.approx(1.0 / dg["valve_Cd"] ** 2, rel=1e-14)
        assert dg["valve_choked"] == 0.0
        if dg["valve_dPt"] < 0.05 * Pt_in:
            # The low-Mach K: the drop over the dynamic head at the valve
            # area, with the jet's first-order compressibility: G^2 = 2 rho
            # dPt (1 - 3 x / (2 gamma)), x = dPt / Pt, so zeta = (G A / m)^2
            # is zeta_inc (1 - 3 x / (2 gamma)) + O(x^2).
            X = cb.mass_to_mole(AIR)
            rho = float(cb.density(T_SUPPLY, Pt_in, X))
            g = float(cb.cp_mass(T_SUPPLY, X)) / float(cb.cv_mass(T_SUPPLY, X))
            q = m * m / (2.0 * rho * VALVE_AREA**2)
            x = dg["valve_dPt"] / Pt_in
            zeta_inc = dg["valve_dPt"] / q
            assert dg["valve_zeta"] == pytest.approx(
                zeta_inc * (1.0 - 1.5 * x / g), rel=2.0 * x * x
            )
            checked_low_mach = True
    assert checked_low_mach == (m_supply < M_SUPPLY)


def test_summarize_collects_the_outlets() -> None:
    net, sp = _network("pressure")
    _, r = _solve(net)
    assert r["__success__"], r["__message__"]
    s = sp.summarize(r["__element_diag__"])
    assert s["reference"] == REFERENCE
    assert s["needs_pump"] == []
    assert s["m_dot"] == pytest.approx(r["feed.m_dot"], rel=1e-12)
    assert set(s["outlets"]) == set(PATHS)
    assert s["outlets"]["o2"]["valve_dPt"] > 0.0


# --- the Jacobian ------------------------------------------------------------


@pytest.mark.parametrize("kind", ["mass", "pressure"])
@pytest.mark.parametrize("at", ["root", "perturbed"])
def test_every_row_is_exact(kind: str, at: str) -> None:
    """Every row against central differences, columns scaled by |x|, at the
    root and at a perturbed point (a fraction row is linear in the flows, so
    it is exact anywhere; the reference row likewise in the pressures)."""
    s, r = _solve(_network(kind)[0])
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"], dtype=float)
    if at == "perturbed":
        x = x * (1.0 + 0.03 * np.random.default_rng(903).standard_normal(x.size))
    J = s._residuals_and_jacobian(x)[1].toarray()

    def F(v: np.ndarray) -> np.ndarray:
        return np.asarray(s._residuals_and_jacobian(v, compute_jacobian=False)[0])

    for j, name in enumerate(s.unknown_names):
        h = 1e-6 * max(abs(x[j]), 1e-3)
        e = np.zeros_like(x)
        e[j] = h
        fd = (F(x + e) - F(x - e)) / (2 * h)
        for i in range(J.shape[0]):
            scale = (np.abs(J[i]) * np.abs(x)).max()
            assert J[i, j] * abs(x[j]) == pytest.approx(
                fd[i] * abs(x[j]), rel=1e-5, abs=1e-7 * scale
            ), (i, name)


def test_a_fraction_row_is_the_linear_controller() -> None:
    """The element's own row and its partials: m_i - f_i sum_j m_j, with
    1 - f_i on its own flow and -f_i on each sibling's."""
    net, _ = _network("mass")
    net.resolve_all_topology()
    e = net.elements["sp__o2"]
    assert isinstance(e, SplitterOutletElement)
    from combaero.network import NetworkMixtureState

    st = NetworkMixtureState(P=2e5, Pt=2e5, Tt=T_SUPPLY, m_dot=0.17, Y=AIR)
    res, jac = e.residuals(st, st, flows={"sp__o1": 0.11, "sp__o3": 0.23})
    assert res[0] == pytest.approx(0.17 - 0.3 * (0.17 + 0.11 + 0.23), rel=1e-15)
    assert jac[0] == {"sp__o2.m_dot": 0.7, "sp__o1.m_dot": -0.3, "sp__o3.m_dot": -0.3}
    with pytest.raises(ValueError, match="sibling outlets"):
        e.residuals(st, st, flows={"sp__o1": 0.11})
    ref = net.elements["sp__o1"]
    res, jac = ref.residuals(st, st)
    assert jac[0] == {"sup.Pt": 1.0, "o1.Pt": -1.0}
    assert e.residual_scale_kind == "mdot" and ref.residual_scale_kind == "p"


# --- refusals ----------------------------------------------------------------


@pytest.mark.parametrize(
    ("outlets", "reference", "match"),
    [
        ({"o1": 0.2, "o2": 0.3, "o3": 0.4}, "o1", "sum to 0.9"),
        ({"o1": 0.2, "o2": 0.3, "o3": 0.6}, "o1", "not renormalised"),
        ({"o1": 0.5, "o2": 0.5, "o3": None}, None, "leaving the remainder"),
        ({"o1": None, "o2": None, "o3": 0.5}, None, "at most one remainder"),
        ({"o1": 0.0, "o2": 0.5, "o3": None}, None, r"in \(0, 1\)"),
        ({"o1": 1.0, "o2": 0.5, "o3": None}, None, r"in \(0, 1\)"),
        ({"o1": -0.2, "o2": 0.5, "o3": None}, None, r"in \(0, 1\)"),
        ({"o1": 0.2, "o2": 0.3, "o3": 0.5}, "o9", "not one of its outlets"),
        ({"o1": 1.0}, "o1", "at least two outlets"),
    ],
)
def test_a_bad_set_is_refused(outlets: dict, reference: str | None, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        FlowSplitter("sp", "sup", outlets, reference=reference)


def test_the_sum_tolerance_is_round_off_only() -> None:
    FlowSplitter("sp", "sup", {"a": 0.5, "b": 0.5 + 0.5 * SPLITTER_FRACTION_SUM_TOL}, reference="a")
    with pytest.raises(ValueError, match="not 1"):
        FlowSplitter("sp", "sup", {"a": 0.5, "b": 0.5 + 2.0 * SPLITTER_FRACTION_SUM_TOL}, "a")


def _hand_built(fractions: list[float], refs: list[bool]) -> FlowNetwork:
    net = FlowNetwork()
    _supply(net, "mass")
    ids = [f"e{i}" for i in range(len(fractions))]
    for i, (f, ref) in enumerate(zip(fractions, refs, strict=True)):
        net.add_node(PlenumNode(f"n{i}"))
        net.add_node(PressureBoundary(f"x{i}", Pt=P_SINK, Tt=300.0, Y=AIR))
        net.add_element(SplitterOutletElement(ids[i], "sup", f"n{i}", f, ids, reference=ref))
        net.add_element(
            PressureLossElement(
                f"l{i}", f"n{i}", f"x{i}", correlation=ConstantHeadLoss(1.0, 1e-4), area=1e-4
            )
        )
    return net


@pytest.mark.parametrize(
    ("fractions", "refs", "match"),
    [
        ([0.5, 0.4], [True, False], "sum to 0.9"),
        ([0.5, 0.5], [True, True], "exactly one reference"),
        ([0.5, 0.5], [False, False], "exactly one reference"),
    ],
)
def test_hand_built_outlets_are_checked_as_one_splitter(fractions, refs, match) -> None:
    """The elements check the splitter as a whole, without FlowSplitter."""
    with pytest.raises(ValueError, match=match):
        _hand_built(fractions, refs).resolve_all_topology()


@pytest.mark.parametrize("kind", ["mass", "pressure"])
@pytest.mark.parametrize("start", ["o2", "o3"])
def test_a_wrong_reference_is_swapped_to_the_right_one(kind: str, start: str) -> None:
    """Started from a less resistive path, the first root is refused (o1
    would need a pump); the solve moves the reference to o1, the outlet
    needing the most pressure, and ends at the root a correct pick gives:
    the same flows, fractions and valve drops, in two solves."""
    _, good = _solve(_network(kind, reference="o1")[0])
    _, r = _solve(_network(kind, reference=start)[0])
    assert good["__success__"] and r["__success__"], r["__message__"]
    assert good["__reference_swap_solves__"] == 1
    assert good["__reference_swaps__"]["sp"]["swapped_from"] == []
    assert r["__reference_swaps__"]["sp"] == {
        "reference": "o1",
        "start": start,
        "swapped_from": [start],
    }
    assert r["__reference_swap_solves__"] == 2
    assert f"reference moved from '{start}' to 'o1'" in r["__message__"]
    d, dg = r["__element_diag__"], good["__element_diag__"]
    for o in PATHS:
        e = f"sp__{o}"
        assert r[f"{e}.m_dot"] == pytest.approx(good[f"{e}.m_dot"], rel=1e-9)
        assert d[e]["fraction_achieved"] == pytest.approx(dg[e]["fraction_achieved"], abs=1e-13)
        assert d[e]["valve_dPt"] == pytest.approx(dg[e]["valve_dPt"], rel=1e-7, abs=1e-4)
        assert d[e]["is_reference"] == (1.0 if o == "o1" else 0.0)
        assert d[e]["reference_swaps"] == 1.0
    assert r["__second_law_failures__"] == []


def test_the_first_guess_without_a_remainder_or_reference() -> None:
    """No remainder, no reference: the start is the path with the largest
    seed resistance at its fraction, here o1, so no swap is needed."""
    net, sp = _network("mass", outlets={"o1": 0.2, "o2": 0.3, "o3": 0.5}, reference=None)
    assert sp.reference == "o1"
    _, r = _solve(net)
    assert r["__success__"], r["__message__"]
    assert r["__reference_swap_solves__"] == 1


def test_opting_out_keeps_the_refusal_and_its_message() -> None:
    """auto_reference=False: o3 as given, the root refused, the message
    naming the outlet that should be the reference."""
    _, r = _solve(_network("mass", reference="o3", auto_reference=False)[0])
    assert not r["__success__"]
    assert "make 'o1' the reference" in r["__message__"]
    assert "__reference_swaps__" not in r
    d = r["__element_diag__"]
    assert d["sp__o1"]["valve_needs_pump"] == 1.0
    assert d["sp__o1"]["valve_dPt"] < 0.0


def test_no_feasible_reference_is_refused_with_the_swaps_tried() -> None:
    """One outlet ends above the supply's pressure: with it as the reference
    the supply runs backwards and the other valve pumps; with the other, it
    pumps. The first refusal is returned, its reference restored, the swap
    tried listed, and the outlet to make the reference still named."""
    net = FlowNetwork()
    net.add_node(PressureBoundary("res", Pt=1.5e5, Tt=300.0, Y=AIR))
    net.add_node(PlenumNode("in"))
    net.add_element(OrificeElement("feed", "res", "in", diameter=0.02, Cd=0.8, correlation="fixed"))
    net.add_node(PressureBoundary("hi", Pt=2.0e5, Tt=300.0, Y=AIR))
    net.add_node(PressureBoundary("lo", Pt=1.0e5, Tt=300.0, Y=AIR))
    FlowSplitter("sp", "in", {"hi": 0.5, "lo": None}).add_to(net)
    _, r = _solve(net)
    assert not r["__success__"]
    assert "found no feasible reference (tried 'sp': lo -> hi)" in r["__message__"]
    assert "make 'hi' the reference" in r["__message__"]
    assert r["__reference_swaps__"]["sp"]["reference"] == "lo"
    assert r["__reference_swap_solves__"] >= 2
    assert net.elements["sp__lo"].reference and not net.elements["sp__hi"].reference


# --- the network round trip --------------------------------------------------


def test_to_dict_round_trips_the_splitter() -> None:
    net, _ = _network("pressure")
    _, a = _solve(net)
    _, b = _solve(FlowNetwork.from_dict(net.to_dict()))
    assert a["__success__"] and b["__success__"]
    for o in PATHS:
        assert b[f"sp__{o}.m_dot"] == pytest.approx(a[f"sp__{o}.m_dot"], rel=1e-12)


# --- the C3X coolant supply through one splitter ------------------------------


def _c3x() -> ModuleType:
    path = Path(__file__).resolve().parents[1] / "examples" / "papers" / "hylton1983_c3x_vane.py"
    spec = importlib.util.spec_from_file_location("c3x_for_splitter_903", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _c3x_coolant(mod: ModuleType, split: str | None) -> tuple[FlowNetwork, dict[str, float]]:
    """The C3X example's 19 coolant paths (each segment's share of its hole,
    a compressible channel to the hole exit plenum) at one supply
    temperature: fed per path by its own MassFlowBoundary (``split`` None,
    as the example), or from ONE supply of the summed flow through a
    splitter at the paths' fractions, ``split`` its reference."""
    run = mod.TUNE_RUN
    holes = mod.coolant(run)
    paths = {}
    for g in mod.segments(run):
        h = holes[g["hole"]]
        paths[f"c{g['hole']}_{g['id']}"] = (g["f"] * h["m_dot"], g["f"], h["D"])
    total = sum(m for m, _, _ in paths.values())
    T = sum(holes[i]["T_avg"] * holes[i]["m_dot"] for i in holes) / sum(
        holes[i]["m_dot"] for i in holes
    )
    Y = mod._air()
    net = FlowNetwork()
    if split is not None:
        net.add_node(MassFlowBoundary("coolant", m_dot=total, Tt=T, Y=Y))
    for c, (m, f, D) in paths.items():
        if split is None:
            net.add_node(MassFlowBoundary(f"{c}_in", m_dot=m, Tt=T, Y=Y))
        else:
            net.add_node(PlenumNode(f"{c}_in"))
        net.add_node(PressureBoundary(f"{c}_out", Pt=mod.P_COOLANT_EXIT, Tt=T, Y=Y))
        net.add_element(
            ChannelElement(
                c,
                from_node=f"{c}_in",
                to_node=f"{c}_out",
                length=mod.SPAN,
                diameter=D * math.sqrt(f),
                Dh=D,
                roughness=0.0,
                regime="compressible",
            )
        )
    if split is not None:
        FlowSplitter(
            "split",
            "coolant",
            {f"{c}_in": m / total for c, (m, _, _) in paths.items()},
            reference=f"{split}_in",
        ).add_to(net)
    return net, {c: m for c, (m, _, _) in paths.items()}


def test_the_c3x_coolant_through_one_splitter_gives_the_per_hole_flows() -> None:
    """One supply and a splitter at the measured fractions reproduce the
    example's per-path flows and the pressures its per-path boundaries float
    to, to solver precision; the reference is the path needing the most."""
    mod = _c3x()
    net, flows = _c3x_coolant(mod, None)
    _, a = _solve(net)
    assert a["__success__"], a["__message__"]
    worst = max(flows, key=lambda c: a[f"{c}_in.Pt"])
    _, b = _solve(_c3x_coolant(mod, worst)[0])
    assert b["__success__"], b["__message__"]
    assert b["__reference_swap_solves__"] == 1
    # Started from the least demanding path instead, the solve finds the same
    # reference in one swap, and the same root.
    easiest = min(flows, key=lambda c: a[f"{c}_in.Pt"])
    _, w = _solve(_c3x_coolant(mod, easiest)[0])
    assert w["__success__"], w["__message__"]
    assert w["__reference_swaps__"]["split"]["reference"] == f"{worst}_in"
    assert w["__reference_swap_solves__"] == 2
    for c in flows:
        assert w[f"{c}_in.Pt"] == pytest.approx(b[f"{c}_in.Pt"], rel=1e-8), c
    assert len(flows) == 19
    for c, m in flows.items():
        assert a[f"{c}.m_dot"] == pytest.approx(m, rel=1e-12), c
        assert b[f"{c}.m_dot"] == pytest.approx(m, rel=1e-10), c
        assert b[f"{c}_in.Pt"] == pytest.approx(a[f"{c}_in.Pt"], rel=1e-8), c
        valve = b["__element_diag__"][f"split__{c}_in"]["valve_dPt"]
        assert valve == pytest.approx(a[f"{worst}_in.Pt"] - a[f"{c}_in.Pt"], rel=1e-7, abs=1e-3)
        assert valve >= -1e-6
    assert b["coolant.Pt"] == pytest.approx(a[f"{worst}_in.Pt"], rel=1e-8)


# --- the collector: many inlets against a hand balance -----------------------


def _species(**mole: float) -> list[float]:
    X = [0.0] * N_SPECIES
    for name, x in mole.items():
        X[cb.species_index_from_name(name)] = x
    s = sum(X)
    return list(cb.mole_to_mass([v / s for v in X]))


#: (Tt [K], Y, supplied m_dot [kg/s]) of each stream into the collector:
#: ambient air, hot combustion products, warm nitrogen, superheated steam.
STREAMS = [
    (300.0, AIR, 0.30),
    (1900.0, _species(N2=0.72, CO2=0.09, H2O=0.17, O2=0.02), 0.12),
    (800.0, _species(N2=1.0), 0.20),
    (550.0, _species(H2O=1.0), 0.05),
]


def _hand_mix() -> tuple[float, list[float]]:
    """The adiabatic mix by hand: Y = sum m_i Y_i / M, and the Tt that holds
    sum m_i h_i(Tt_i) = M h(Tt, Y) (absolute enthalpies, root by brentq)."""
    M = sum(m for _, _, m in STREAMS)
    Y = [sum(m * y[k] for _, y, m in STREAMS) / M for k in range(N_SPECIES)]
    H = sum(m * float(cb.h_mass(T, cb.mass_to_mole(y))) for T, y, m in STREAMS)
    X = cb.mass_to_mole(Y)
    T = brentq(lambda t: M * float(cb.h_mass(t, X)) - H, 200.0, 3000.0, xtol=1e-12)
    return T, Y


@pytest.mark.parametrize("collector", ["plenum", "merge_chamber"])
def test_a_many_inlet_collector_closes_a_hand_energy_and_species_balance(collector: str) -> None:
    """Four streams 300-1900 K of four compositions, each fed at a fixed flow
    (the first through a duct, the others through orifices) into a plenum or
    a merge chamber (the duct its main inlet, the others side streams), out
    through a duct to ambient."""
    net = FlowNetwork()
    main = "f0"
    if collector == "plenum":
        net.add_node(PlenumNode("mix"))
    else:
        net.add_node(MomentumChamberNode("mix", main_inlet=main))
    for i, (T, Y, m) in enumerate(STREAMS):
        net.add_node(MassFlowBoundary(f"s{i}", m_dot=m, Tt=T, Y=Y))
        if f"f{i}" == main:
            # The main inlet: a duct of the chamber's bore.
            net.add_element(
                ChannelElement(main, f"s{i}", "mix", length=0.5, diameter=0.1, roughness=0.0)
            )
            continue
        net.add_element(
            OrificeElement(f"f{i}", f"s{i}", "mix", diameter=0.03, Cd=0.7, correlation="fixed")
        )
    net.add_node(PressureBoundary("amb", Pt=1.0e5, Tt=300.0, Y=AIR))
    net.add_element(ChannelElement("out", "mix", "amb", length=0.5, diameter=0.1, roughness=0.0))
    _, r = _solve(net)
    assert r["__success__"], r["__message__"]
    T, Y = _hand_mix()
    assert r["mix.Tt"] == pytest.approx(T, rel=1e-9)
    for k in range(N_SPECIES):
        assert r[f"mix.Y[{k}]"] == pytest.approx(Y[k], abs=1e-12), cb.species_name(k)
    assert r["out.m_dot"] == pytest.approx(sum(m for _, _, m in STREAMS), rel=1e-12)
    led = r["__ledger__"]
    assert abs(led["energy_closure_rel"]) < 1e-10
    assert abs(led["mass_closure_rel"]) < 1e-12


def test_a_stalled_first_solve_also_moves_the_reference() -> None:
    """A wrong reference can stall the cold start before it reaches the
    refused root (stress splitter seed 9, incompressible, fixed supply flow,
    started from the least demanding path). The best iterate already shows
    the metered outlets flowing against their drop; the reference moves as
    at a root, and the second solve succeeds."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))
    from networks import FAMILIES

    net, _ = FAMILIES["splitter"](9)
    outs = [e for e in net.elements.values() if isinstance(e, SplitterOutletElement)]
    need = {
        e.id: net.get_downstream_elements(e.to_node)[0]._seed_resistance(1.2, 1e5, 0.1)
        * e.fraction**2
        for e in outs
    }
    easiest = min(need, key=lambda k: need[k])
    for e in outs:
        e.reference = e.id == easiest
    s = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s._in_reference_swap = True  # the plain solve, no swap: it stalls
        stalled = s.solve()
        s._in_reference_swap = False
        r = s.solve()
    assert not stalled["__converged__"]
    assert r["__success__"], r["__message__"]
    assert r["__reference_swaps__"]["sp"]["swapped_from"] == [net.elements[easiest].to_node]
