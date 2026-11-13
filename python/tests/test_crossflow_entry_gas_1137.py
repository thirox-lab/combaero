"""A crossflow segment's reservoir entry takes the reservoir's gas (#1137).

``entry_K`` makes the from-node a reservoir: the segment carries the entry
drop (1 + K_in) m^2 / (2 rho A_to^2) (C++ ``channel_entry_drop``). Its density
was the TO-node's gas, which is a mixed or burnt gas when that node also takes
another stream. The duct carries the reservoir's gas: its T0 and composition,
at the duct's static P at the to-end. Both regimes read the term (the
incompressible row, the compressible duct ends' K/(1 + K) share), and the
incompressible chamber-exit head moves with it so the two heads still cancel
to K q (#1139).

Idelchik (1966), Eq. 1-64/1-66 (PDF p. 38, printed 30), defines a resistance
coefficient on the dynamic head of the fluid in the section it refers to, at
that fluid's own operating density: an entrance coefficient (diagram 3-1)
refers to the duct, so its loss is K_in m^2 / (2 rho_duct A^2). The test of
that below is a check of the CONVENTION (the density the coefficient is
applied at), not an accuracy score: no data are compared.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    CrossflowSegmentElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y = list(cb.mole_to_mass(cb.species.dry_air()))
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
A_DUCT = 2.0e-4
K_IN = 0.5
P0 = 1.0e5


def _merge(
    regime: str,
    dp: float,
    node: str = "plenum",
    T_hot: float = 900.0,
    rev: bool = False,
    fed: bool = False,
    fuel: float = 0.0,
) -> FlowNetwork:
    """A 300 K supply through one entry segment (no station) into ``node``
    ('plenum', 'chamber' or 'combustor'), which also takes a hot side stream
    (or CH4 that burns), then a channel to the exit. ``fed``: the reservoir is
    a PlenumNode fed through an orifice, so its state is the solver's."""
    p_sup, p_exit = (P0, P0 + dp) if rev else (P0 + dp, P0)
    g = FlowNetwork()
    g.add_node(PressureBoundary("sup", Pt=p_sup, Tt=300.0, Y=Y))
    g.add_node(PressureBoundary("exit", Pt=p_exit, Tt=300.0, Y=Y))
    res = "sup"
    if fed:
        res = "res"
        g.add_node(PlenumNode("res"))
        g.add_element(
            OrificeElement(
                "o", "sup", "res", diameter=0.03, Cd=0.8, correlation="fixed", in_pipe=False
            )
        )
    if node == "chamber":
        g.add_node(MomentumChamberNode("mc", area=A_DUCT, regime=regime, main_inlet="s1x"))
    elif node == "combustor":
        g.add_node(CombustorNode("mc"))
    else:
        g.add_node(PlenumNode("mc"))
    g.add_element(
        CrossflowSegmentElement(
            "s1x", res, "mc", length=0.05, area=A_DUCT, entry_K=K_IN, regime=regime
        )
    )
    if fuel > 0.0:
        g.add_node(MassFlowBoundary("fuel", m_dot=fuel, Tt=300.0, Y=CH4))
        g.add_element(LosslessConnectionElement("inj", "fuel", "mc"))
    else:
        g.add_node(PressureBoundary("hot", Pt=p_sup, Tt=T_hot, Y=Y))
        g.add_element(OrificeElement("h", "hot", "mc", diameter=0.008, Cd=0.6, correlation="fixed"))
    D = math.sqrt(4.0 * A_DUCT / math.pi)
    g.add_element(
        ChannelElement("tail", "mc", "exit", length=0.05, diameter=1.3 * D, regime=regime)
    )
    return g


def _solve(g: FlowNetwork) -> dict:
    r = NetworkSolver(g).solve()
    assert r["__success__"], r.get("__message__")
    return r


def _gap(build) -> float:
    inc = float(_solve(build("incompressible"))["s1x.m_dot"])
    com = float(_solve(build("compressible"))["s1x.m_dot"])
    return (inc - com) / com


def _segment(regime: str, node: str) -> CrossflowSegmentElement:
    g = _merge(regime, 2000.0, node, fed=True)
    for el in g.elements.values():
        el.resolve_topology(g)
    return g.elements["s1x"]


def _states(T_out: float, Y_out: list[float], m: float) -> tuple:
    st_in = NetworkMixtureState(1.02e5, 1.02e5, Tt=300.0, m_dot=m, Y=Y)
    st_out = NetworkMixtureState(1.0e5, 1.0e5, Tt=T_out, m_dot=m, Y=Y_out)
    return st_in, st_out


#: A to-node composition other than the reservoir's (5 % CH4 by mass).
MIXED = [0.95 * a + 0.05 * f for a, f in zip(Y, CH4, strict=True)]


@pytest.mark.parametrize("m", [0.01, -0.01])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_entry_term_is_the_reservoirs_gas(regime: str, m: float) -> None:
    """The entry term and its partials do not depend on the to-node's gas (its
    T0 or composition), and its T0 partial is the reservoir's. Reversed the
    term is 0 (the reservoir receives at static continuity)."""
    seg = _segment(regime, "plenum")
    terms = {}
    cases = {"cold": (300.0, Y), "hot": (900.0, Y), "hot_mixed": (1800.0, MIXED)}
    for key, (T_out, Y_out) in cases.items():
        st_in, st_out = _states(T_out, Y_out, m)
        t = next(t for t in seg._momentum_terms(st_in, st_out, {}) if t[5] == "entry")
        terms[key] = (t[0].dP, seg._term_partials(t))
    assert terms["hot"] == terms["cold"] == terms["hot_mixed"]
    dP, partials = terms["hot"]
    if m < 0.0:
        assert dP == 0.0 and all(v == 0.0 for v in partials.values())
        return
    assert dP > 0.0
    assert partials["res.Tt"] != 0.0
    assert not any(k.startswith("mc.") and not k.endswith(".P") for k in partials), partials


def test_entry_loss_is_idelchiks_definition() -> None:
    """Incompressible into a moving chamber: the entry term less the chamber
    exit head is the entry loss alone, K_in m^2 / (2 rho A^2) at the duct gas's
    density, the reservoir's T0 at the duct static P (Idelchik Eq. 1-66), with
    the chamber holding a 900 K gas. A convention check, not an accuracy score.
    """
    seg = _segment("incompressible", "chamber")
    m = 0.01
    st_in, st_out = _states(900.0, Y, m)
    kinds = {t[5]: t[0].dP for t in seg._momentum_terms(st_in, st_out, {})}
    rho = float(cb.density(300.0, st_out.P, st_in.X))
    q = m * m / (2.0 * rho * A_DUCT * A_DUCT)
    assert kinds["entry"] == pytest.approx((1.0 + K_IN) * q, rel=1e-13)
    assert kinds["entry"] + kinds["exit_head"] == pytest.approx(K_IN * q, rel=1e-12)


@pytest.mark.parametrize("T_hot", [900.0, 1500.0])
@pytest.mark.parametrize("dp", [500.0, 2000.0, 8000.0])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_hot_merge_plenum_regimes_agree(dp: float, T_hot: float) -> None:
    """A 300 K entry into a plenum that also takes a hot stream: within
    0.02 %; with the plenum's mixed gas in the entry term it was -5.1 %
    (900 K) and -8.4 % (1500 K) at every drop."""
    g = _gap(lambda r: _merge(r, dp, "plenum", T_hot))
    assert abs(g) < 5e-4, g


@pytest.mark.parametrize("node", ["plenum", "chamber"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_hot_merge_compressible_flow_takes_the_carried_gas(node: str) -> None:
    """Compressible, the entry loss K/(1+K) of the entry drop is charged at the
    carried gas: the cold supply's flow rises 4.4 % (chamber) and 2.6 %
    (plenum) at 2 kPa over the mixed-gas charge. Pinned to the measured
    values, and above the old ones."""
    m = float(_solve(_merge("compressible", 2000.0, node))["s1x.m_dot"])
    old, new = {"chamber": (9.5445e-3, 9.9678e-3), "plenum": (9.1450e-3, 9.3848e-3)}[node]
    assert m == pytest.approx(new, rel=1e-4)
    assert m > old * 1.02


def _fd_error(build, x: np.ndarray) -> float:
    s = NetworkSolver(build())
    s._build_x0()
    J = s._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(s, x, floor=1e-10)
    return float(np.max(scaled_error(J, fd, x)))


JAC_CASES = [
    ("plenum", False, 0.0),
    ("plenum", True, 0.0),
    ("chamber", False, 0.0),
    ("chamber", True, 0.0),
    # Burning, forward only: reversed, the exit_in term reads the burnt gas
    # without a composition partial (#1146).
    ("combustor", False, 2e-4),
]


@pytest.mark.parametrize("point", ["root", "perturbed", "m+", "m-"])
@pytest.mark.parametrize(("node", "rev", "fuel"), JAC_CASES)
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_entry_gas_jacobian(regime: str, node: str, rev: bool, fuel: float, point: str) -> None:
    """Fresh solvers (#1142), the reservoir a solved PlenumNode so its gas's
    partials are columns of their own: the entry term's T0 partial lands on
    the reservoir, not on the mixed node. At the root, perturbed, and either
    side of m = 0."""
    dp = 8000.0 if node == "combustor" else 2000.0

    def build() -> FlowNetwork:
        return _merge(regime, dp, node, rev=rev, fed=True, fuel=fuel)

    s = NetworkSolver(build())
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point in ("m+", "m-"):
        x[s.unknown_names.index("s1x.m_dot")] = 1e-7 if point == "m+" else -1e-7
    assert _fd_error(build, x) < 1e-5


@pytest.mark.parametrize("node", ["plenum", "chamber", "combustor"])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_entry_gas_second_law(regime: str, node: str) -> None:
    fuel = 2e-4 if node == "combustor" else 0.0
    r = _solve(_merge(regime, 8000.0, node, fuel=fuel))
    assert not r["__second_law_failures__"], r["__second_law__"]
