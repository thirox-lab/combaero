"""A bleed station's kappa term follows the ARRIVING flow (#1178).

Over a station the side stream's axial velocity is kappa times the arriving
streams' mass-weighted mean velocity, so

    P_a - P_b = (m_b^2 - m_a^2 + kappa (m_a - m_b) (max(m_a, 0) - max(-m_b, 0))) / (rho A^2)

with m_a, m_b the station's from- and to-side flows, signed along the
channel. Forward it is Bassett's kappa (m_a - m_b) m_a, unchanged. Before
#1178 that m_a form was used for every sign: in reversed flow m_a is the
LEAVING stream, which added kappa (m_a - m_b)^2 / (rho A^2) of static rise
along the flow (393 Pa and 180 Pa at s4 and s5 of the #1157 reversed root).

* A bleed chain built backwards (every segment flow negative) is its own
  mirror: the same flows and station pressures. Under the m_a form the
  reversed chain ran a different flow.
* Two reservoirs at one pressure feeding a bleed station from both sides
  (a meeting station, outside Bassett's dividing-flow data) split the hole
  flow equally. Under the m_a form the from-side stream was favoured.
* The analytic Jacobian meets FD on every row of these networks and of the
  #1157 reversed network (FD in ``J_ij |x_j|``, fresh solver, warm starts
  dropped, pressure steps capped off the 1 Pa Fanno floor), either side of
  the kink the term has where a through-flow changes sign (C0, not C1).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from test_station_chain_jacobian_1157 import ROOT_X, TOL, _fixed_x, _fresh, _worst

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
BLEED = cb.STATION_KAPPA_BLEED_BASSETT


def _boundaries(g: FlowNetwork, pts: dict[str, float]) -> None:
    for n, pt in pts.items():
        b = PressureBoundary(n)
        b.Pt, b.Tt, b.Y = pt, 300.0, Y_AIR
        g.add_node(b)


def _holes(g: FlowNetwork, stations: tuple[str, ...], regime: str) -> None:
    for st in stations:
        g.add_element(
            OrificeElement(
                f"h_{st}", st, "sink", Cd=0.7, diameter=0.01, correlation="fixed", regime=regime
            )
        )


def _bleed_chain(mirror: bool, regime: str) -> FlowNetwork:
    """res -> b1 -> b2 -> out with a bleed hole to ``sink`` at b1 and b2. The
    mirror builds the same chain from ``out`` to ``res``: every segment flow
    negative, the physics identical (the pattern of #815's merge chain)."""
    g = FlowNetwork()
    _boundaries(g, {"res": 1.06e5, "out": 1.0e5, "sink": 0.97e5})
    g.add_node(PlenumNode("b1"))
    g.add_node(PlenumNode("b2"))
    kw = {"length": 0.05, "area": 2e-3, "regime": regime}
    if not mirror:
        g.add_element(
            CrossflowSegmentElement(
                "s0", "res", "b1", entry_K=0.0, to_kappa=BLEED, next_seg="s1", **kw
            )
        )
        g.add_element(
            CrossflowSegmentElement(
                "s1", "b1", "b2", from_kappa=BLEED, prev_seg="s0", to_kappa=BLEED,
                next_seg="s2", **kw,
            )
        )  # fmt: skip
        g.add_element(
            CrossflowSegmentElement("s2", "b2", "out", from_kappa=BLEED, prev_seg="s1", **kw)
        )
    else:
        g.add_element(
            CrossflowSegmentElement("s2", "out", "b2", to_kappa=BLEED, next_seg="s1", **kw)
        )
        g.add_element(
            CrossflowSegmentElement(
                "s1", "b2", "b1", from_kappa=BLEED, prev_seg="s2", to_kappa=BLEED,
                next_seg="s0", **kw,
            )
        )  # fmt: skip
        g.add_element(
            CrossflowSegmentElement("s0", "b1", "res", from_kappa=BLEED, prev_seg="s1", **kw)
        )
    _holes(g, ("b1", "b2"), regime)
    return g


def _meeting(regime: str) -> FlowNetwork:
    """Two reservoirs at one Pt feeding a bleed station st from both sides."""
    g = FlowNetwork()
    _boundaries(g, {"left": 1.03e5, "right": 1.03e5, "sink": 1.0e5})
    g.add_node(PlenumNode("st"))
    kw = {"length": 0.03, "area": 1e-3, "regime": regime}
    g.add_element(
        CrossflowSegmentElement(
            "s0", "left", "st", entry_K=0.0, to_kappa=BLEED, next_seg="s1", **kw
        )
    )
    g.add_element(
        CrossflowSegmentElement("s1", "st", "right", from_kappa=BLEED, prev_seg="s0", **kw)
    )
    _holes(g, ("st",), regime)
    return g


def _solve(g: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(g)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r.get("__message__")
    return s, r


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_a_reversed_bleed_chain_is_its_mirror(regime: str) -> None:
    """To round-off in both regimes (measured 1e-13 and 2e-16). Under the m_a
    form the mirror's flows differed by 4.7e-4 (compressible) and 5.1e-4
    (incompressible), its stations by 4.5 Pa."""
    _, fwd = _solve(_bleed_chain(False, regime))
    _, mir = _solve(_bleed_chain(True, regime))
    tol = 1e-9
    for i in range(3):
        assert float(fwd[f"s{i}.m_dot"]) > 0.0
        assert -float(mir[f"s{i}.m_dot"]) == pytest.approx(float(fwd[f"s{i}.m_dot"]), rel=tol)
    for n in ("b1", "b2"):
        assert float(mir[f"{n}.P"]) == pytest.approx(float(fwd[f"{n}.P"]), rel=tol * 1e-2)


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_a_meeting_bleed_station_splits_equally(regime: str) -> None:
    """Both reservoirs at one Pt: the station takes half its hole flow from
    each side (exactly, both regimes). Under the m_a form the from-side flow
    REVERSED (s0 -5.5e-3 kg/s): the right reservoir's stream ran through the
    station into the left one."""
    _, r = _solve(_meeting(regime))
    m0, m1 = float(r["s0.m_dot"]), float(r["s1.m_dot"])
    assert m0 > 0.0 > m1
    tol = 1e-9
    assert -m1 == pytest.approx(m0, rel=tol)
    assert {v["verdict"] for v in r["__second_law__"].values()} == {"pass"}


def _fd_worst(g_builder, regime: str, x: np.ndarray | None = None) -> tuple[float, str]:  # noqa: ANN001
    solver, r = _solve(g_builder(regime))
    x0 = np.array(r["__x_solution__"], dtype=float) if x is None else x
    s = NetworkSolver(g_builder(regime))
    s._build_x0()  # resolves the topology: the model a solve uses (#1148)
    return _worst(s, x0)


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("case", ["mirror", "meeting"])
def test_every_row_meets_fd_at_the_root(case: str, regime: str) -> None:
    build = (lambda reg: _bleed_chain(True, reg)) if case == "mirror" else _meeting
    e, row = _fd_worst(build, regime)
    assert e <= TOL, f"{case} {regime}: {row} scaled {e:.2e}"


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("m_dot", [2e-5, -2e-5])
def test_every_row_meets_fd_either_side_of_the_kink(m_dot: float, regime: str) -> None:
    """The meeting root with s1's flow moved to either side of zero: the
    station st is meeting (m_b < 0) or forward (m_b > 0), the kappa term's
    slope in m_b differing by kappa m_a / (rho A^2) across the joint."""
    _, r = _solve(_meeting(regime))
    s = NetworkSolver(_meeting(regime))
    s._build_x0()
    x = np.array(r["__x_solution__"], dtype=float)
    x[s.unknown_names.index("s1.m_dot")] = m_dot
    e, row = _worst(s, x)
    assert e <= TOL, f"s1.m_dot={m_dot:+.0e} {regime}: {row} scaled {e:.2e}"


@pytest.mark.slow
@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
def test_every_row_meets_fd_on_the_1157_reversed_root(regime: str) -> None:
    """The #1157 reversed network at the root this change moved it to
    (``ROOT_X``: every supply segment forward). The pre-#1178 point (s3 a
    meeting station, s4 and s5 reversed) is checked by the #1157 module."""
    solver = _fresh("reversed", "compressible")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = solver.solve(x0=_fixed_x(solver, "reversed", ROOT_X))
    assert r["__converged__"], r.get("__message__")
    s = _fresh("reversed", regime)
    # Holes p1 and p3 sit within a pascal of zero flow at this root
    # (NEAR_ZERO_HOLE_TOL in the #1157 module).
    e, row = _worst(s, np.array(r["__x_solution__"], dtype=float), near_zero_holes=True)
    assert e <= TOL, f"{regime}: {row} scaled {e:.2e}"


def test_the_1157_reversed_stations_rise_as_their_mirror() -> None:
    """At the #1157 reversed point, each reversed station's static rise along
    the flow is the forward kernel's fed in flow orientation (mirror
    symmetry). The m_a form gave +392.7 Pa (s4) and +180.1 Pa (s5) more."""
    from combaero.network.components import _moving_static_T

    s = _fresh("reversed", "compressible")
    x = _fixed_x(s, "reversed")
    s._propagate_states(x)
    m = {f"s{i}x": float(x[s._unknown_indices[f"s{i}x"][0]]) for i in range(6)}
    A = 1.2e-4
    for i in (4, 5):
        a, b = m[f"s{i - 1}x"], m[f"s{i}x"]
        assert a < 0.0 and b < 0.0
        st = s._get_node_state(s.network.nodes[f"s{i}"], x)
        T = _moving_static_T(st, 0.5 * (a + b), A, False).T
        X = list(st.X)
        rev = cb.station_half_drop(a, b, float(st.P), T, X, A, BLEED).dP
        flow = cb.station_half_drop(-b, -a, float(st.P), T, X, A, BLEED).dP
        assert rev == pytest.approx(-flow, rel=1e-12)
