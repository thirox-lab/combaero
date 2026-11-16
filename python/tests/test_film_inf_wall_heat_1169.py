"""A film's undiluted gas takes its node's wall heat at the current x (#1169).

``_film_unmix`` mixes a film discharge node's inflows, without the film's
own coolant, into the T_inf the film's effectiveness is defined against
(#694). The node's wall EnergyBoundary is among the energy boundaries it
mixes, so the wall heat into the node is included (the owner's option B:
heat a walled duct gives up before the holes cools the gas reaching them).

The mix used to run before the node's walls were evaluated, so it read the
wall boundary's Q of the previous pass (or of the previous evaluation, or 0
on a fresh solver), with no lagged-read flag; and its relay took no wall
heat, so d(wall Q)/dx was missing from T_inf and from every row the film
reaches. The converged T_inf was right only because another lagged read (an
outlet boundary reading its feeder ahead of the order) forced a repeat pass.
The mix now runs after the walls and relays their heat.

The network (the issue's probe): hot gas into plenum "g" through a walled
duct "h0" (wall "w0" to a coolant duct "c0"), a FilmCoolingElement "film"
discharging into "g", and a filmed wall "w" on the duct "hot" out of "g",
its backside fed from "c0p". A network whose ONLY lagged read is this one
cannot be built: the order puts every flow boundary first, so any outlet
boundary reads its feeder ahead of the order. A single pass (the cap at 1)
isolates it instead: nothing upstream of "g" is read lagged there.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FilmCoolingElement,
    FlowNetwork,
    MassFlowBoundary,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.solver import _KERNEL_WARM_STARTS

# Over 20 s on one worker: runs in CI's parallel slow job (see
# validation-tests.yml), still before every merge.
pytestmark = pytest.mark.slow

Y_AIR = cb.mole_to_mass(cb.species.dry_air())
SHAPED = {"hole_type": "shaped", "lateral_pitch": 6e-3, "area_ratio": 3.5, "breakout_width": 3e-3}
KEY = "__film_inf__film"
T_HOT = 1500.0
PT_HOT = 2.2e5
#: The analytic Jacobian against FD, scaled (``scaled_error``): measured
#: <= 1.3e-7; without the wall heat in T_inf's relay 1.1e-6 at the root and
#: 1.7e-4 at the perturbed point.
JAC_TOL = 1e-6

#: The solved root (repr round-trips each float exactly).
X_ROOT = {
    "c0_in.P": 201069.8781620898,
    "c0_in.Pt": 201069.8781620898,
    "g.P": 219232.47349148977,
    "g.Pt": 219232.47349148977,
    "hp.P": 219057.20276683135,
    "hp.Pt": 219057.20276683135,
    "cp.P": 200655.60353924325,
    "cp.Pt": 200655.60353924325,
    "c0p.P": 200731.3031197055,
    "c0p.Pt": 200731.3031197055,
    "h0.m_dot": 0.463077879028738,
    "c0.m_dot": 0.05,
    "c0o.m_dot": 0.026414580438097808,
    "film.m_dot": 0.006086249093475303,
    "hot.m_dot": 0.4691641281222133,
    "ho.m_dot": 0.46916412812221336,
    "cold.m_dot": 0.023585419561902188,
    "co.m_dot": 0.023585419561902188,
}


def _net() -> FlowNetwork:
    g = FlowNetwork()
    g.add_node(PressureBoundary("hot_in", Pt=PT_HOT, Tt=T_HOT, Y=Y_AIR))
    g.add_node(PressureBoundary("out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cold_out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(MassFlowBoundary("c0_in", m_dot=0.05, Tt=400.0, Y=Y_AIR))
    g.add_node(PressureBoundary("c0_out", Pt=2.0e5, Tt=300.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=2.3e5, Tt=700.0, Y=Y_AIR))
    for n in ("g", "hp", "cp", "c0p"):
        g.add_node(PlenumNode(n))
    area = math.pi * 0.08 * 0.1

    def duct(eid: str, a: str, b: str, d: float) -> ChannelElement:
        surf = ConvectiveSurface(area=area)
        return ChannelElement(eid, a, b, length=0.1, diameter=d, roughness=0.0, surface=surf)

    def orifice(eid: str, a: str, b: str, d: float, **kw: bool) -> OrificeElement:
        return OrificeElement(eid, a, b, Cd=0.8, diameter=d, correlation="fixed", **kw)

    layers = [WallLayer(0.002, 25.0)]
    # The hot feed into g is a walled duct: g receives wall heat.
    g.add_element(duct("h0", "hot_in", "g", 0.06))
    g.add_element(duct("c0", "c0_in", "c0p", 0.03))
    g.add_element(orifice("c0o", "c0p", "c0_out", 0.03, in_pipe=False))
    g.add_wall(ThermalWall(id="w0", element_a="h0", element_b="c0", layers=layers))
    film = FilmCoolingElement(
        "film",
        "cool",
        "g",
        hole_diameter=1e-3,
        holes_per_row=30,
        wall_thickness=2e-3,
        n_rows=2,
        row_pitch=0.01,
        **SHAPED,
    )
    g.add_element(film)
    g.add_element(duct("hot", "g", "hp", 0.08))
    g.add_element(orifice("ho", "hp", "out", 0.07))
    g.add_element(duct("cold", "c0p", "cp", 0.03))
    g.add_element(orifice("co", "cp", "cold_out", 0.03, in_pipe=False))
    wall = ThermalWall(
        id="w",
        element_a="hot",
        element_b="cold",
        layers=layers,
        film="film",
        film_start=0.0,
        film_end=0.05,
    )
    g.add_wall(wall)
    return g


def _fresh() -> NetworkSolver:
    s = NetworkSolver(_net())
    s.network.resolve_all_topology()
    s._build_x0()
    return s


def _x(s: NetworkSolver) -> np.ndarray:
    return np.array([X_ROOT[n] for n in s.unknown_names], dtype=float)


def _cold(s: NetworkSolver) -> None:
    """Drop the kernels' warm starts: F is then a function of x alone."""
    for el in s.network.elements.values():
        for k in _KERNEL_WARM_STARTS:
            if k in vars(el):
                delattr(el, k)


def _F(s: NetworkSolver, x: np.ndarray) -> np.ndarray:
    _cold(s)
    return np.asarray(s._residuals(x), dtype=float)


def test_the_undiluted_gas_takes_the_wall_heat_into_its_node() -> None:
    """T_inf is the hot stream mixed with the wall heat g receives (h0 loses
    3.4 kW to its coolant), at this x: below the 1500 K feed."""
    s = _fresh()
    x = _x(s)
    _F(s, x)
    Q = s._wall_ebs["g"].Q
    assert Q < -1000.0
    m = float(x[s.unknown_names.index("h0.m_dot")])
    mix = _solver_tools.mixer_from_streams_and_jacobians(
        [cb.MassStream(m, T_HOT, PT_HOT, Y_AIR)], Q=Q, fraction=0.0
    )
    T_inf = s._film_inf_states[KEY][0]
    assert T_inf < T_HOT - 1.0
    assert T_inf == float(mix.Tt_mix)


def test_one_pass_reads_the_wall_heat_at_this_x() -> None:
    """A single pass (cap 1) gives T_inf, and F, of x alone: on a fresh
    solver and after an evaluation at another x alike. The mix read the
    previous evaluation's wall Q: 1500 K fresh, 1492.6 K after x2."""
    s = _fresh()
    s._MAX_PROPAGATION_PASSES = 1
    x = _x(s)
    x2 = x.copy()
    for i, n in enumerate(s.unknown_names):
        if n.endswith("m_dot"):
            x2[i] *= 1.3
    with warnings.catch_warnings():
        # One pass does not settle the outlets' lagged reads: warned.
        warnings.simplefilter("ignore", RuntimeWarning)
        F1 = _F(s, x)
        T1 = s._film_inf_states[KEY][0]
        _F(s, x2)
        F3 = _F(s, x)
        T3 = s._film_inf_states[KEY][0]
    assert T1 == T3
    assert np.array_equal(F1, F3)
    # The single pass already holds the settled T_inf.
    full = _fresh()
    _F(full, x)
    T_settled = full._film_inf_states[KEY][0]
    assert T_settled == T1


def _jac_err(s: NetworkSolver, x: np.ndarray) -> float:
    _cold(s)
    J = s._residuals_and_jacobian(x)[1]
    orig = s._residuals

    def F(z: np.ndarray) -> np.ndarray:
        _cold(s)
        return orig(z)

    s._residuals = F
    try:
        fd = fd_jacobian(s, x)
    finally:
        del s._residuals
    return float(scaled_error(J, fd, x).max())


def test_the_jacobian_reaches_the_wall_heat_in_the_undiluted_gas() -> None:
    """FD over every row, at the root and a perturbed point: d(wall Q)/dx
    reaches T_inf, the film's T_aw and the backside duct's rows. Dropping
    the wall heat from T_inf's relay (its value kept) fails the check."""
    x = _x(_fresh())
    xp = x * (1.0 + 0.03 * np.sin(np.arange(len(x)) + 1.0))
    for point in (x, xp):
        assert _jac_err(_fresh(), point) < JAC_TOL

    s = _fresh()
    unmix = s._film_unmix

    def no_wall_relay(nid, pid, stream_info, _walls, *rest):  # noqa: ANN001, ANN202
        return unmix(nid, pid, stream_info, [], *rest)

    s._film_unmix = no_wall_relay
    assert _jac_err(s, xp) > 10.0 * JAC_TOL
