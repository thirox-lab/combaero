"""An in-pipe orifice whose approach pipe is past its sonic limit (#958).

The orifice drives from its approach pipe's static pressure, p1 =
duct_static_state(Tt, Pt, |m| / A_pipe).P (#651). Past the pipe's sonic
limit that is the sonic pressure P*(Tt, Pt), flat in the flux. It was
returned as critical_pressure_ratio(Tt, Pt) Pt with dP/dTt = 0: with a
variable cp P*/Pt moves with Tt, and that ratio comes from a golden-section
search on G(P), good to ~1e-8, whose noise central differences read as a
derivative. A channel choked into a plenum ahead of the orifice (#788's
network) put the orifice row's Tt column -- carried into the air and fuel
flow columns through the combustor's mixing temperature -- off by a factor
of two (7.2 analytic against 3.4 by central difference). P* is now the
sonic point on the isentrope through (Tt, Pt), closed form in T*, with
analytic Tt and Pt partials.
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
    LosslessConnectionElement,
    MassFlowBoundary,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

X_AIR = cb.species.dry_air()
Y_AIR = list(cb.mole_to_mass(X_AIR))
_X_CH4 = [0.0] * cb.num_species()
_X_CH4[cb.species_index_from_name("CH4")] = 1.0
Y_CH4 = list(cb.mole_to_mass(_X_CH4))
PT = 40e5
G_CHOKED = 1e5


@pytest.mark.parametrize("Tt", [400.0, 1900.0])
def test_beyond_choke_the_duct_state_is_the_sonic_point_of_its_isentrope(Tt: float) -> None:
    s = cb.duct_static_state(Tt, PT, G_CHOKED, X_AIR)
    assert s.choked and s.M == 1.0
    a = cb.speed_of_sound(s.T, X_AIR)
    # M = 1 by energy, and the static state's stagnation pressure is Pt.
    assert cb.h_mass(Tt, X_AIR) - cb.h_mass(s.T, X_AIR) == pytest.approx(0.5 * a * a, rel=1e-10)
    assert cb.Pt_from_static(s.P, s.T, 1.0, X_AIR) == pytest.approx(PT, rel=1e-11)
    assert s.rho == pytest.approx(s.P / (cb.specific_gas_constant(X_AIR) * s.T), rel=1e-14)


@pytest.mark.parametrize("Tt", [400.0, 1900.0])
def test_beyond_choke_the_duct_state_partials_match_central_differences(Tt: float) -> None:
    s = cb.duct_static_state(Tt, PT, G_CHOKED, X_AIR)
    assert s.choked
    for arg, step in (("Tt", 1e-4 * Tt), ("Pt", 1e-5 * PT)):

        def at(d: float, arg: str = arg) -> cb.DuctStaticState:
            kw = {"Tt": Tt, "Pt": PT}
            kw[arg] += d
            return cb.duct_static_state(kw["Tt"], kw["Pt"], G_CHOKED, X_AIR)

        sp, sm = at(step), at(-step)
        for q in ("T", "P", "rho"):
            fd = (getattr(sp, q) - getattr(sm, q)) / (2 * step)
            assert getattr(s, f"d{q}_d{arg}") == pytest.approx(fd, rel=1e-6, abs=1e-10), (q, arg)


def test_the_choked_approach_static_meets_the_subsonic_branch() -> None:
    """P(G) is continuous through the sonic flux: the subsonic root just
    below G* lies within its sqrt(1 - G/G*) approach of the choked value."""
    Tt = 1900.0
    choked = _solver_tools.approach_static_pressure(PT, Tt, G_CHOKED, Y_AIR, False)
    s = cb.duct_static_state(Tt, PT, G_CHOKED, X_AIR)
    G_star = s.rho * cb.speed_of_sound(s.T, X_AIR)
    below = _solver_tools.approach_static_pressure(PT, Tt, (1.0 - 1e-8) * G_star, Y_AIR, False)
    assert choked.choked and not below.choked
    assert pytest.approx(choked.P, rel=1e-3) == below.P
    assert below.P > choked.P


def _net(bore: float) -> FlowNetwork:
    """#788's network with a PlenumNode junction: the 1 m x 0.1 m channel
    chokes into the junction, so the in-pipe exit orifice reads its approach
    pipe past its sonic limit at the junction's Pt."""
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("air", m_dot=32.0, Tt=897.0, Y=Y_AIR))
    n.add_node(MassFlowBoundary("fuel", m_dot=1.0, Tt=300.0, Y=Y_CH4))
    n.add_node(PlenumNode("pl"))
    n.add_node(CombustorNode("cmb"))
    n.add_node(PlenumNode("jct"))
    n.add_node(PressureBoundary("out", Pt=20e5, Tt=300.0, Y=Y_AIR))
    n.add_element(LosslessConnectionElement("ea", "air", "pl"))
    n.add_element(OrificeElement("of", "fuel", "pl", Cd=0.6, diameter=0.02, correlation="fixed"))
    n.add_element(LosslessConnectionElement("ep", "pl", "cmb"))
    n.add_element(ChannelElement("ch", "cmb", "jct", length=1.0, diameter=0.1))
    n.add_element(
        OrificeElement("ox", "jct", "out", Cd=0.6, diameter=bore, correlation="fixed", in_pipe=True)
    )
    return n


def test_the_network_jacobian_past_the_approach_pipes_sonic_limit() -> None:
    s = NetworkSolver(_net(0.099))
    s.homotopy_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=30.0)
    assert r["__success__"]
    assert r["ch.mach_out"] == pytest.approx(1.0, abs=1e-9)
    x_root = np.array(r["__x_solution__"])
    rng = np.random.default_rng(958)
    for x in (x_root, x_root * (1.0 + 0.003 * rng.standard_normal(x_root.size))):
        F, J = s._residuals_and_jacobian(x)
        J = J.toarray() if hasattr(J, "toarray") else np.asarray(J)
        Jfd = fd_jacobian(s, x)
        rows = s._row_kinds(x)
        for i in range(J.shape[0]):
            scale = max(np.max(np.abs(J[i])), 1e-12)
            assert np.max(np.abs(J[i] - Jfd[i])) / scale < 1e-5, rows[i]
