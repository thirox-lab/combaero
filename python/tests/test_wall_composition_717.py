"""A wall side's h, T_aw and capacity rate move with its composition (#717).

The wall relay carried a side's sensitivity to its upstream temperature,
flows and pressure, not to its composition. Behind a combustor the products'
transport properties move with the fuel/air split, and the walled network's
Jacobian was off by up to 18% (the issue's audit) -- 2e-4 to 1.2e-3 on the
networks below. The C++ results now carry dh/dY and dT_aw/dY (channel_smooth,
the duct static state, the secant cp), and the relay dots them with the
side's upstream dY/dx.
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
import combaero.network.solver as solver_mod
from combaero.network import (
    ChannelElement,
    CombustorNode,
    ConvectiveSurface,
    FlowNetwork,
    LosslessConnectionElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SmoothModel

N_SPECIES = len(cb.species.dry_air())
IX = cb.species_index_from_name


def _products() -> np.ndarray:
    Y = np.zeros(N_SPECIES)
    for name, v in (("N2", 0.72), ("O2", 0.08), ("CO2", 0.11), ("H2O", 0.09)):
        Y[IX(name)] = v
    return Y / Y.sum()


def _X(Y: np.ndarray) -> list[float]:
    return list(cb.mass_to_mole(list(Y)))


# Tangent directions in Y (sum zero), one into a species absent from the
# mixture: its derivative is the one-sided one as it enters.
DIRECTIONS = [("CO2", "O2"), ("H2O", "N2"), ("CH4", "N2")]


def _along(f, Y: np.ndarray, a: str, b: str) -> float:
    d = np.zeros(N_SPECIES)
    d[IX(a)], d[IX(b)] = 1.0, -1.0
    h = 1e-6
    if Y[IX(a)] > h:
        return (f(Y + h * d) - f(Y - h * d)) / (2 * h)
    return (f(Y + h * d) - f(Y)) / h


def _dir(a: str, b: str) -> np.ndarray:
    d = np.zeros(N_SPECIES)
    d[IX(a)], d[IX(b)] = 1.0, -1.0
    return d


@pytest.mark.parametrize(("a", "b"), DIRECTIONS)
def test_channel_smooth_composition_gradient(a: str, b: str) -> None:
    """At fixed mdot, T and P: Re via mu, Pr via cp mu/k, h via k, T_aw via
    cp, rho and the recovery factor."""
    T, P, D, A, m = 1100.0, 1.5e5, 0.02, 3.14e-4, 0.05
    Y = _products()

    def ch(Yv: np.ndarray):
        X = _X(Yv)
        return cb.channel_smooth(T, P, X, m / (cb.density(T, P, X) * A), D, 0.3, flow_area=A)

    r = ch(Y)
    d = _dir(a, b)
    for field, grad, tol in (("h", r.dh_dY, 1e-6), ("T_aw", r.dT_aw_dY, 1e-5)):
        fd = _along(lambda v, f=field: getattr(ch(v), f), Y, a, b)
        assert float(np.dot(grad, d)) == pytest.approx(fd, rel=tol)


@pytest.mark.parametrize("G_frac", [0.15, 0.6, 0.9992])
@pytest.mark.parametrize(("a", "b"), DIRECTIONS[:2])
def test_duct_static_state_composition_gradient(G_frac: float, a: str, b: str) -> None:
    """The implicit energy/entropy balances, including the bend near choke
    where the sonic flux moves with the composition too."""
    T0, P0 = 1500.0, 2e5
    Y = _products()
    X = _X(Y)
    Gs = cb.mass_flux_isentropic(T0, P0, cb.critical_pressure_ratio(T0, P0, X) * P0, X)
    G = G_frac * Gs
    sd = cb.duct_htc_static_state_dY(T0, P0, G, X, cb.duct_htc_static_state(T0, P0, G, X))
    d = _dir(a, b)
    for field, grad in (
        ("T", sd.dT_dY),
        ("P", sd.dP_dY),
        ("rho", sd.drho_dY),
        ("G_eff", sd.dGeff_dY),
    ):
        fd = _along(
            lambda v, f=field: getattr(cb.duct_htc_static_state(T0, P0, G, _X(v)), f), Y, a, b
        )
        assert float(np.dot(grad, d)) == pytest.approx(fd, rel=1e-5, abs=1e-9 * max(1.0, abs(fd)))


@pytest.mark.parametrize(("T1", "T2"), [(1500.0, 400.0), (900.0, 900.0)])
def test_secant_cp_composition_gradient(T1: float, T2: float) -> None:
    Y = _products()
    g = cb.cp_secant_mass_dY(T1, T2, _X(Y))
    for a, b in DIRECTIONS:
        fd = _along(lambda v: cb.cp_secant_mass_and_jacobian(T1, T2, _X(v))[0], Y, a, b)
        assert float(np.dot(g, _dir(a, b))) == pytest.approx(fd, rel=1e-6)


def _network(method: str, wall_on: str) -> FlowNetwork:
    """Air and pure CH4 burn in a combustor; a cooling channel shares a wall
    with the combustor itself or with its exhaust duct."""
    g = FlowNetwork()
    Yf = [0.0] * N_SPECIES
    Yf[IX("CH4")] = 1.0
    node_surface = ConvectiveSurface(model=SmoothModel(), area=0.05) if wall_on == "node" else None
    g.add_node(PressureBoundary("ai", Pt=1.6e5, Tt=700.0))
    g.add_node(PressureBoundary("fi", Pt=1.6e5, Tt=300.0, Y=Yf))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    g.add_node(CombustorNode("hn", method=method, area=1e-3, Dh=0.02, surface=node_surface))
    g.add_node(PlenumNode("hx"))
    g.add_element(
        ChannelElement("oa", "ai", "hn", length=0.3, diameter=0.02, regime="compressible")
    )
    g.add_element(
        ChannelElement("of", "fi", "hn", length=0.3, diameter=0.003, regime="compressible")
    )
    duct_surface = ConvectiveSurface(model=SmoothModel()) if wall_on == "duct" else None
    g.add_element(
        ChannelElement(
            "hout",
            "hn",
            "hx",
            length=0.4,
            diameter=0.02,
            regime="compressible",
            surface=duct_surface,
        )
    )
    g.add_element(LosslessConnectionElement("hx2", "hx", "ho"))
    g.add_node(PressureBoundary("ci", Pt=1.6e5, Tt=400.0))
    g.add_node(PressureBoundary("co", Pt=1.0e5, Tt=300.0))
    g.add_node(PlenumNode("cp"))
    g.add_element(LosslessConnectionElement("co1", "ci", "cp"))
    g.add_element(
        ChannelElement(
            "cch",
            "cp",
            "co",
            length=0.5,
            diameter=0.015,
            regime="compressible",
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    side = "hn" if wall_on == "node" else "hout"
    g.add_wall(ThermalWall("w", side, "cch", layers=[WallLayer(0.001, 20.0)]))
    return g


def _jacobian_error(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j|.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


@pytest.mark.parametrize("method", ["complete", "equilibrium"])
@pytest.mark.parametrize("wall_on", ["node", "duct"])
def test_a_wall_behind_a_combustor_has_an_exact_jacobian(
    method: str, wall_on: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    s = NetworkSolver(_network(method, wall_on))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    # Measured 7e-8 to 8e-8.
    assert _jacobian_error(s, x) < 1e-6
    # Falsified: without the composition relay it is 6e-4 to 8e-4.
    monkeypatch.setattr(
        solver_mod.NetworkSolver, "_wall_side_dQ_dY", staticmethod(lambda *a, **k: None)
    )
    assert _jacobian_error(s, x) > 1e-4
