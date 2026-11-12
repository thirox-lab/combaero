"""A merge chamber's face state moves with both compositions (#743).

``chamber_merge_face_state`` gave partials in the flows, P and the
temperatures but none in the chamber gas's or the main stream's composition.
A side stream of another gas moves both with every flow, and the main
inlet's row was off by up to 8.8e-3 (CH4 into air), 1.3e-3 (O2), 4.3e-4 (N2);
a temperature-only difference was exact. The chamber closure's own
composition partials, a central difference since #517, are now analytic too,
and a wall on the chamber carries its static T's composition term.
"""

from __future__ import annotations

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    LosslessConnectionElement,
    MomentumChamberNode,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)
from combaero.network.components import SmoothModel

N_SPECIES = len(cb.species.dry_air())
IX = cb.species_index_from_name
AIR = np.array(cb.mole_to_mass(cb.species.dry_air()))


def _X(Y: np.ndarray) -> list[float]:
    return list(cb.mass_to_mole(list(Y)))


@pytest.mark.parametrize("incompressible", [False, True])
@pytest.mark.parametrize("side", ["outlet", "main"])
@pytest.mark.parametrize(("a", "b"), [("CH4", "N2"), ("CO2", "O2")])
def test_face_state_composition_gradient(incompressible: bool, side: str, a: str, b: str) -> None:
    mix = AIR.copy()
    mix[IX("CH4")] = 0.03
    mix /= mix.sum()
    d = np.zeros(N_SPECIES)
    d[IX(a)], d[IX(b)] = 1.0, -1.0
    h = 1e-6

    def face(s: float):
        Ym = AIR + s * h * d if side == "main" else AIR
        Yo = mix + s * h * d if side == "outlet" else mix
        return cb.chamber_merge_face_state(
            0.2, 700.0, _X(Ym), 0.23, 1.5e5, 650.0, _X(Yo), 8.0, 3.14e-4, incompressible
        )

    r = face(0.0)
    base = AIR if side == "main" else mix
    one_sided = base[IX(a)] <= h
    lo = face(0.0) if one_sided else face(-1.0)
    hi = face(1.0)
    div = h if one_sided else 2 * h
    suffix = "_main" if side == "main" else ""
    for field, grad in (("P_face", "dPf_dY"), ("Pt_face", "dPtf_dY")):
        fd = (getattr(hi, field) - getattr(lo, field)) / div
        an = float(np.dot(getattr(r, grad + suffix), d))
        assert an == pytest.approx(fd, rel=2e-6)


def _network(side_gas: str, walled: bool) -> FlowNetwork:
    g = FlowNetwork()
    Ys = [0.0] * N_SPECIES
    Ys[IX(side_gas)] = 1.0
    g.add_node(PressureBoundary("ai", Pt=1.6e5, Tt=700.0))
    g.add_node(PressureBoundary("fi", Pt=1.6e5, Tt=300.0, Y=Ys))
    g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
    surface = ConvectiveSurface(model=SmoothModel(), area=0.05) if walled else None
    g.add_node(MomentumChamberNode("hn", Dh=0.02, surface=surface, main_inlet="oa"))
    g.add_element(
        ChannelElement("oa", "ai", "hn", length=0.3, diameter=0.02, regime="compressible")
    )
    g.add_element(
        ChannelElement("of", "fi", "hn", length=0.3, diameter=0.003, regime="compressible")
    )
    g.add_element(
        ChannelElement("hout", "hn", "ho", length=0.3, diameter=0.02, regime="compressible")
    )
    if walled:
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
        g.add_wall(ThermalWall("w", "hn", "cch", layers=[WallLayer(0.001, 20.0)]))
    return g


def _jacobian_error(s: NetworkSolver, x: np.ndarray) -> float:
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    # RESIDUAL-ROW-SCALING: compare J_ij |x_j|.
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    return float(np.max(np.abs(Js - fds) / scale))


@pytest.mark.parametrize(
    ("side_gas", "walled"), [("CH4", False), ("O2", False), ("N2", False), ("CH4", True)]
)
def test_a_merge_of_two_gases_has_an_exact_jacobian(
    side_gas: str, walled: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    s = NetworkSolver(_network(side_gas, walled))
    r = s.solve()
    assert r["__success__"], r["__message__"]
    x = np.array(r["__x_solution__"])
    # Measured 1e-8 to 9e-8; before: 4.3e-4 (N2) to 8.8e-3 (CH4).
    assert _jacobian_error(s, x) < 1e-6
    if not walled:
        # Falsified: drop the face's composition partials again.
        real = cb.chamber_merge_face_state

        class _NoY:
            def __init__(self, r) -> None:
                self._r = r
                self.dPf_dY = self.dPf_dY_main = self.dPtf_dY = self.dPtf_dY_main = []

            def __getattr__(self, name: str):
                return getattr(self._r, name)

        monkeypatch.setattr(cb, "chamber_merge_face_state", lambda *a: _NoY(real(*a)))
        assert _jacobian_error(s, x) > 1e-4
