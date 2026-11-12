"""A ribbed ChannelElement's drop against measured rib friction (#1106).

`validation/cooling/runner.py` scores the correlation; these score the
ELEMENT, whose drop had Han's four-sided Fanning f applied as Darcy and
ignored the wall count -- 59% low on Han's own two-wall data while every
correlation-level score passed.

The measured fbar is recovered from Han's plotted R by definition only (Eq. 7
and Eq. 6 inverted; see rib_element_friction_runner). Bands are pinned on
the scored numbers with headroom, falsified once: f_multiplier 1.2 moves
the pooled fidelity bias to +19%, 0.8 to -21%.
"""

from __future__ import annotations

import math
import statistics

import pytest

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network.components import (
    ChannelElement,
    ConvectiveSurface,
    MassFlowBoundary,
    PressureBoundary,
    RibbedModel,
)
from validation.cooling import rib_element_friction_runner as rer
from validation.cooling import runner, schema


@pytest.fixture(scope="module")
def rows() -> list[rer.Row]:
    return rer.run()


def test_fidelity_on_hans_own_two_wall_data(rows) -> None:
    """Han's own R data (Figs. 4.46-4.48), n = 2: bias -0.9%, MAE 5.9%."""
    errs = [r.predicted / r.measured - 1.0 for r in rows if r.basis == "fidelity"]
    assert len(errs) >= 100
    assert abs(statistics.fmean(errs)) < 0.03
    assert statistics.fmean(abs(e) for e in errs) < 0.10


def test_every_fidelity_series_is_two_wall(rows) -> None:
    """No Han paper measured one or four ribbed walls: those are accuracy."""
    assert {r.n_walls for r in rows if r.basis == "fidelity"} == {2}
    assert any(r.n_walls == 1 and r.basis == "accuracy" for r in rows)


def test_a_measured_two_wall_point_through_a_solved_network() -> None:
    """One of Han's own points (Fig. 4.46, e/D 0.063, P/e 10, W/H 1), solved.

    A mass flow at the point's Re through the ribbed duct into a pressure
    boundary; the solved friction drop, read back as a passage Fanning fbar,
    must land on the measured one within the source's 9% friction
    uncertainty (Han and Park 1988, p. 187).
    """
    series = next(
        s for s in schema.load_dataset() if s.path.name == "fig4.46_R_eD0.063_pe10_wh1.csv"
    )
    rib_set = runner.SETS[series.scores]()
    geom = rer._geometry(series, rib_set)
    pts = schema.load_points(series)
    p = pts[len(pts) // 2]
    r_abs = p.y / runner._normalised_R(rib_set, geom, 1.0)
    f4 = rer._four_sided_from_R(r_abs, geom)
    re = p.x / (geom.e_D * math.sqrt(f4 / 2.0))
    fbar_measured = (f4 * geom.W_H + rer._f_s(re)) / (geom.W_H + 1.0)

    T, P, D, L = 300.0, 2.0e5, 0.05, 1.0
    X = cb.species.dry_air()
    Y = cb.mole_to_mass(X)
    area = math.pi / 4.0 * D * D
    mu = cb.complete_state(T, P, X).transport.mu
    m_dot = re * area * mu / D

    g = FlowNetwork()
    g.add_node(MassFlowBoundary("in", m_dot=m_dot, Tt=T, Y=Y))
    out = PressureBoundary("out")
    out.Pt, out.Tt, out.Y = P, T, Y
    g.add_node(out)
    g.add_element(
        ChannelElement(
            id="ch",
            from_node="in",
            to_node="out",
            length=L,
            diameter=D,
            roughness=0.0,
            surface=ConvectiveSurface(
                area=0.1,
                model=RibbedModel(
                    correlation_set=rib_set,
                    e_D=geom.e_D,
                    p_e=geom.p_e,
                    W_H=geom.W_H,
                    n_ribbed_walls=2,
                ),
            ),
        )
    )
    sol = NetworkSolver(g).solve()
    assert sol["__success__"], sol.get("__message__")

    rho = cb.density(T, P, X)
    q_exit = m_dot * m_dot / (2.0 * rho * area * area)
    dP = sol["in.Pt"] - P - q_exit
    fbar_network = dP * 2.0 * rho * area * area / (4.0 * (L / D) * m_dot * m_dot)
    assert fbar_network == pytest.approx(fbar_measured, rel=0.09)
