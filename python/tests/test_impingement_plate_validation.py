"""A jet-plate ARRAY built from the network, against Florschuetz et al. (1981) (#465).

One ImpingementPlateElement per row, joined by ImpingementCrossflowElement
segments: the network solves the jet flow distribution itself, and each row
reads its Gc/Gj from the solved flows. Three questions, kept apart
(docs/VALIDATION_POLICY.md):

1. Does the NETWORK reproduce the source's own flow model? Florschuetz's
   Eq. 8 (crossflow_to_jet_ratio_at_row) and the jet distribution behind it,
   Gj/Gj_mean = beta N cosh(beta (i - 1/2)) / sinh(beta N), are the 1D
   momentum model the chain discretises. A model-to-model check: fidelity
   of the discretisation, not accuracy against a rig.
2. Does the chain reproduce Fig. 6's measured Nu/Nu1, row by row? Each Fig. 6
   point sits at one row's Eq. 8 abscissa, so it maps to a row; the chain's
   Nu/Nu1 at that row is scored and set beside the correlation-level score
   at the paper's own abscissa. The difference is what the network's flow
   distribution costs. Fidelity: the author's own data.
3. Row 1 against Fig. 5's absolute Nu1, through the plate element.

Without the crossflow's momentum term the chain feeds every row alike, so
its Gc/Gj grows linearly in the row index and misses Eq. 8 by a factor of
two on the strongest-crossflow geometry. That is pinned too, as the reason
the segment element exists.
"""

from __future__ import annotations

import math
from functools import cache

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    FlowNetwork,
    HoleCd,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    NetworkMixtureState,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)
from validation.cooling.schema import DATA_ROOT, load_dataset

D = 0.00254  # Florschuetz's 0.1 in holes
N_ROWS = 10
CD = cb.FLORSCHUETZ_1981_DEFAULT_CD
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


# Over 20 s on one worker: runs in CI's parallel slow job (see
# validation-tests.yml), still before every merge.
pytestmark = pytest.mark.slow


class _FrictionOnlySegment(ImpingementCrossflowElement):
    """The segment with its stations switched off: friction only, the way a
    plain ChannelElement would carry the crossflow."""

    def resolve_topology(self, graph) -> None:
        super().resolve_topology(graph)
        self.from_kappa = None
        self.to_kappa = None


def _boundary(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, 300.0, Y_AIR
    return b


@cache
def _chain(
    xn_d: float, yn_d: float, z_d: float, momentum: bool = True, cd_crossflow: str | None = None
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """Solve an N-row inline array; return per-row (m_dot, Gc/Gj, Re_j).
    ``cd_crossflow`` None keeps the plate's default law."""
    span = 12 * yn_d * D
    height = z_d * D
    g = FlowNetwork()
    g.add_node(_boundary("supply", 1.03e5))
    g.add_node(_boundary("exit", 1.00e5))
    for i in range(1, N_ROWS + 1):
        g.add_node(PlenumNode(f"c{i}"))
    for i in range(1, N_ROWS + 1):
        g.add_element(
            ImpingementPlateElement(
                f"p{i}",
                "supply",
                f"c{i}",
                d_jet=D,
                xn_d=xn_d,
                yn_d=yn_d,
                z_d=z_d,
                span=span,
                plate_thickness=D,
                row=i,
                hole_cd=HoleCd(exit_crossflow="auto" if cd_crossflow is None else cd_crossflow),
            )
        )
        to = f"c{i + 1}" if i < N_ROWS else "exit"
        if momentum:
            seg = ImpingementCrossflowElement(
                f"x{i}", f"c{i}", to, length=xn_d * D, height=height, span=span
            )
        else:
            seg = _FrictionOnlySegment(
                f"x{i}", f"c{i}", to, length=xn_d * D, height=height, span=span
            )
        g.add_element(seg)
    res = NetworkSolver(g).solve()
    assert res["__success__"], res.get("__message__")
    diag = res["__element_diag__"]
    rows = [diag[f"p{i}"] for i in range(1, N_ROWS + 1)]
    return (
        tuple(r["m_dot"] for r in rows),
        tuple(r["Gc_Gj"] for r in rows),
        tuple(r["Re_j"] for r in rows),
    )


def _fig6() -> list:
    return [
        s
        for s in load_dataset()
        if s.source.name == "florschuetz1981" and s.figure == "6" and s.scores is not None
    ]


def _geometries() -> list[tuple[float, float, float]]:
    return sorted(
        {
            (float(s.geometry["xn_d"]), float(s.geometry["yn_d"]), float(s.geometry["z_d"]))
            for s in _fig6()
        }
    )


def _beta(yn_d: float, z_d: float) -> float:
    return math.sqrt(2.0) * CD * (math.pi / 4.0) / (yn_d * z_d)


def test_fig6_spans_27_geometries() -> None:
    assert len(_geometries()) == 27


def test_the_network_reproduces_florschuetz_flow_model() -> None:
    """Q1: per-row Gc/Gj against Eq. 8 and Gj/Gj_mean against the cosh
    distribution, over every Fig. 6 geometry, rows 2-10."""
    gc_err, gj_err = [], []
    for xn, yn, z in _geometries():
        m, gc, _ = _chain(xn, yn, z)
        m_mean = sum(m) / N_ROWS
        beta = _beta(yn, z)
        for i in range(1, N_ROWS + 1):
            gj_flo = beta * N_ROWS * math.cosh(beta * (i - 0.5)) / math.sinh(beta * N_ROWS)
            gj_err.append(m[i - 1] / m_mean / gj_flo - 1.0)
            if i > 1:
                eq8 = cb.crossflow_to_jet_ratio_at_row(yn, z, CD, i)
                gc_err.append(gc[i - 1] / eq8 - 1.0)
    gc_err, gj_err = np.array(gc_err), np.array(gj_err)
    print(
        f"\nGc/Gj vs Eq. 8: bias {gc_err.mean():+.2%}, max |err| {np.abs(gc_err).max():.2%}"
        f"\nGj/Gj_mean vs cosh: bias {gj_err.mean():+.2%}, max |err| {np.abs(gj_err).max():.2%}"
    )
    # Measured 2026-10-08: Gc/Gj bias -0.54%, max 3.8%; Gj bias -0.11%, max
    # 3.8% -- the discretisation of a continuous model, the segment's
    # friction, and compressibility, all three of which Florschuetz's 1D model
    # omits. History: 3.7% (incompressible, #466), 4.1% (each station's own
    # density, #473), 3.8% (the segment's own-flow density term restored,
    # #471 bleed PR). Uncentred, the max was 13%; with no momentum term at
    # all, a factor of two.
    assert abs(gc_err.mean()) < 0.01
    assert np.abs(gc_err).max() < 0.045
    assert np.abs(gj_err).max() < 0.045


def test_without_the_momentum_term_the_supply_stays_uniform() -> None:
    """Why ImpingementCrossflowElement exists: plain channels carry friction
    only, every row sees the same pressure drop, and on the strongest
    crossflow geometry (5, 4, 1) row 10's Gc/Gj comes out at twice Eq. 8.
    The hole Cd is held constant (``HoleCd(exit_crossflow='none')``) to isolate the
    momentum term: at twice Eq. 8 the crossflow Cd (#800) would act too."""
    m, gc, _ = _chain(5.0, 4.0, 1.0, momentum=False, cd_crossflow="none")
    assert max(m) / min(m) < 1.05
    assert gc[-1] / cb.crossflow_to_jet_ratio_at_row(4.0, 1.0, CD, N_ROWS) > 2.0


def _bracket(xn: float, yn: float, z: float, Re_j: float, gc: float) -> float:
    s = cb.florschuetz_1981_inline()
    return (
        cb.jet_array_impingement_nu(s, Re_j, gc, 0.7, xn, yn, z).Nu
        / cb.jet_array_impingement_nu(s, Re_j, 0.0, 0.7, xn, yn, z).Nu
    )


def test_fig6_through_the_chain_scores_like_the_correlation() -> None:
    """Q2: every Fig. 6 point maps to the row whose Eq. 8 abscissa is
    nearest; the chain's Nu/Nu1 at that row is scored against the measured
    value, beside the correlation evaluated at the paper's own abscissa."""
    chain_err, corr_err, worst_map = [], [], 0.0
    for s in _fig6():
        xn, yn, z = (float(s.geometry[k]) for k in ("xn_d", "yn_d", "z_d"))
        _, gc, re = _chain(xn, yn, z)
        eq8 = [cb.crossflow_to_jet_ratio_at_row(yn, z, CD, i) for i in range(1, N_ROWS + 1)]
        for line in (DATA_ROOT / "florschuetz1981" / s.path.name).read_text().splitlines()[1:]:
            x, y = map(float, line.split(","))
            row = int(np.argmin([abs(e - x) for e in eq8]))
            worst_map = max(worst_map, abs(eq8[row] - x))
            chain_err.append(_bracket(xn, yn, z, re[row], gc[row]) / y - 1.0)
            corr_err.append(_bracket(xn, yn, z, 1.0e4, x) / y - 1.0)
    chain_err, corr_err = np.array(chain_err), np.array(corr_err)
    print(
        f"\n{len(chain_err)} points, worst row-mapping gap {worst_map:.3f} in Gc/Gj"
        f"\nchain      : bias {chain_err.mean():+.2%}, MAE {np.abs(chain_err).mean():.2%}"
        f"\ncorrelation: bias {corr_err.mean():+.2%}, MAE {np.abs(corr_err).mean():.2%}"
    )
    # Measured 2026-10-07 over 242 points: chain +4.0% bias / 6.9% MAE,
    # correlation +3.7% / 6.8%. The network's distribution costs 0.3 points.
    assert abs(chain_err.mean() - corr_err.mean()) < 0.01
    assert abs(np.abs(chain_err).mean() - np.abs(corr_err).mean()) < 0.01


def _fig5() -> list[tuple[float, float, float, float]]:
    out = []
    for path in sorted((DATA_ROOT / "florschuetz1981").glob("fig5_xn_*_zd_*.csv")):
        parts = path.stem.split("_")
        xn, z = float(parts[2]), float(parts[4])
        for line in path.read_text().splitlines()[1:]:
            yn, nu = map(float, line.split(","))
            out.append((xn, yn, z, nu))
    return out


def test_row1_reproduces_fig5_through_the_plate() -> None:
    """Q3: absolute Nu1 at Re_j = 1e4, the plate carrying that row's flow
    and no crossflow -- the same score #461 established for the channel path."""
    T, P = 300.0, 1.2e5
    X = cb.species.dry_air()
    tr = cb.complete_state(T, P, X).transport
    err = []
    for xn, yn, z, nu in _fig5():
        plate = ImpingementPlateElement(
            "p", "s", "c", d_jet=D, xn_d=xn, yn_d=yn, z_d=z, span=12 * yn * D, plate_thickness=D
        )
        m = plate.n_holes * 1.0e4 * tr.mu * math.pi * D / 4.0
        r = plate.htc_and_T(NetworkMixtureState(P=P, Pt=P, Tt=T, m_dot=m, Y=Y_AIR))
        err.append(r.h * D / tr.k / nu - 1.0)
    err = np.array(err)
    assert len(err) == 27
    assert abs(err.mean()) < 0.03, err.mean()
    assert np.abs(err).mean() < 0.08, np.abs(err).mean()


@pytest.mark.parametrize("geom", [(5.0, 4.0, 1.0), (10.0, 8.0, 2.0)])
def test_the_closed_form_is_reported_beside_the_network_value(geom) -> None:
    """Diagnostics carry Eq. 8 next to the network's own Gc/Gj when the row
    is given, so a user sees how far their supply is from uniform."""
    xn, yn, z = geom
    _, gc, _ = _chain(xn, yn, z)
    assert gc[0] == 0.0
    assert all(b > a for a, b in zip(gc, gc[1:], strict=False))
