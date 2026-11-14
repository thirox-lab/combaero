"""The #1142 station chain's Jacobian meets FD on every row (#1157).

The network: a 5-row ImpingementPlateElement array fed by a
CrossflowSegmentElement chain src -> s1..s5 -> sink with Bassett bleed
stations, its gap an ImpingementCrossflowElement chain x1..x5 to the exit.
Fresh solvers on it read up to 3.8e-3 (``s2x`` / ``s3.P``, the supply
reversing at s3) against ``_fd_jacobian.fd_jacobian``.

The analytic Jacobian was right; the FD was not. Each flagged row's DUCT
drive (the Fanno kernel's Pt_in - Pt_out between the segment's own end
states, after the station half-drops) sits just above
``FANNO_FLOW_DRIVE_FLOOR`` (1 Pa): s2x 1.06 Pa, x3 1.19 Pa, forward x4
1.28 Pa. Below the floor the C++ kernel blends the flux into a cubic that
is C1 at the floor, so G'' jumps there. ``fd_jacobian`` steps a pressure
column by 1 % of its smallest pressure GAP (0.24 Pa on s3.P, 0.48 Pa at the
Richardson 2h), which straddles that joint: the extrapolation assumes a
smooth h^2 error and amplifies the jump instead. Refining the step, the
s2x / s3.P entry converges to J as h^2 (7e-3 at 0.1 Pa, 8e-5 at 0.03 Pa,
9e-6 at 0.01 Pa central; 1e-8 Richardson at 0.01 Pa). Every other flagged
entry does the same. The step the pressure unknowns show is not the scale the
row curves on (the ``_fd_jacobian`` docstring's "a drop the unknowns do not
show").

So a compressible segment's row is differenced here at a step on its own
duct's scale: at most 1/50 of its drive (the flux ~ sqrt(drive), whose
Richardson error is ~0.2 (h / drive)^4) and 1/5 of its distance to the floor
(h and 2h on one side of the joint). Every other row keeps
``fd_steps``' step: a hole row's inner solve repeats only to round-off, which
a 0.01 Pa step lifts to 1.5e-7. The kernels' warm starts are cleared before
each evaluation, so F is a deterministic function of x (with them, one x read
twice differs by ~1e-11, the noise in a structurally zero column).

Measured on the fixed points, the solved roots and perturbations, both supply
directions: compressible <= 7.1e-8, incompressible <= 1e-10. The fixed points
are the quiet-machine roots of #1157, so the check does not depend on which
root a loaded machine's solve lands on (under load the reversed case lands on
another root; #1157).
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from _fd_jacobian import fd_steps, is_pressure, scaled_error

import combaero as cb
from combaero.network import (
    CrossflowSegmentElement,
    FlowNetwork,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
)

pytestmark = pytest.mark.slow

D = 0.00254
SPAN = 0.122
Z_D = 2.0
XN_D = 5.0
N_ROWS = 5
TOL = 1e-7
#: The kernels' warm starts (``solver._KERNEL_WARM_STARTS``).
WARM_STARTS = ("_last_G", "_last_M_exit")
#: A compressible segment row's step: at most this fraction of its duct
#: drive and of the drive's distance to the floor (module docstring).
DRIVE_FRACTION = 1.0 / 50.0
FLOOR_FRACTION = 1.0 / 5.0
#: Supply-duct area [m^2]: the supply reverses at a station / runs forward.
A_SUP = {"reversed": 1.2e-4, "forward": 1e-3}

#: The roots of #1157 on a quiet machine (compressible), by unknown; each
#: PlenumNode's Pt unknown is its P.
FIXED_X = {
    "reversed": {
        "s1.P": 101947.932693,
        "c1.P": 100485.360575,
        "s2.P": 102656.168413,
        "c2.P": 100441.367459,
        "s3.P": 102632.282588,
        "c3.P": 100343.062678,
        "s4.P": 101866.733778,
        "c4.P": 100203.861532,
        "s5.P": 100819.780339,
        "c5.P": 100064.587175,
        "s0x.m_dot": 0.00724152047452,
        "s1x.m_dot": 0.00443261598258,
        "p1.m_dot": 0.00280890449194,
        "x1.m_dot": 0.00280890449194,
        "s2x.m_dot": 0.000981516565898,
        "p2.m_dot": 0.00345109941668,
        "x2.m_dot": 0.00626000390862,
        "s3x.m_dot": -0.00252431448381,
        "p3.m_dot": 0.00350583104971,
        "x3.m_dot": 0.00976583495833,
        "s4x.m_dot": -0.00551278836964,
        "p4.m_dot": 0.00298847388584,
        "x4.m_dot": 0.0127543088442,
        "s5x.m_dot": -0.00752823457069,
        "p5.m_dot": 0.00201544620105,
        "x5.m_dot": 0.0147697550452,
    },
    "forward": {
        "s1.P": 101318.281054,
        "c1.P": 100395.226145,
        "s2.P": 101458.323421,
        "c2.P": 100370.807626,
        "s3.P": 101605.285434,
        "c3.P": 100316.522969,
        "s4.P": 101758.901872,
        "c4.P": 100224.457101,
        "s5.P": 101918.554335,
        "c5.P": 100084.561247,
        "s0x.m_dot": 0.0641074282322,
        "s1x.m_dot": 0.0618743696744,
        "p1.m_dot": 0.00223305855776,
        "x1.m_dot": 0.00223305855776,
        "s2x.m_dot": 0.0594513924415,
        "p2.m_dot": 0.00242297723286,
        "x2.m_dot": 0.00465603579062,
        "s3x.m_dot": 0.0568151232272,
        "p3.m_dot": 0.00263626921428,
        "x3.m_dot": 0.00729230500491,
        "s4x.m_dot": 0.0539406455196,
        "p4.m_dot": 0.00287447770767,
        "x4.m_dot": 0.0101667827126,
        "s5x.m_dot": 0.0508012421842,
        "p5.m_dot": 0.00313940333535,
        "x5.m_dot": 0.0133061860479,
    },
}


def _boundary(name: str, Pt: float) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, 300.0, cb.mole_to_mass(cb.species.dry_air())
    return b


def _build(supply: str, regime: str) -> FlowNetwork:
    """The #1142 network: bleed-station supply chain, plate rows, gap chain."""
    net = FlowNetwork()
    net.add_node(_boundary("src", 1.03e5))
    net.add_node(_boundary("sink", 1.02e5))
    net.add_node(_boundary("exit", 1.0e5))
    bleed = cb.STATION_KAPPA_BLEED_BASSETT
    dx = XN_D * D
    area = A_SUP[supply]
    for i in range(1, N_ROWS + 1):
        net.add_node(PlenumNode(f"s{i}"))
        net.add_node(PlenumNode(f"c{i}"))
    net.add_element(
        CrossflowSegmentElement(
            "s0x", "src", "s1", length=0.5 * dx, area=area, entry_K=0.0, to_kappa=bleed,
            next_seg="s1x", regime=regime,
        )
    )  # fmt: skip
    for i in range(1, N_ROWS + 1):
        last = i == N_ROWS
        net.add_element(
            CrossflowSegmentElement(
                f"s{i}x", f"s{i}", "sink" if last else f"s{i + 1}",
                length=(0.5 if last else 1.0) * dx, area=area, from_kappa=bleed,
                prev_seg=f"s{i - 1}x", to_kappa=None if last else bleed,
                next_seg=None if last else f"s{i + 1}x", regime=regime,
            )
        )  # fmt: skip
        net.add_element(
            ImpingementPlateElement(
                f"p{i}", f"s{i}", f"c{i}", d_jet=D, xn_d=XN_D, yn_d=4.0, z_d=Z_D, span=SPAN,
                plate_thickness=D, row=i,
            )
        )  # fmt: skip
        net.add_element(
            ImpingementCrossflowElement(
                f"x{i}", f"c{i}", "exit" if last else f"c{i + 1}", length=dx,
                height=Z_D * D, span=SPAN, regime=regime,
            )
        )  # fmt: skip
    return net


def _fresh(supply: str, regime: str) -> NetworkSolver:
    s = NetworkSolver(_build(supply, regime))
    s._build_x0()
    return s


def _fixed_x(s: NetworkSolver, supply: str) -> np.ndarray:
    point = FIXED_X[supply]
    x = []
    for name in s.unknown_names:
        if name.endswith(".Pt"):
            name = name[: -len("Pt")] + "P"
        x.append(point[name])
    return np.array(x, dtype=float)


def _residuals(s: NetworkSolver, x: np.ndarray) -> np.ndarray:
    """F at x with the kernels started cold: a function of x alone."""
    for el in s.network.elements.values():
        for key in WARM_STARTS:
            if key in vars(el):
                delattr(el, key)
    return np.asarray(s._residuals(x))


def _duct_drives(s: NetworkSolver, x: np.ndarray) -> dict[str, float]:
    """Each compressible segment's duct drive |Pt_in - Pt_out| [Pa]: the
    Fanno kernel's, between the segment's own end states (``_duct_ends``)."""
    s._propagate_states(x)
    out = {}
    for eid, el in s.network.elements.items():
        if not isinstance(el, CrossflowSegmentElement) or el.regime != "compressible":
            continue
        st_in = s._get_node_state(s.network.nodes[el.from_node], x)
        st_out = s._get_node_state(s.network.nodes[el.to_node], x)
        st_in.m_dot = st_out.m_dot = float(x[s._unknown_indices[eid][0]])
        terms = el._momentum_terms(st_in, st_out, s._network_flows(el, x))
        src, dst, _, _ = el._duct_ends(st_in, st_out, terms)
        out[eid] = abs(float(src.Pt) - float(dst.Pt))
    return out


def _row_steps(s: NetworkSolver, x: np.ndarray) -> np.ndarray:
    """Per row, the largest pressure step its duct allows (inf: no duct)."""
    floor = cb.FANNO_FLOW_DRIVE_FLOOR
    caps = {
        eid: min(DRIVE_FRACTION * d, FLOOR_FRACTION * abs(d - floor))
        for eid, d in _duct_drives(s, x).items()
    }
    steps = np.full(len(x), np.inf)
    for i, (owner, _) in enumerate(s._row_kinds(x)):
        if not owner.startswith("node"):
            steps[i] = caps.get(owner.rsplit("'", 2)[1], np.inf)
    return steps


def _fd(s: NetworkSolver, x: np.ndarray, J: np.ndarray) -> np.ndarray:
    """Central FD of F at x; pressure columns Richardson-extrapolated, each
    row at ``fd_steps``' step capped by its own duct (``_row_steps``)."""
    h0 = fd_steps(s, x)
    row_cap = _row_steps(s, x)

    def central(j: int, h: float) -> np.ndarray:
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        return (_residuals(s, xp) - _residuals(s, xm)) / (xp[j] - xm[j])

    fd = np.zeros_like(J)
    for j, name in enumerate(s.unknown_names):
        if not is_pressure(name):
            fd[:, j] = central(j, h0[j])
            continue
        # Never below 1e-8 |P|, where round-off takes over (P_STEP_FLOOR).
        h_rows = np.maximum(np.minimum(h0[j], row_cap), 1e-8 * abs(x[j]))
        # Rows this column does not enter are zero at any step.
        h_rows[J[:, j] == 0.0] = h0[j]
        for h in np.unique(h_rows):
            d = (4.0 * central(j, h) - central(j, 2.0 * h)) / 3.0
            rows = h_rows == h
            fd[rows, j] = d[rows]
    return fd


def _worst(s: NetworkSolver, x: np.ndarray) -> tuple[float, str]:
    _residuals(s, x)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        J = s._residuals_and_jacobian(x)[1].toarray()
    e = scaled_error(J, _fd(s, x, J), x)
    i = int(np.argmax(e))
    return float(e[i]), str(s._row_kinds(x)[i][0])


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("supply", ["reversed", "forward"])
def test_every_row_meets_fd_at_the_1157_point(supply: str, regime: str) -> None:
    """All rows at #1157's root (both regimes evaluate the compressible root:
    the check needs a point, not a root)."""
    s = _fresh(supply, regime)
    e, row = _worst(s, _fixed_x(s, supply))
    assert e <= TOL, f"{supply} {regime}: {row} scaled {e:.2e}"


def test_the_reversed_point_puts_s2x_just_above_the_floor() -> None:
    """The point exercises #1157's case: s2x's duct drive within 0.1 Pa
    above the 1 Pa floor, where the default FD step straddles the joint."""
    s = _fresh("reversed", "compressible")
    d = _duct_drives(s, _fixed_x(s, "reversed"))["s2x"]
    assert cb.FANNO_FLOW_DRIVE_FLOOR < d < cb.FANNO_FLOW_DRIVE_FLOOR + 0.1


@pytest.mark.parametrize("regime", ["compressible", "incompressible"])
@pytest.mark.parametrize("m_dot", [0.0, 1e-7, -1e-7])
@pytest.mark.parametrize("segment", ["s2x", "s3x"])
def test_every_row_meets_fd_through_a_segments_zero_flow(
    segment: str, m_dot: float, regime: str
) -> None:
    """The reversed point with s2x (forward, slow) or s3x (the first reversed
    leg) at zero flow and either side of it."""
    s = _fresh("reversed", regime)
    x = _fixed_x(s, "reversed")
    x[s.unknown_names.index(f"{segment}.m_dot")] = m_dot
    e, row = _worst(s, x)
    assert e <= TOL, f"{segment}.m_dot={m_dot:+.0e} {regime}: {row} scaled {e:.2e}"


@pytest.mark.parametrize("supply", ["reversed", "forward"])
def test_every_row_meets_fd_at_the_solved_root_and_off_it(supply: str) -> None:
    """Whatever root this machine's solve lands on, and a perturbation of it
    (30 Pa on pressures, 3 % on flows)."""
    s = _fresh(supply, "compressible")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(_build(supply, "compressible")).solve()
    x0 = np.array(r["__x_solution__"], dtype=float)
    rng = np.random.default_rng(1157)
    xp = x0.copy()
    for j, name in enumerate(s.unknown_names):
        if is_pressure(name):
            xp[j] += 30.0 * rng.normal()
        else:
            xp[j] *= 1.0 + 0.03 * rng.normal()
    for tag, x in (("root", x0), ("perturbed", xp)):
        e, row = _worst(s, x)
        assert e <= TOL, f"{supply} {tag}: {row} scaled {e:.2e}"
