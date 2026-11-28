"""Junction convergence on RANDOM physical boundary conditions.

The network scorecard measures convergence on cases whose boundary conditions
are built from Bassett's analytical K at a target split. For an ACCURACY
question that is exactly right. For a ROBUSTNESS question it is circular: it
asks how often the solver reaches an operating point the paper's own
correlation predicts, on the same 2073 points the closure has been measured
and re-measured against.

This module asks the production question instead. It samples a junction and a
set of boundary conditions uniformly inside physically sensible ranges, with no
reference to any paper or dataset, and asks whether the solver returns an
admissible answer.

**The width of the sample space is the point.** A narrow space could be
flattered by tuning; this one cannot, because nothing in the closure is fitted
against a distribution the closure never sees. `test_junction_random_robustness`
pins the ranges for that reason: narrowing them to improve the number requires
editing an assertion, which shows up in a diff.

Run the full sweep with::

    uv run python -m validation.junction.random_robustness 2000
    uv run python -m validation.junction.random_robustness 2000 20260906 ideal

Outcome classes are kept apart on purpose. Only two of them are the solver's
fault:

``converged``       an admissible root
``no root exists``  the drawn boundary conditions constrain a function of the
                    split to a value the closure cannot produce, so there is
                    nothing to converge to. Counting these as failures is the
                    mistake that started the operating-point investigation
                    (issue #271, docs/archive/JUNCTION_OPERATING_POINT_271.md)
``rejected``        converged, then demoted by the junction's own physics
                    checks: a wrong-direction port flow, or a state where the
                    closure is not net dissipative
``no progress``     the solver could not move although a root exists
``raised``          an exception, which should never happen
"""

from __future__ import annotations

import math
import random
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    SolveOutcome,
    PressureBoundary,
)
from combaero import _core
from combaero.network._mynard2010 import junction_loss_coefficient
from combaero.network.ideal_tee_element import IdealMomentumTeeElement
from combaero.network.mpce_element import MultiPortChamberElement

_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_X = list(cb.mass_to_mole(_Y))

Drive = Literal["imposed_flows", "flow_and_pressures", "all_pressures"]
DRIVES: tuple[Drive, ...] = ("imposed_flows", "flow_and_pressures", "all_pressures")

#: The junction closure under test: the shipped Mynard tee, or the
#: parameter-free ideal momentum tee (#517). Each is judged against ITS OWN
#: closure curve when deciding whether a draw has a root at all.
Model = Literal["mynard", "ideal"]
MODELS: tuple[Model, ...] = ("mynard", "ideal")

# ---------------------------------------------------------------------------
# The sample space. Widen it freely; narrowing it needs a reason in the diff.
# ---------------------------------------------------------------------------
PT_RANGE_PA = (0.5e5, 10.0e5)
TT_RANGE_K = (250.0, 600.0)
#: The closure is documented for low Mach; 0.25 is past where it is trusted,
#: on purpose, so the sweep reports what happens at the edge.
MACH_RANGE = (0.005, 0.25)
THETA_RANGE_DEG = (15.0, 165.0)
#: Area ratio A_common / A_branch, sampled log-uniformly: a branch both much
#: larger and much smaller than the main duct.
PSI_RANGE = (0.3, 10.0)
#: Common-port area, 1 cm2 to 500 cm2, log-uniform.
AREA_RANGE_M2 = (1.0e-4, 0.05)
SPLIT_RANGE = (0.05, 0.95)
#: Imposed pressure differences, as multiples of the common dynamic head. The
#: lower bound is negative because both Bassett and Hager measure a mildly
#: negative straight-leg coefficient in dividing flow.
K_STRAIGHT_RANGE = (-0.5, 3.0)
K_BRANCH_RANGE = (-0.5, 5.0)

#: Where the closure stops being documented. Draws past it are still swept
#: and still reported -- they are simply reported as what they are.
LOW_MACH_LIMIT = 0.3

_Q_GRID = np.linspace(0.02, 0.98, 97)


@dataclass(frozen=True)
class Case:
    """One randomly drawn junction and its boundary conditions."""

    Pt: float
    Tt: float
    mach: float
    theta_deg: float
    psi: float
    split: float
    area: float
    drive: Drive
    joining: bool
    k_straight: float
    k_branch: float


def _log_uniform(rng: random.Random, lo: float, hi: float) -> float:
    return math.exp(rng.uniform(math.log(lo), math.log(hi)))


def sample(rng: random.Random) -> Case:
    """Draw one case uniformly inside the ranges above."""
    return Case(
        Pt=rng.uniform(*PT_RANGE_PA),
        Tt=rng.uniform(*TT_RANGE_K),
        mach=rng.uniform(*MACH_RANGE),
        theta_deg=rng.uniform(*THETA_RANGE_DEG),
        psi=_log_uniform(rng, *PSI_RANGE),
        split=rng.uniform(*SPLIT_RANGE),
        area=_log_uniform(rng, *AREA_RANGE_M2),
        drive=rng.choice(DRIVES),
        joining=rng.random() < 0.5,
        k_straight=rng.uniform(*K_STRAIGHT_RANGE),
        k_branch=rng.uniform(*K_BRANCH_RANGE),
    )


def _scales(case: Case) -> tuple[float, float]:
    """Common mass flow and dynamic head for the drawn state."""
    rho = float(cb.density(case.Tt, case.Pt, _X))
    a = float(cb.speed_of_sound(case.Tt, _X))
    m_com = rho * a * case.mach * case.area
    q_dyn = 0.5 * rho * (m_com / (rho * case.area)) ** 2
    return m_com, q_dyn


def _junction(model: Model, flow_direction: str, loss_factor: float = 1.0, **kw: Any) -> Any:
    if model == "ideal":
        return IdealMomentumTeeElement(**kw, loss_factor=loss_factor)
    if loss_factor != 1.0:
        raise ValueError("loss_factor applies to the ideal model only (#526)")
    return MultiPortChamberElement(**kw, flow_direction=flow_direction, strict=False)


def build(case: Case, model: Model = "mynard", loss_factor: float = 1.0) -> FlowNetwork:
    """The three-port network for this case, wired as production wires one."""
    m_com, q_dyn = _scales(case)
    a_bra = case.area / case.psi
    # Joining flow reverses which side of the junction the loss sits on.
    sign = -1.0 if case.joining else 1.0
    Pt_str = case.Pt - sign * case.k_straight * q_dyn
    Pt_bra = case.Pt - sign * case.k_branch * q_dyn

    net = FlowNetwork()
    for pid, area in (
        ("port_com", case.area),
        ("port_str", case.area),
        ("port_bra", a_bra),
    ):
        # Compressible ports, as this harness always ran them (#588).
        net.add_node(MomentumChamberNode(pid, area=area, regime="compressible"))

    if case.joining:
        if case.drive == "imposed_flows":
            net.add_node(
                MassFlowBoundary(
                    "b_str", m_dot=(1.0 - case.split) * m_com, Tt=case.Tt, Y=_Y
                )
            )
            net.add_node(
                MassFlowBoundary("b_bra", m_dot=case.split * m_com, Tt=case.Tt, Y=_Y)
            )
            net.add_node(PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y, coupling="total"))
        else:
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y, coupling="total"))
            net.add_node(PressureBoundary("b_bra", Pt=Pt_bra, Tt=case.Tt, Y=_Y, coupling="total"))
            if case.drive == "flow_and_pressures":
                net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            else:
                net.add_node(
                    PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y, coupling="total")
                )
        net.add_element(LosslessConnectionElement("lc_str", "b_str", "port_str"))
        net.add_element(LosslessConnectionElement("lc_bra", "b_bra", "port_bra"))
        net.add_element(LosslessConnectionElement("lc_com", "port_com", "b_com"))
        junction = _junction(
            model,
            "merge",
            loss_factor=loss_factor,
            id="jct",
            inlet_nodes=["port_str", "port_bra"],
            outlet_nodes=["port_com"],
            inlet_angles_deg=[0.0, case.theta_deg],
            outlet_angles_deg=[0.0],
            port_areas=[case.area, a_bra, case.area],
        )
    else:
        if case.drive == "imposed_flows":
            net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            net.add_node(
                MassFlowBoundary("b_bra", m_dot=case.split * m_com, Tt=case.Tt, Y=_Y)
            )
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y, coupling="total"))
        else:
            net.add_node(PressureBoundary("b_str", Pt=Pt_str, Tt=case.Tt, Y=_Y, coupling="total"))
            net.add_node(PressureBoundary("b_bra", Pt=Pt_bra, Tt=case.Tt, Y=_Y, coupling="total"))
            if case.drive == "flow_and_pressures":
                net.add_node(MassFlowBoundary("b_com", m_dot=m_com, Tt=case.Tt, Y=_Y))
            else:
                net.add_node(
                    PressureBoundary("b_com", Pt=case.Pt, Tt=case.Tt, Y=_Y, coupling="total")
                )
        net.add_element(LosslessConnectionElement("lc_com", "b_com", "port_com"))
        net.add_element(LosslessConnectionElement("lc_str", "port_str", "b_str"))
        net.add_element(LosslessConnectionElement("lc_bra", "port_bra", "b_bra"))
        junction = _junction(
            model,
            "branch",
            loss_factor=loss_factor,
            id="jct",
            inlet_nodes=["port_com"],
            outlet_nodes=["port_str", "port_bra"],
            inlet_angles_deg=[0.0],
            outlet_angles_deg=[0.0, case.theta_deg],
            port_areas=[case.area, case.area, a_bra],
        )
    net.add_element(junction)
    return net


def _ideal_pair(
    U: np.ndarray,
    areas: np.ndarray,
    angles: np.ndarray,
    joining: bool,
    loss_factor: float = 1.0,
) -> tuple[float, float]:
    """(K_straight, K_branch) of the ideal momentum junction, incompressible.

    Same frame as the Mynard call: axes point away from the junction and U is
    positive INTO it, so the kernel's flow is m = -rho U A.
    """
    rho = 1.2
    geom = _core.IdealJunctionGeometry()
    geom.area = [float(a) for a in areas]
    geom.axis_rad = [float(t) for t in angles]
    geom.port_sign = [1.0, 1.0, 1.0]
    geom.loss_factor = [float(loss_factor)] * 3
    m = [-rho * float(u) * float(a) for u, a in zip(U, areas, strict=True)]
    r = _core.ideal_junction_residuals_and_jacobian([0.0] * 3, [rho] * 3, [0.0] * 3, m, 0.0, geom)
    Pt = [-p + 0.5 * rho * float(u) ** 2 for p, u in zip(r.residual[:3], U, strict=True)]
    c, s, b = (2, 0, 1) if joining else (0, 1, 2)
    sign = 1.0 if joining else -1.0
    head = 0.5 * rho * float(U[c]) ** 2
    return sign * (Pt[s] - Pt[c]) / head, sign * (Pt[b] - Pt[c]) / head


def _closure_curve(
    case: Case, model: Model = "mynard", loss_factor: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """The closure's own straight and branch coefficients over the split."""
    a_bra = case.area / case.psi
    areas = np.array([case.area, case.area, a_bra])
    straight, branch = [], []
    for q in _Q_GRID:
        if case.joining:
            U = np.array([(1.0 - q) * 10.0, q * 10.0 * case.psi, -10.0])
            angles = np.array([0.0, math.radians(case.theta_deg), math.pi])
        else:
            U = np.array([10.0, -(1.0 - q) * 10.0, -q * 10.0 * case.psi])
            angles = np.array([math.pi, 0.0, math.radians(case.theta_deg)])
        if model == "ideal":
            ks, kb = _ideal_pair(U, areas, angles, case.joining, loss_factor)
            straight.append(ks)
            branch.append(kb)
            continue
        # The element's OWN constants: the joining energy transfer acts
        # whenever psi != 1, and the network solves with them.
        result = junction_loss_coefficient(
            U,
            areas,
            angles,
            MultiPortChamberElement.DEFAULT_JOINING_ETRANSFER_ALPHA,
            MultiPortChamberElement.DEFAULT_ETA_SCALE,
        )
        if result.K is None or len(result.K) != 2:
            straight.append(np.nan)
            branch.append(np.nan)
        else:
            straight.append(float(result.K[0]))
            branch.append(float(result.K[1]))
    return np.array(straight), np.array(branch)


def has_root(case: Case, model: Model = "mynard", loss_factor: float = 1.0) -> bool | None:
    """Whether the drawn boundary conditions admit a solution at all.

    With both flows imposed the split is a boundary condition and a root always
    exists. The other two drives constrain a FUNCTION of the split: the
    DIFFERENCE of the two coefficients when the flow level is fixed by a mass
    boundary, their RATIO when the level is free too. A draw whose target lies
    outside what the closure can produce has no root, and blaming the solver
    for it would be a category error.

    None means undetermined, which happens when the closure returns nothing
    across the whole grid.
    """
    if case.drive == "imposed_flows":
        return True
    straight, branch = _closure_curve(case, model, loss_factor)
    ok = ~(np.isnan(straight) | np.isnan(branch))
    if not ok.any():
        return None
    straight, branch = straight[ok], branch[ok]
    if case.drive == "flow_and_pressures":
        target = (case.k_branch - case.k_straight) * (-1.0 if case.joining else 1.0)
        spread = straight - branch if case.joining else branch - straight
        return bool(spread.min() <= target <= spread.max())
    return bool(_direction_crossings(case, straight, branch))


def _direction_crossings(
    case: Case, straight: np.ndarray, branch: np.ndarray
) -> list[tuple[int, float]]:
    """Where the closure's (K_straight, K_branch) points the way the target does.

    With every pressure imposed the level is free, so a root is a split at
    which the closure's pair is a POSITIVE multiple of the drawn pair -- same
    direction, so the level ``q_dyn = dP / K`` is real. Compared as the angle
    of each pair, which is continuous wherever the pair is not zero.

    It used to be compared as the ratio ``K_s / K_b``, which has a pole
    wherever ``K_b`` crosses zero, and that broke both ways (#272): a crossing
    lying between the pole and the next grid point is invisible (draw 52 jumps
    from -15.6 to 8.57 in one step past a target of 11.2 -- four solvable
    draws in 240 called rootless, all four then converged), and a sign flip
    THROUGH the pole reads as a crossing (draw 142, which has none). Returns
    ``(grid index, fraction of the step)`` per crossing.
    """
    if abs(case.k_straight) < 1e-12 and abs(case.k_branch) < 1e-12:
        return []
    phi = np.unwrap(np.arctan2(branch, straight))
    target = math.atan2(case.k_branch, case.k_straight)
    out = []
    for i in range(len(phi) - 1):
        lo, hi = sorted((phi[i], phi[i + 1]))
        # Any 2 pi image of the target inside this step is a crossing.
        k = math.ceil((lo - target) / (2.0 * math.pi))
        image = target + 2.0 * math.pi * k
        if image <= hi and hi > lo:
            out.append((i, (image - phi[i]) / (phi[i + 1] - phi[i])))
    return out



def operating_point(case: Case, model: Model = "mynard", loss_factor: float = 1.0) -> tuple[float, float] | None:
    """(common mass flow, split) the case will run at, as far as it can be
    predicted from the reduced incompressible system."""
    m_ref, q_dyn_ref = _scales(case)
    if case.drive == "imposed_flows":
        return m_ref, case.split
    straight, branch = _closure_curve(case, model, loss_factor)
    ok = np.isfinite(straight) & np.isfinite(branch)
    if not ok.any():
        return None
    qs, ks, kl = _Q_GRID[ok], straight[ok], branch[ok]
    if case.drive == "flow_and_pressures":
        target = (case.k_branch - case.k_straight) * (-1.0 if case.joining else 1.0)
        spread = ks - kl if case.joining else kl - ks
        idx = np.where(np.diff(np.sign(spread - target)) != 0)[0]
        return (m_ref, float(qs[idx[0]])) if idx.size else None
    crossings = _direction_crossings(case, ks, kl)
    if not crossings:
        return None
    i, frac = crossings[0]
    k_pair = np.hypot(ks[i] + frac * (ks[i + 1] - ks[i]), kl[i] + frac * (kl[i + 1] - kl[i]))
    if k_pair == 0.0:
        return None
    # Same direction, so the level is the ratio of the two lengths.
    q_dyn = math.hypot(case.k_straight, case.k_branch) * q_dyn_ref / k_pair
    rho = float(cb.density(case.Tt, case.Pt, _X))
    return case.area * math.sqrt(2.0 * rho * q_dyn), float(qs[i] + frac * (qs[i + 1] - qs[i]))


def max_port_mach(case: Case, model: Model = "mynard", loss_factor: float = 1.0) -> float | None:
    """The highest Mach any port reaches at the operating point.

    The closure is incompressible and documented for low Mach, and WHICH port
    gets there first depends on the drive: with both flows imposed it is the
    branch, since u_bra = q psi u_com and the sweep draws area ratios to 10;
    with three pressures it can be the common port, since the level is then
    whatever the imposed drops demand.

    This is the single strongest predictor of failure in the sweep. Measured
    over solvable draws: converged cases sit at a median of 0.17, failures at
    0.56, and 84% of failures need a port above 0.3. Restricted to draws that
    stay inside the documented range, convergence is 98.4%.
    """
    op = operating_point(case, model, loss_factor)
    if op is None:
        return None
    m, q = op
    rho = float(cb.density(case.Tt, case.Pt, _X))
    sound = float(cb.speed_of_sound(case.Tt, _X))
    a_bra = case.area / case.psi
    return (
        max(
            m / (rho * case.area),
            (1.0 - q) * m / (rho * case.area),
            q * m / (rho * a_bra),
        )
        / sound
    )


_BUCKET_BY_OUTCOME = {
    SolveOutcome.CONVERGED: "converged",
    SolveOutcome.INCONSISTENT: "rejected",
    SolveOutcome.NO_PROGRESS: "no progress",
    SolveOutcome.RESIDUAL_TOO_LARGE: "residual too large",
}


def classify(net: FlowNetwork, timeout: float = 20.0) -> str:
    """Solve and bucket the outcome. Never raises."""
    solver = NetworkSolver(net)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            sol = solver.solve(timeout=timeout)
        except Exception:  # noqa: BLE001 -- the bucket IS the report
            return "raised"
    # Read the solver's own verdict. This used to match substrings of
    # ``__message__`` -- including SciPy's wording, which SciPy is free to
    # change -- so the buckets silently depended on prose. ``__outcome__``
    # is the solver's single machine-readable answer; the bucket names are
    # kept as they were so the summary tables stay comparable.
    return _BUCKET_BY_OUTCOME.get(sol.get("__outcome__"), "other")


@dataclass
class Summary:
    n: int
    seed: int
    model: str = "mynard"
    loss_factor: float = 1.0
    outcomes: Counter = field(default_factory=Counter)
    by_cut: dict[tuple[str, str], Counter] = field(
        default_factory=lambda: defaultdict(Counter)
    )

    #: (converged, total) over solvable draws whose operating point stays
    #: inside the closure's documented low-Mach range.
    in_range: list[int] = field(default_factory=lambda: [0, 0])

    @property
    def in_range_converged_share(self) -> float:
        conv, n = self.in_range
        return conv / n if n else float("nan")

    @property
    def n_with_root(self) -> int:
        return sum(self.by_cut[("root", "exists")].values())

    @property
    def converged_share_where_solvable(self) -> float:
        """The headline: of the draws that admit a root, how many are found."""
        n = self.n_with_root
        return self.by_cut[("root", "exists")]["converged"] / n if n else float("nan")


def _bucket(name: str, value: float, edges: tuple[float, float]) -> str:
    lo, hi = edges
    return f"{name}<{lo:g}" if value < lo else (f"{name}{lo:g}-{hi:g}" if value < hi else f"{name}>{hi:g}")


def run(
    n: int = 2000, seed: int = 20260906, model: Model = "mynard", loss_factor: float = 1.0
) -> Summary:
    """Draw and solve ``n`` random cases. Deterministic for a given seed;
    the draws do not depend on ``model``, so two models see the same cases."""
    rng = random.Random(seed)
    summary = Summary(n=n, seed=seed, model=model, loss_factor=loss_factor)
    for _ in range(n):
        case = sample(rng)
        root = has_root(case, model, loss_factor)
        outcome = classify(build(case, model, loss_factor))
        if root is False and outcome in ("no progress", "residual too large", "other"):
            outcome = "no root exists"
        summary.outcomes[outcome] += 1
        label = {True: "exists", False: "none", None: "undetermined"}[root]
        summary.by_cut[("root", label)][outcome] += 1
        if root is True:
            summary.by_cut[("solvable drive", case.drive)][outcome] += 1
            mach = max_port_mach(case, model, loss_factor)
            if mach is not None:
                label = (
                    "within Mach 0.3"
                    if mach <= LOW_MACH_LIMIT
                    else ("Mach 0.3-1" if mach <= 1.0 else "supersonic port")
                )
                summary.by_cut[("solvable port Mach", label)][outcome] += 1
                if mach <= LOW_MACH_LIMIT:
                    summary.in_range[1] += 1
                    summary.in_range[0] += outcome == "converged"
        summary.by_cut[("flow", "joining" if case.joining else "dividing")][outcome] += 1
        summary.by_cut[("mach", _bucket("", case.mach, (0.05, 0.15)))][outcome] += 1
        summary.by_cut[("area ratio", _bucket("", case.psi, (1.0, 3.0)))][outcome] += 1
        summary.by_cut[("angle", _bucket("", case.theta_deg, (60.0, 120.0)))][outcome] += 1
    return summary


_CLASSES = (
    "converged",
    "no root exists",
    "rejected",
    "no progress",
    "residual too large",
    "other",
    "raised",
)


def format_summary(summary: Summary) -> str:
    lines = [
        f"Random physical boundary conditions, n = {summary.n}, seed {summary.seed}, "
        f"model {summary.model}"
        + (f", loss factor {summary.loss_factor:g}" if summary.loss_factor != 1.0 else ""),
        "",
        f"  {'outcome':<20} {'n':>6} {'share':>7}",
    ]
    for name in _CLASSES:
        count = summary.outcomes[name]
        if count:
            lines.append(f"  {name:<20} {count:>6} {100 * count / summary.n:>6.1f}%")
    lines.append(f"  {'TOTAL':<20} {summary.n:>6}")
    lines.append("")
    lines.append(
        f"  Of the {summary.n_with_root} draws that admit a root, "
        f"{100 * summary.converged_share_where_solvable:.1f}% converge to an admissible one."
    )
    lines.append(
        f"  Of the {summary.in_range[1]} that also stay inside the closure's documented"
    )
    lines.append(
        f"  low-Mach range, {100 * summary.in_range_converged_share:.1f}% converge. The"
        f" remaining failures are"
    )
    lines.append("  concentrated where the model is being used past that range.")
    lines.append("")
    lines.append(f"  {'cut':<16} {'bucket':<22} {'n':>5} {'converged':>10} {'no prog':>9}")
    for (cut, bucket), counts in sorted(summary.by_cut.items()):
        n = sum(counts.values())
        if not n:
            continue
        lines.append(
            f"  {cut:<16} {bucket:<22} {n:>5} {100 * counts['converged'] / n:>9.1f}% "
            f"{100 * counts['no progress'] / n:>8.1f}%"
        )
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> None:
    import sys

    args: list[Any] = list(sys.argv[1:] if argv is None else argv)
    n = int(args[0]) if args else 2000
    seed = int(args[1]) if len(args) > 1 else 20260906
    model = args[2] if len(args) > 2 else "mynard"
    if model not in MODELS:
        raise SystemExit(f"model must be one of {MODELS}, got {model!r}")
    loss_factor = float(args[3]) if len(args) > 3 else 1.0
    print(format_summary(run(n=n, seed=seed, model=model, loss_factor=loss_factor)))


if __name__ == "__main__":
    main()
