"""
Shared helpers for the network-mode validation adapters.

The network mode tests the *wrapped* junction element (TeeJunctionElement,
MultiPortChamberBase) as it behaves inside a NetworkSolver, not the raw
K math.

Three separating-flow topologies are supported. The same junction element
is exercised under all three; only the boundary set changes.

  IMPOSED_Q  (2x MFB + 1x PB)
        mb_in (MFB, m_in)  ->  port_com -> junction -> port_str -> pb_str (PB)
                                                    \\-> port_bra -> mb_lat (MFB, m_lat)
        Mass flows imposed; pressure datum at the straight outlet.

  THREE_PB  (3x PB)
        pb_com (PB) -> port_com -> junction -> port_str -> pb_str (PB)
                                            \\-> port_bra -> pb_bra (PB)
        Pressures imposed at every boundary; junction determines split.
        BC pressures are derived from the analytical K so the target q is
        a converged root (tests solver robustness, not K accuracy).

  MFB_TWO_PB  (1x MFB + 2x PB)
        mb_in (MFB) -> port_com -> junction -> port_str -> pb_str (PB)
                                            \\-> port_bra -> pb_bra (PB)
        Inlet flow imposed; outlet pressures set the split.
        Production-realistic BC set.
"""

from __future__ import annotations

import math
import time
import warnings
from dataclasses import dataclass
from collections.abc import Callable
from typing import Literal

import combaero as cb
from combaero.network import (
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkSolver,
    PressureBoundary,
)

_DRY_AIR_Y = list(cb.mole_to_mass(cb.species.dry_air()))
_F_C = 0.01  # 100 cm^2 -- standard area for the validation networks
_M_DOT_REF = 0.1  # kg/s -- low-Mach inlet flow
_TT_REF = 300.0  # K
_PT_REF = 1.0e5  # Pa


Topology = Literal["imposed_q", "three_pb", "mfb_two_pb"]
ALL_TOPOLOGIES: tuple[Topology, ...] = ("imposed_q", "three_pb", "mfb_two_pb")


@dataclass
class NetworkResult:
    """Output of a single network-mode evaluation."""

    converged: bool
    K_lateral: float | None = None
    K_straight: float | None = None
    residual_norm: float = math.inf
    wall_time_s: float = 0.0
    message: str = ""
    q_converged: float | None = None  # actual q the solver found, may drift from target


def _x_air() -> list[float]:
    """Dry-air mole fractions (cached call to cb.mass_to_mole)."""
    return list(cb.mass_to_mole(_DRY_AIR_Y))


def _achieved_q(
    net: FlowNetwork,
    sol: dict[str, float],
    *,
    common_node: str,
    lateral_node: str,
    junction_id: str = "jct",
) -> float | None:
    """Lateral mass-flow fraction at the converged state, on the NETWORK axis.

    In `imposed_q` the split is a boundary condition and this returns the
    target. In the pressure-driven topologies it is an OUTCOME: the solve
    settles wherever the model's coefficients reproduce the imposed pressure
    differences, which is generally not the target q (issue #271,
    docs/archive/JUNCTION_OPERATING_POINT_271.md). Comparing an extracted K
    against the paper at the target q is then a comparison at the wrong
    operating point, so every record carries the q it actually reached.

    Callers that re-index q onto a paper's own axis must apply the same
    transform to this value -- see `MPCENetwork`.
    """

    def _mdot_at(node_id: str) -> float | None:
        for elem in net.elements.values():
            if elem.id == junction_id:
                continue
            frm = getattr(elem, "from_node", None)
            to = getattr(elem, "to_node", None)
            if node_id not in (frm, to):
                continue
            key = f"{elem.id}.m_dot"
            if key in sol:
                return abs(float(sol[key]))
        return None

    m_com = _mdot_at(common_node)
    m_lat = _mdot_at(lateral_node)
    if m_com is None or m_lat is None or m_com <= 0.0:
        return None
    return m_lat / m_com


def _extract_K(
    sol: dict[str, float],
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float,
) -> tuple[float, float]:
    """Extract K_lateral and K_straight from a converged solution.

    K = (Pt_common - Pt_branch) / (0.5 * rho_common * u_common^2),
    where Pt is the stagnation pressure at the junction face of each branch.
    """
    P_com = sol[f"{common_node}.P"]
    Pt_com = sol[f"{common_node}.Pt"]
    Pt_str = sol[f"{straight_node}.Pt"]
    Pt_lat = sol[f"{lateral_node}.Pt"]
    X = _x_air()
    rho_com = float(cb.density(_TT_REF, P_com, X))
    u_com = m_dot_ref / (rho_com * area)
    q_dyn_com = 0.5 * rho_com * u_com * u_com
    return (
        (Pt_com - Pt_lat) / q_dyn_com,
        (Pt_com - Pt_str) / q_dyn_com,
    )


# ---------------------------------------------------------------------------
# Network builders
# ---------------------------------------------------------------------------


#: The junction data is gas measured up to Ma ~0.3 (Wang 2014): the harness
#: validates the COMPRESSIBLE port closure, so it pins it rather than taking
#: the library's default regime, which is incompressible until #585 (#588).
_REGIME = "compressible"


def _q_dyn_ref(m_dot: float = _M_DOT_REF, area: float = _F_C, Pt: float = _PT_REF) -> float:
    """Reference dynamic head at the low-Mach inlet (rho*u^2/2)."""
    X = _x_air()
    rho = float(cb.density(_TT_REF, Pt, X))
    u = m_dot / (rho * area)
    return 0.5 * rho * u * u


def build_separating_network_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    m_lateral: float | None = None,
    Pt_ref: float = _PT_REF,
) -> FlowNetwork:
    """IMPOSED_Q skeleton: MFB on common inlet + MFB on lateral outlet + PB on straight.

    Mass flows imposed; PB at the straight sets the pressure datum.
    Caller adds the junction element + per-port loss + lossless wiring.
    """
    if m_lateral is None:
        m_lateral = 0.5 * m_in
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=_TT_REF, Y=Y))
    net.add_node(MassFlowBoundary("mb_lat", m_dot=m_lateral, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_ref, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_in", "mb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    return net


def build_separating_three_pb_skeleton(
    *,
    K_lateral_target: float,
    K_straight_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """THREE_PB skeleton: PB on every branch.

    Pressure drops imply the target K relations: Pt_com - Pt_lat = K_lat * q_dyn
    and Pt_com - Pt_str = K_str * q_dyn. The network is solved blind; if the
    junction element's K differs from the target, the solver finds the q where
    the model's K matches the imposed Pt drops.

    BC pressures are based on the reference dynamic head, which is fixed
    rather than re-derived per (m, A) -- the absolute value drops out of the
    K extraction.
    """
    q_dyn = _q_dyn_ref(Pt=Pt_com)
    Pt_lat = Pt_com - K_lateral_target * q_dyn
    Pt_str = Pt_com - K_straight_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_com", Pt=Pt_com, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_lat, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_com", "pb_com", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    return net


def build_separating_mfb_two_pb_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    K_lateral_target: float,
    K_straight_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """MFB_TWO_PB skeleton: MFB on common inlet + PB on each outlet.

    Inlet flow imposed at m_in; outlet pressures sized from analytical K
    so the target q is a converged root. The most production-realistic BC
    set for steady-state design.
    """
    q_dyn = _q_dyn_ref(m_dot=m_in, Pt=Pt_com)
    Pt_lat = Pt_com - K_lateral_target * q_dyn
    Pt_str = Pt_com - K_straight_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_in", m_dot=m_in, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_lat, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_in", "mb_in", "port_com"))
    net.add_element(LosslessConnectionElement("lc_str", "port_str", "pb_str"))
    return net


def build_joining_imposed_q_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    m_lateral: float | None = None,
    Pt_ref: float = _PT_REF,
) -> FlowNetwork:
    """IMPOSED_Q joining skeleton: MFB on each inlet + PB on common outlet.

    m_in is the combined outflow at the common; m_lateral is the branch inlet
    mass flow; the straight inlet flow is m_in - m_lateral. PB at the common
    outlet sets the pressure datum.
    """
    if m_lateral is None:
        m_lateral = 0.5 * m_in
    m_straight = m_in - m_lateral
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_str", m_dot=m_straight, Tt=_TT_REF, Y=Y))
    net.add_node(MassFlowBoundary("mb_bra", m_dot=m_lateral, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_com", Pt=Pt_ref, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_str", "mb_str", "port_str"))
    net.add_element(LosslessConnectionElement("lc_com", "port_com", "pb_com"))
    return net


def build_joining_three_pb_skeleton(
    *,
    K11_target: float,
    K12_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """THREE_PB joining skeleton: PB on every branch.

    BC pressures derived from analytical K11/K12 so the target q is a
    converged root. Suppliers (str, bra) sit at HIGHER Pt than the common
    collector since K11, K12 > 0 means pressure drops INTO the junction.
    """
    q_dyn = _q_dyn_ref(Pt=Pt_com)
    Pt_str = Pt_com + K11_target * q_dyn
    Pt_bra = Pt_com + K12_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(PressureBoundary("pb_com", Pt=Pt_com, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_bra, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_str", "pb_str", "port_str"))
    net.add_element(LosslessConnectionElement("lc_com", "port_com", "pb_com"))
    return net


def build_joining_mfb_two_pb_skeleton(
    *,
    m_in: float = _M_DOT_REF,
    K11_target: float,
    K12_target: float,
    Pt_com: float = _PT_REF,
) -> FlowNetwork:
    """MFB_TWO_PB joining skeleton: PB on each inlet + MFB on common outlet.

    Common-outlet flow imposed at m_in; supplier pressures sized from
    analytical K11/K12 so target q is a converged root. The joining analog
    of the production-realistic separating mfb_two_pb.
    """
    q_dyn = _q_dyn_ref(m_dot=m_in, Pt=Pt_com)
    Pt_str = Pt_com + K11_target * q_dyn
    Pt_bra = Pt_com + K12_target * q_dyn
    net = FlowNetwork()
    Y = _DRY_AIR_Y
    net.add_node(MassFlowBoundary("mb_com", m_dot=m_in, Tt=_TT_REF, Y=Y))
    net.add_node(PressureBoundary("pb_str", Pt=Pt_str, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(PressureBoundary("pb_bra", Pt=Pt_bra, Tt=_TT_REF, Y=Y, coupling="total"))
    net.add_node(MomentumChamberNode("port_com", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_str", area=_F_C, regime=_REGIME))
    net.add_node(MomentumChamberNode("port_bra", area=_F_C, regime=_REGIME))
    net.add_element(LosslessConnectionElement("lc_str", "pb_str", "port_str"))
    net.add_element(LosslessConnectionElement("lc_com", "port_com", "mb_com"))
    return net


def _extract_K_joining(
    sol: dict[str, float],
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float,
) -> tuple[float, float]:
    """Extract joining flow K coefficients.

    Returns ``(K_lateral, K_straight)`` to match the separating
    ``_extract_K`` return order so callers unpack the same way:
      K_lateral = K12 = (Pt_branch   - Pt_common) / q_dyn_common
      K_straight = K11 = (Pt_straight - Pt_common) / q_dyn_common
    Both > 0 in normal joining flow (suppliers above the chamber Pt).
    """
    P_com = sol[f"{common_node}.P"]
    Pt_com = sol[f"{common_node}.Pt"]
    Pt_str = sol[f"{straight_node}.Pt"]
    Pt_lat = sol[f"{lateral_node}.Pt"]
    X = _x_air()
    rho_com = float(cb.density(_TT_REF, P_com, X))
    u_com = m_dot_ref / (rho_com * area)
    q_dyn_com = 0.5 * rho_com * u_com * u_com
    return (
        (Pt_lat - Pt_com) / q_dyn_com,  # K12
        (Pt_str - Pt_com) / q_dyn_com,  # K11
    )



def mass_flow_for_mach(
    mach: float,
    area: float = _F_C,
    Pt: float = _PT_REF,
    Tt: float = _TT_REF,
) -> float:
    """Mass flow that puts a duct of ``area`` at ``mach``, isentropically.

    Used to place a validation network at a source's own Mach number instead of
    the fixture's fixed low-speed reference. Goes through the library's own
    thermodynamics (`mass_flux_isentropic`) rather than a perfect-gas formula so
    the property model is the one the solver will use.
    """
    X = _x_air()
    gamma = float(cb.isentropic_expansion_coefficient(Tt, X))
    ratio = (1.0 + 0.5 * (gamma - 1.0) * mach * mach) ** (gamma / (gamma - 1.0))
    G = float(cb.mass_flux_isentropic(Tt, Pt, Pt / ratio, X))
    return G * area


def _extract_K_joining_dynamic_pressure(
    sol: dict[str, float],
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float,
) -> tuple[float, float]:
    """Joining K normalised on the common port's TOTAL MINUS STATIC pressure.

    Wang 2014 defines K_13 = (p0_1 - p0_3) / (p0_3 - p_3) and K_23 likewise,
    where 3 is the common branch. The usual `_extract_K_joining` divides by
    1/2 rho u^2 evaluated at a fixed reference flow, which agrees with that only
    in the incompressible limit and only when the flow is at the reference
    level: p0 - p exceeds 1/2 rho u^2 by about M^2/4 in relative terms, so the
    two definitions differ by roughly 0.2% at Mach 0.1 and 9% at Mach 0.6.

    Extracting what the source actually measured keeps the comparison honest at
    every Mach, rather than folding a definitional difference into what looks
    like model error.
    """
    Pt_com = sol[f"{common_node}.Pt"]
    P_com = sol[f"{common_node}.P"]
    q_dyn = Pt_com - P_com
    if q_dyn <= 0.0:
        raise KeyError("common port has no dynamic pressure to normalise on")
    return (
        (sol[f"{lateral_node}.Pt"] - Pt_com) / q_dyn,
        (sol[f"{straight_node}.Pt"] - Pt_com) / q_dyn,
    )


def solve_and_extract(
    net: FlowNetwork,
    *,
    common_node: str,
    straight_node: str,
    lateral_node: str,
    m_dot_ref: float,
    area: float = _F_C,
    flow_direction: Literal["separating", "joining"] = "separating",
    verifier: Callable[[dict[str, float]], bool] | None = None,
    extractor: Callable[..., tuple[float, float]] | None = None,
) -> NetworkResult:
    """Run NetworkSolver, return convergence + K diagnostics."""
    solver = NetworkSolver(net)
    t0 = time.perf_counter()
    try:
        # Suppress the per-call non-convergence UserWarnings -- the scorecard
        # already captures success/residual norm explicitly.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sol = solver.solve()
    except Exception as e:
        return NetworkResult(
            converged=False, message=str(e)[:120], wall_time_s=time.perf_counter() - t0
        )
    wall_time = time.perf_counter() - t0
    converged = bool(sol.get("__success__", False))
    # `__final_norm__` is the key the solver sets; `__residual_norm__` never
    # existed, so every record's residual norm was silently infinity.
    res_norm = float(sol.get("__final_norm__", math.inf))
    if not converged:
        return NetworkResult(
            converged=False,
            residual_norm=res_norm,
            wall_time_s=wall_time,
            message=sol.get("__message__", "")[:120],
        )
    if verifier is not None and not verifier(sol):
        return NetworkResult(
            converged=False,
            residual_norm=res_norm,
            wall_time_s=wall_time,
            message="post-solve verification failed (soft barrier landed off-physical)",
        )
    try:
        if extractor is None:
            extractor = _extract_K_joining if flow_direction == "joining" else _extract_K
        K_lat, K_str = extractor(
            sol,
            common_node=common_node,
            straight_node=straight_node,
            lateral_node=lateral_node,
            m_dot_ref=m_dot_ref,
            area=area,
        )
    except KeyError as e:
        return NetworkResult(
            converged=False,
            residual_norm=res_norm,
            wall_time_s=wall_time,
            message=f"missing solution key {e}",
        )
    return NetworkResult(
        converged=True,
        K_lateral=K_lat,
        K_straight=K_str,
        residual_norm=res_norm,
        wall_time_s=wall_time,
        q_converged=_achieved_q(
            net, sol, common_node=common_node, lateral_node=lateral_node
        ),
    )
