"""Second-law check per element and a control-volume ledger (#545).

Every converged solve is checked twice after the fact:

* **Per element: no flow work manufactured.** Over each element's control
  volume the entropy generation must not be negative::

      sigma = sum_out m s - sum_in m s - integral dQ / T  >=  0

  with ``s = cb.s_mass(Tt, X, Pt)`` (stagnation state; it includes the ideal
  mixing term) and the heat taken at the fluid's static temperature, which
  makes ``sigma`` the INTERNAL entropy generation -- a necessary condition
  that needs no wall temperature. Both sums and the heat integral are C++
  (``cb.entropy_generation``, ``cb.heat_entropy_flux``). For an adiabatic
  two-port of fixed composition sigma is exactly ``-m R ln(Pt_out / Pt_in)``,
  compressible flow included: "total pressure does not rise". Heating at Mach
  M costs total pressure and cooling may raise it (Rayleigh); the heat term
  allows exactly that, through the static temperature.

  Inside an element only its end Mach numbers are known (its own ``mach_in``
  / ``mach_out`` where it reports them), so the heat is taken along a linear
  Mach path and bracketed by constant-Mach paths at either end:
  ``sigma_norm_low <= sigma_norm <= sigma_norm_high``. ``sigma_norm = sigma /
  (R sum_in m)`` is dimensionless; for a two-port it is ``-ln(Pt_out /
  Pt_in)``. A verdict is ``"pass"``, ``"fail"`` (``sigma_norm_high <
  -SIGMA_TOL``: even the most favourable path destroys entropy) or
  ``"not_checked"`` with the reason -- never a silent pass. A stream on the
  incompressible closure widens the band by that closure's resolution,
  order M^4 (``sigma_norm_resolution``; see ``_sigma``).

* **The network: mass and energy close over its control volume.** The
  streams crossing at the flow boundaries (``NetworkNode.is_flow_boundary``),
  the heat given to the network, and the heat leaving it (``Q_wall_out``, a
  wall's heat landing on a boundary; ``Q_loss``, a wall's heat to a HeatLoss
  side's external temperature, #970) or withheld (``Q_withheld``, heat a
  stagnant node cannot take). Enthalpy is absolute (formation included), so a
  combustor closes the balance as it is.

Both demote a converged solve that fails them to INCONSISTENT: the ledger
since #552, an element's verdict since #598. A ThermalWall's own check
(``wall_verdicts``) stays a report.

Control volumes
---------------
* A two-port element. The stream enters at its upstream node's state and
  leaves at the downstream node's total pressure (a merge chamber's main
  inlet: at the chamber's main face, which is what it delivers), with the upstream
  composition and the enthalpy its walls added (the solver deposits wall
  heat on the downstream node, so the node's own temperature is a mix).
* A multi-port junction. Each port's stream at its port node's state.
* A chain of crossflow segments (a liner's backside duct). A segment's
  momentum spans half of each station, so one segment is not a control
  volume; the chain is. It is built from the segments' ``prev_seg`` /
  ``next_seg`` links; its stations are the nodes they declare
  (``from_kappa`` / ``to_kappa``), less any boundary. Streams crossing it:
  the chain's two ends, and every other element at a station (the holes), at
  the station's ``Pt`` (which a station's Pt = P closure makes its static
  pressure). Heat on a station node is taken as raising
  the node's mixed stream to the node's temperature. A hole whose jet brings
  axial momentum into the chain (``injection_angle_deg`` below 90, #813) is a
  member instead: its stream enters at its supply's total pressure. The
  chain's pressure recovery is the jet's, so closing the hole at the
  station's Pt would charge the chain with momentum it was given.

Not checked: elements without a solved flow, the legacy
``TeeJunctionElement``, and a stagnant control volume (no flow through it).
"""

from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

import combaero as cb

#: Below this, sigma_norm counts as negative. Chosen from a sweep of 238 element
#: instances across the suite's networks plus 50 random junctions: the most
#: negative best-case sigma_norm_high seen was 0 (a lossless connection sits
#: exactly at 0) and the smallest real loss was 1.5e-5 (an iterative junction),
#: so nothing came within three orders of magnitude of this bound from below.
#: 1e-8 is 100x tighter than the former 1e-6 -- it now catches a sub-ppm
#: manufactured total pressure -- while keeping a ~10x margin over the ~1e-9
#: solver round-off. Global, not per element: no element class showed a noise
#: floor that would need a looser per-element bound (#545).
SIGMA_TOL = 1e-8

#: Above this relative mass or energy closure of the network's control volume,
#: a converged solve is demoted to INCONSISTENT (#552). Over the suite's
#: 25 147 accepted solves the worst closures were 4.2e-8 (mass) and 1.1e-7
#: (energy); the roots #552 found accepted were open by 4e-3 to 8e-3 in mass
#: (an LM fallback) and 2e-5 in energy (a stale reversed state). Both
#: reproducers converge closed on today's solver; this is the guard that would
#: have refused them.
LEDGER_CLOSURE_TOL = 1e-6


@dataclass
class _Stream:
    m: float  # > 0
    Tt: float
    Pt: float
    Y: Any
    #: (Pt - P) / P of the node whose Pt the stream carries when that node is
    #: on the incompressible closure, else 0 (see ``_sigma``, #598).
    x_inc: float = 0.0


class _Probe:
    """The (Tt, Y) of a stream, as a turbomachine's outlet_stagnation reads them."""

    def __init__(self, stream: Any) -> None:
        self.Tt = stream.Tt
        self.Y = stream.Y


@dataclass
class _Volume:
    inflows: list[_Stream] = field(default_factory=list)
    outflows: list[_Stream] = field(default_factory=list)
    #: Entropy the heat brought in [W/K], one entry per heat input, as
    #: (Mach linear along the path, lowest, highest): inside an element only
    #: the end Mach numbers are known, so the path is bracketed.
    heat: list[tuple[float, float, float]] = field(default_factory=list)


def _X(Y: Any) -> np.ndarray:
    return np.asarray(cb.mass_to_mole(list(Y)), dtype=float)


def _sigma(cv: _Volume) -> dict[str, float]:
    """sigma [W/K], sigma_norm [-] and its band over the heat's Mach path, in C++."""

    def streams(sts: list[_Stream]) -> list[Any]:
        return [cb.MassStream(st.m, st.Tt, st.Pt, list(st.Y)) for st in sts]

    s_in, s_out = streams(cv.inflows), streams(cv.outflows)
    S_lin = sum(h[0] for h in cv.heat)
    S_low = sum(min(h) for h in cv.heat)
    S_high = sum(max(h) for h in cv.heat)
    central = cb.entropy_generation(s_in, s_out, S_lin)
    # Less heat entropy, more generation: the least favourable path is S_high.
    best = cb.entropy_generation(s_in, s_out, S_low) if cv.heat else central
    worst = cb.entropy_generation(s_in, s_out, S_high) if cv.heat else central
    # THE INCOMPRESSIBLE CLOSURE'S RESOLUTION (#598). sigma_norm carries
    # -m ln(Pt) per stream over sum_in m. A stream on the incompressible
    # closure has its energy as (Pt - p_ref) / rho_i, each port at its own
    # constant density, where the gas has R T ln(Pt / p_ref): the two differ
    # by x_i |P_i - p_ref| / P + x_i^2 / 2 in ln(Pt), x = (Pt - P) / P its
    # head. With the static spread inside a control volume taken as at most
    # its largest head, each stream is resolved to x_i (x_max + x_i / 2).
    # Order gamma^2 M^4 / 4: an element lossless in its own model (an ideal
    # tee's loss_factor = 0 split) sits on the bound, and measured on that
    # split at port Mach 0.033 the gap was -3.0e-7 against a bound of
    # 1.1e-6. Zero on a compressible node.
    m_in = sum(st.m for st in cv.inflows)
    every = cv.inflows + cv.outflows
    x_max = max((st.x_inc for st in every), default=0.0)
    res = (
        sum(st.m * st.x_inc * (x_max + 0.5 * st.x_inc) for st in every) / m_in
        if m_in > 0.0
        else 0.0
    )
    out = {
        "sigma": float(central.sigma),
        "sigma_norm": float(central.sigma_norm),
        "sigma_norm_low": float(worst.sigma_norm) - res,
        "sigma_norm_high": float(best.sigma_norm) + res,
    }
    if res:
        out["sigma_norm_resolution"] = res
    return out


class SecondLawReport:
    """Evaluates both checks on a solver's converged vector ``x``."""

    def __init__(self, solver: Any, x: np.ndarray, sol: dict[str, Any]) -> None:
        self.solver = solver
        self.net = solver.network
        self.x = x
        self.sol = sol
        self._states: dict[str, Any] = {}

    # -- states and flows -------------------------------------------------

    def _state(self, nid: str) -> Any:
        if nid not in self._states:
            self._states[nid] = self.solver._get_node_state(self.net.nodes[nid], self.x)
        return self._states[nid]

    def _flow(self, eid: str) -> float | None:
        idx = self.solver._unknown_indices.get(eid, [])
        if not idx or not self.solver.unknown_names[idx[0]].endswith(".m_dot"):
            return None
        return float(self.x[idx[0]])

    def _mach(self, nid: str) -> float:
        return float(self.sol.get(f"{nid}.mach", 0.0) or 0.0)

    def _wall_heat_into(self, eid: str) -> float:
        """Heat the thermal walls give element ``eid`` [W] (Q > 0 runs A -> B)."""
        q = 0.0
        for wid, wall in self.net.walls.items():
            Q = self.sol.get(f"{wid}.Q")
            if Q is None:
                continue
            if wall.element_a == eid:
                q -= float(Q)
            if wall.element_b == eid:
                q += float(Q)
        return q

    def _element_mach(self, eid: str, up: str, dn: str, forward: bool) -> tuple[float, float]:
        """The duct's Mach at the element's inlet and outlet, in flow order.

        An element that reports ``mach_in`` / ``mach_out`` knows its duct's;
        the end nodes may be plenums at Mach 0.
        """
        diag = (self.sol.get("__element_diag__") or {}).get(eid) or {}
        m_from = diag.get("mach_in")
        m_to = diag.get("mach_out")
        if m_from is None or m_to is None:
            return self._mach(up), self._mach(dn)
        a, b = abs(float(m_from)), abs(float(m_to))
        return (a, b) if forward else (b, a)

    def _heat_entropy(
        self, m: float, T0_in: float, T0_out: float, Y: Any, M_in: float, M_out: float
    ) -> tuple[float, float, float]:
        Y = list(Y)
        lin = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_in, M_out))
        if M_in == M_out:
            return lin, lin, lin
        a = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_in, M_in))
        b = float(cb.heat_entropy_flux(m, T0_in, T0_out, Y, M_out, M_out))
        return lin, min(a, b, lin), max(a, b, lin)

    def _incompressible_head(self, nid: str) -> float:
        """(Pt - P) / P at a node on the incompressible closure, else 0 (#598)."""
        if getattr(self.net.nodes[nid], "regime", None) != "incompressible":
            return 0.0
        st = self._state(nid)
        P, Pt = float(st.P), float(st.Pt)
        return max(Pt - P, 0.0) / P if P > 0.0 else 0.0

    def _stream_at(
        self, nid: str, m: float, Pt: float | None = None, pt_node: str | None = None
    ) -> _Stream:
        """The stream at ``nid``'s state; its Pt from ``Pt`` (taken at node
        ``pt_node``) when given."""
        st = self._state(nid)
        return _Stream(
            m,
            float(st.Tt),
            float(st.Pt if Pt is None else Pt),
            st.Y,
            self._incompressible_head(nid if Pt is None else (pt_node or nid)),
        )

    def _heated(self, inflow: _Stream, Q: float) -> float:
        """Total temperature after ``Q`` [W] is added to ``inflow``."""
        if not Q:
            return inflow.Tt
        X = _X(inflow.Y)
        h = float(cb.h_mass(inflow.Tt, X)) + Q / inflow.m
        # An enthalpy no positive T reaches (a wall taking 3.5 MJ/kg from a
        # stream holding 1.45): the inverse raises (#606).
        try:
            return float(cb.calc_T_from_h_mass(h, X, inflow.Tt))
        except RuntimeError as err:
            raise ValueError(
                f"the heat ({Q:.4g} W) takes more than the stream's enthalpy "
                f"({inflow.m:.4g} kg/s at {inflow.Tt:.1f} K)"
            ) from err

    def _outlet_Pt(self, eid: str, dn: str, m: float, forward: bool) -> float:
        """Total pressure the element delivers at its outlet node.

        A merge chamber's declared main inlet arrives at the chamber's main
        FACE, not at its mixed outlet state (#471): the side jets can lift the
        chamber's Pt above the main duct's own inlet Pt without the duct
        manufacturing anything, so closing the duct's control volume at the
        chamber Pt would charge it with the jets' momentum. The face is the
        one the solver hands the element (NetworkSolver._main_face).
        """
        node = self.net.nodes[dn]
        st = self._state(dn)
        if (
            forward
            and getattr(node, "main_inlet", None) == eid
            and getattr(node, "has_side_streams", lambda: False)()
        ):
            face_in = copy.copy(st)
            face_in.m_dot = m
            return float(node.main_face_state(face_in, m_main=m)[1])
        return float(st.Pt)

    # -- control volumes ---------------------------------------------------

    def _drive_pair(self, element: Any) -> Any:
        """A hole's drive pair (#1087, _HoleCdChain.drive_pair), or None
        where its streams are at its nodes' Pt."""
        chain = getattr(element, "_hole_cd", None)
        if chain is None:
            return None
        flows = self.solver._network_flows(element, self.x)
        return chain.drive_pair(self._state(element.from_node), self._state(element.to_node), flows)

    def _two_port(self, eid: str, element: Any) -> _Volume | str:
        m = self._flow(eid)
        if m is None:
            return "no solved flow"
        up, dn = (
            (element.from_node, element.to_node)
            if m >= 0.0
            else (
                element.to_node,
                element.from_node,
            )
        )
        m = abs(m)
        if m < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        forward = up == element.from_node
        inflow = self._stream_at(up, m)
        Q = self._wall_heat_into(eid)
        Tt_out = self._heated(inflow, Q)
        if getattr(element, "outlet_stagnation", None) is not None and forward:
            # A compressor or turbine (#893): shaft work enters as enthalpy and
            # carries no entropy, so the outlet is the machine's own Tt_out.
            # sigma >= 0 then holds for an efficiency <= 1.
            Tt_out = float(element.outlet_stagnation(_Probe(inflow))[0])
        outflow = _Stream(
            m,
            Tt_out,
            self._outlet_Pt(eid, dn, m, forward),
            inflow.Y,
            self._incompressible_head(dn),
        )
        pair = self._drive_pair(element)
        if pair is not None:
            # A hole between its drive pair (#1087): the stream leaves its
            # supply at the drive (a duct-fed hole's duct total) and its jet
            # discharges at the pressure it drives against (a flowing
            # chamber's static); reversed, it draws from that pressure and
            # discharges against the drive.
            inflow.Pt, outflow.Pt = (pair.p_up, pair.p_dn) if forward else (pair.p_dn, pair.p_rev)
        cv = _Volume([inflow], [outflow])
        if Q:
            M_in, M_out = self._element_mach(eid, up, dn, forward)
            cv.heat.append(self._heat_entropy(m, inflow.Tt, Tt_out, inflow.Y, M_in, M_out))
        return cv

    def _junction(self, element: Any) -> _Volume | str:
        cv = _Volume()
        for i, pid in enumerate(element.port_nodes):
            outer = element._port_element_ids[i]
            m = self._flow(outer) if outer else None
            if m is None:
                return f"port {pid} has no solved flow"
            m_out = float(element._port_signs[i]) * m  # > 0 leaves the junction
            if m_out == 0.0:
                continue
            st = self._stream_at(pid, abs(m_out))
            (cv.outflows if m_out > 0.0 else cv.inflows).append(st)
        if sum(st.m for st in cv.inflows) < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        return cv

    def _momentum_jets(self, members: list[str]) -> list[str]:
        """Holes at the chain's stations whose jets join it with axial
        momentum (#813): a non-normal ``injection_angle_deg``, flowing in."""
        stations = self._stations(members)
        out = []
        for eid, element in self.net.elements.items():
            if eid in members or element.to_node not in stations:
                continue
            angle = getattr(element, "injection_angle_deg", None)
            if angle is None or abs(np.cos(np.radians(float(angle)))) < 1e-12:
                continue
            m = self._flow(eid)
            if m is not None and m > 0.0:
                out.append(eid)
        return out

    def _chain(self, members: list[str], jets: list[str] | None = None) -> _Volume | str:
        """A crossflow chain: its stations are the nodes its members declare.
        ``jets`` enter at their supply's total pressure (``_momentum_jets``)."""
        els = self.net.elements
        stations = self._stations(members)
        jets = jets or []
        cv = _Volume()
        for eid, element in els.items():
            touches = [n for n in (element.from_node, element.to_node) if n in stations]
            if not touches:
                continue
            m = self._flow(eid)
            if m is None:
                return f"element {eid} at a station has no solved flow"
            if eid in members:
                # a chain end: the member's far node is outside the chain
                outside = [n for n in (element.from_node, element.to_node) if n not in stations]
                if not outside:
                    continue
                far = outside[0]
                into_chain = (m >= 0.0) == (element.from_node == far)
                st_node = far if into_chain else touches[0]
                stream = self._stream_at(st_node, abs(m), Pt=self._state(far).Pt, pt_node=far)
                if not into_chain:
                    # its wall heat lands on the far node: add it on the way out
                    Q = self._wall_heat_into(eid)
                    if Q:
                        T0_in = stream.Tt
                        stream.Tt = self._heated(stream, Q)
                        forward = element.to_node == far
                        M_in, M_out = self._element_mach(eid, touches[0], far, forward)
                        cv.heat.append(
                            self._heat_entropy(stream.m, T0_in, stream.Tt, stream.Y, M_in, M_out)
                        )
            else:
                # a hole (or any other element) at a station
                station = touches[0]
                into_chain = (m >= 0.0) == (element.to_node == station)
                src = element.from_node if m >= 0.0 else element.to_node
                at = src if eid in jets else station
                Pt = self._state(at).Pt
                pair = self._drive_pair(element) if at == element.from_node else None
                if pair is not None and pair.supply:
                    # A duct-fed hole at its drive pair (#1087): it draws at
                    # the duct's total, and a reversed jet discharges at
                    # p_rev (the total, or the station static where a
                    # converted static-referenced Cd recovers no head).
                    Pt = pair.p_up if m >= 0.0 else pair.p_rev
                stream = self._stream_at(
                    src if into_chain else station,
                    abs(m),
                    Pt=Pt,
                    pt_node=at,
                )
            if stream.m:
                (cv.inflows if into_chain else cv.outflows).append(stream)
        for nid in stations:
            node = self.net.nodes[nid]
            q = sum(float(eb.heat) for eb in getattr(node, "energy_boundaries", None) or [])
            q -= float(self.sol.get(f"{nid}.Q_withheld", 0.0))
            m_thru = sum(m for m, _, _ in self.solver._node_inflows(nid, self.x))
            if q and m_thru > 0.0:
                # the node's heat raised its mixed stream to the node's Tt
                st = self._state(nid)
                X = _X(st.Y)
                T0_out = float(st.Tt)
                h_before = float(cb.h_mass(T0_out, X)) - q / m_thru
                T0_in = float(cb.calc_T_from_h_mass(h_before, X, T0_out))
                M = self._mach(nid)
                cv.heat.append(self._heat_entropy(m_thru, T0_in, T0_out, st.Y, M, M))
        if sum(st.m for st in cv.inflows) < cb.MIXER_HEAT_MDOT_FLOOR:
            return "no flow"
        return cv

    def _stations(self, members: list[str]) -> set[str]:
        """Nodes a member declares a station (``from_kappa`` / ``to_kappa``).

        A boundary is never inside a control volume: an impingement array's
        last segment declares its exit (a PressureBoundary) a station whose
        onward flow is zero, and the stream into it leaves the chain.
        """
        out: set[str] = set()
        for eid in members:
            e = self.net.elements[eid]
            if e.from_kappa is not None:
                out.add(e.from_node)
            if e.to_kappa is not None:
                out.add(e.to_node)
        return {nid for nid in out if not self.net.nodes[nid].is_boundary}

    def _chains(self) -> list[list[str]]:
        from .components import CrossflowSegmentElement

        segs = {
            eid: e for eid, e in self.net.elements.items() if isinstance(e, CrossflowSegmentElement)
        }
        seen: set[str] = set()
        chains: list[list[str]] = []
        for eid in segs:
            if eid in seen:
                continue
            chain, todo = [], [eid]
            while todo:
                cur = todo.pop()
                if cur in seen or cur not in segs:
                    continue
                seen.add(cur)
                chain.append(cur)
                todo += [n for n in (segs[cur].prev_seg, segs[cur].next_seg) if n]
            chains.append(sorted(chain))
        return chains

    # -- reports -----------------------------------------------------------

    def element_verdicts(self) -> dict[str, dict[str, Any]]:
        from .components import MultiPortChamberBase, TeeJunctionElement

        out: dict[str, dict[str, Any]] = {}

        def record(eids: list[str], build: Any, group: str | None) -> None:
            # A REPORT must never take the solve down: a control volume that
            # cannot be evaluated (a converged but clamped state, heat beyond
            # what a stream holds) is not checked, and says why.
            try:
                cv = build()
                verdict_or_reason: dict[str, Any] | str = cv if isinstance(cv, str) else _sigma(cv)
            except (ValueError, RuntimeError) as exc:
                verdict_or_reason = f"cannot be evaluated: {exc}"
            if isinstance(verdict_or_reason, str):
                # The band keys stay, NaN: a result has the same keys whatever
                # its operating point (zero flow has no sigma per unit flow).
                verdict: dict[str, Any] = {
                    "verdict": "not_checked",
                    "reason": verdict_or_reason,
                    "sigma_norm": math.nan,
                    "sigma_norm_low": math.nan,
                    "sigma_norm_high": math.nan,
                }
            else:
                verdict = verdict_or_reason
                # Fail only where even the most favourable Mach path does.
                ok = verdict["sigma_norm_high"] >= -SIGMA_TOL
                verdict["verdict"] = "pass" if ok else "fail"
            if group:
                verdict["group"] = group
            for eid in eids:
                out[eid] = dict(verdict)

        chained: set[str] = set()
        for chain in self._chains():
            if not self._stations(chain):
                continue  # no side streams: each segment is a plain two-port
            jets = self._momentum_jets(chain)
            group = sorted(chain + jets)
            chained.update(group)
            record(
                group,
                lambda chain=chain, jets=jets: self._chain(chain, jets),
                group="+".join(group),
            )
        for eid, element in self.net.elements.items():
            if eid in chained:
                continue
            if isinstance(element, TeeJunctionElement):
                record([eid], lambda: "legacy TeeJunctionElement", None)
            elif isinstance(element, MultiPortChamberBase):
                record([eid], lambda element=element: self._junction(element), None)
            else:
                record([eid], lambda eid=eid, element=element: self._two_port(eid, element), None)
        return out

    def _wall_side_entropy(self, eid: str, Q_in: float) -> tuple[float, float, float]:
        """Entropy [W/K] the heat ``Q_in`` brings into element ``eid``'s stream,
        integrated along the stream's own temperature path (lin, low, high)."""
        element = self.net.elements.get(eid)
        if element is None:
            raise ValueError(f"'{eid}' is a node, not a flowing element")
        n_walls = sum(1 for w in self.net.walls.values() if eid in w.heated_endpoints)
        if n_walls > 1:
            raise ValueError(f"'{eid}' carries {n_walls} walls; their heat paths are shared")
        m = self._flow(eid)
        if m is None or abs(m) < cb.MIXER_HEAT_MDOT_FLOOR:
            raise ValueError(f"'{eid}' has no solved flow")
        up, dn = (
            (element.from_node, element.to_node)
            if m >= 0.0
            else (element.to_node, element.from_node)
        )
        inflow = self._stream_at(up, abs(m))
        Tt_out = self._heated(inflow, Q_in)
        M_in, M_out = self._element_mach(eid, up, dn, up == element.from_node)
        return self._heat_entropy(abs(m), inflow.Tt, Tt_out, inflow.Y, M_in, M_out)

    def wall_verdicts(self) -> dict[str, dict[str, Any]]:
        """Each ThermalWall as its own control volume (#551).

        The heat leaves side A's stream and enters side B's, each at that
        stream's own temperature along its path, so the wall generates
        ``sigma_wall = S_A + S_B`` [W/K] with ``S = integral dQ / T`` (negative
        on the side giving heat). It is negative when heat runs from the
        stream that is colder on average to the hotter one -- a wall model
        with no capacity limit can do that (a hot stream cooled below the cold
        stream's inlet). Evaluating it at the temperatures the wall model sees
        instead would be vacuous: ``Q = U (T_A - T_B)`` makes ``Q (1/T_B -
        1/T_A)`` non-negative by construction.

        A necessary condition, not a sufficient one: the path means can stay
        ordered while the end temperatures cross. ``sigma_wall_norm`` is
        ``sigma_wall`` over the mean ``|S|`` of the two sides. A REPORT only.
        """
        out: dict[str, dict[str, Any]] = {}
        for wid, wall in self.net.walls.items():
            Q = self.sol.get(f"{wid}.Q")
            if Q is None:
                out[wid] = {"verdict": "not_checked", "reason": "no heat solved for this wall"}
                continue
            if wall.adiabatic:
                out[wid] = {"verdict": "pass", "sigma_wall": 0.0, "sigma_wall_norm": 0.0}
                continue
            try:
                # Q > 0 runs A -> B.
                S_a = self._wall_side_entropy(wall.element_a, -float(Q))
                if wall.heat_loss is not None:
                    # A reservoir at T_ext takes the heat (#970): Q / T_ext.
                    s_ext = float(Q) / float(wall.heat_loss.T_ext)
                    S_b = (s_ext, s_ext, s_ext)
                else:
                    S_b = self._wall_side_entropy(wall.element_b, float(Q))
            except (ValueError, RuntimeError) as exc:
                out[wid] = {"verdict": "not_checked", "reason": f"cannot be evaluated: {exc}"}
                continue
            scale = 0.5 * (abs(S_a[0]) + abs(S_b[0]))
            if scale <= 0.0:
                out[wid] = {"verdict": "pass", "sigma_wall": 0.0, "sigma_wall_norm": 0.0}
                continue
            sigma = S_a[0] + S_b[0]
            # The most favourable paths: the largest entropy on both sides.
            sigma_high = S_a[2] + S_b[2]
            out[wid] = {
                "verdict": "pass" if sigma_high / scale >= -SIGMA_TOL else "fail",
                "sigma_wall": sigma,
                "sigma_wall_norm": sigma / scale,
                "sigma_wall_norm_low": (S_a[1] + S_b[1]) / scale,
                "sigma_wall_norm_high": sigma_high / scale,
            }
        return out

    def ledger(self) -> dict[str, Any]:
        """Mass and energy over the network's control volume."""
        rows: list[dict[str, Any]] = []

        def add(node: str, via: str, m_in: float, Tt: float, Y: Any) -> None:
            H = m_in * float(cb.h_mass(Tt, _X(Y)))
            rows.append({"node": node, "via": via, "m_dot_in": m_in, "H_in": H})

        for nid, node in self.net.nodes.items():
            if not node.is_flow_boundary:
                continue
            ups = self.net.get_upstream_elements(nid)
            downs = self.net.get_downstream_elements(nid)
            if ups and downs:
                # an injection: only its own stream crosses here
                Y_inj = node.Y if node.Y is not None else self.solver._default_Y
                add(nid, nid, float(node.m_dot), float(node.Tt), Y_inj)
                continue
            for element in {e.id: e for e in [*ups, *downs]}.values():
                ind = self.solver._unknown_indices.get(element.id, [])
                if not ind:
                    continue
                # signed source -> sink, as the element itself reports it here
                m = float(element.flow_at_node(nid, self.x, ind))
                at_source = nid in element.all_source_nodes()
                if at_source == (m >= 0.0):
                    st = self._state(nid)
                    add(nid, element.id, abs(m), float(st.Tt), st.Y)
                    continue
                # leaving, at the mixed specific enthalpy of what feeds it
                far = element.all_sink_nodes() if at_source else element.all_source_nodes()
                m_far = H_far = 0.0
                for src in far:
                    if src == nid:
                        continue
                    m_src = abs(float(element.flow_at_node(src, self.x, ind)))
                    st = self._state(src)
                    m_far += m_src
                    H_far += m_src * float(cb.h_mass(float(st.Tt), _X(st.Y)))
                H = abs(m) * H_far / m_far if m_far > 0.0 else 0.0
                if getattr(element, "outlet_stagnation", None) is not None and m_far > 0.0:
                    # Leaving through a compressor or turbine (#893): at its
                    # outlet Tt, not at the enthalpy of what feeds it.
                    st = self._state(far[0])
                    Tt_out = float(element.outlet_stagnation(st)[0])
                    H = abs(m) * float(cb.h_mass(Tt_out, _X(st.Y)))
                rows.append({"node": nid, "via": element.id, "m_dot_in": -abs(m), "H_in": -H})

        q_given = 0.0
        for nid, node in self.net.nodes.items():
            for eb in getattr(node, "energy_boundaries", None) or []:
                if eb.id != f"_wall_{nid}":
                    q_given += float(eb.heat)  # the SOLVED heat of a Tt_target (#792)
            q_given += float(self.sol.get(f"{nid}.Q_fraction", 0.0))
        # Shaft work into the flow (#893): m dh of every compressor (> 0) and
        # turbine (< 0), an energy input like a given heat, with no entropy.
        w_shaft = 0.0
        for eid, element in self.net.elements.items():
            if getattr(element, "outlet_stagnation", None) is None:
                continue
            m = self._flow(eid)
            dh = self.sol.get(f"{eid}.dh")
            if m is not None and dh is not None and math.isfinite(float(dh)):
                w_shaft += abs(float(m)) * float(dh)
        q_out = sum(float(v) for k, v in self.sol.items() if k.endswith(".Q_wall_out"))
        q_withheld = sum(float(v) for k, v in self.sol.items() if k.endswith(".Q_withheld"))
        # Heat a wall's HeatLoss side takes out of the network (#970).
        q_loss = sum(float(self.sol.get(f"{wid}.Q_loss", 0.0)) for wid in self.net.walls)

        m_in = sum(r["m_dot_in"] for r in rows if r["m_dot_in"] > 0.0)
        m_out = -sum(r["m_dot_in"] for r in rows if r["m_dot_in"] < 0.0)
        H_in = sum(r["H_in"] for r in rows if r["m_dot_in"] > 0.0)
        H_out = -sum(r["H_in"] for r in rows if r["m_dot_in"] < 0.0)
        energy = H_in + q_given + w_shaft - H_out - q_out - q_withheld - q_loss
        h_scale = max(abs(H_in), abs(H_out), sum(abs(r["H_in"]) for r in rows), 1.0)
        return {
            "streams": rows,
            "m_dot_in": m_in,
            "m_dot_out": m_out,
            "mass_closure": m_in - m_out,
            "mass_closure_rel": (m_in - m_out) / max(m_in, m_out, cb.MIXER_HEAT_MDOT_FLOOR),
            "H_in": H_in,
            "H_out": H_out,
            "Q_given": q_given,
            "W_shaft": w_shaft,
            "Q_wall_out": q_out,
            "Q_withheld": q_withheld,
            "Q_loss": q_loss,
            "energy_closure": energy,
            "energy_closure_rel": energy / h_scale,
        }
