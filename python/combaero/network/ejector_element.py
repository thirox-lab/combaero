"""EjectorElement: 3-port supersonic ejector on the momentum-CV junction
topology (``MultiPortChamberBase``).

Ports: ``primary_node`` (motive, high-pressure inlet), ``secondary_node``
(suction/entrained inlet), ``outlet_node``. Reuses the validated closed-form
physics in ``combaero.network._ejector_huang1999`` (Huang 1999's entrainment
ratio, Eqs. 1-8; the (P_c*, P_b0) anchors by ``anchor_closure``: Huang's
printed P_c* chain with a sudden-expansion dead head by default, or Kracik &
Dvorak 2016's lossless mixing, their Eqs. 7-13) -- see that module's docstring and
``validation/ejector/README.md`` for the full derivation, paper references,
and the accuracy comparison against three alternative closures that were
also evaluated.

Operating regimes. This element resolves all three physical regimes in one
Newton-solvable residual system (design + provenance in
``validation/ejector/OPERATING_REGIMES_DESIGN.md``):

  * CRITICAL (double-choked): primary choked, entrained flow critical, so
    ``omega`` is back-pressure-independent. The outlet floats freely below
    the ejector's critical back pressure ``P_c*``.
  * SUBCRITICAL droop: primary still choked but the outlet is pushed above
    ``P_c*``; the outlet pressure sets the operating point.
  * UNCHOKED jet pump: the forced primary mass flow is below the choke
    threshold ``mdot_choked(P_back)``, so the primary Pt floats to just above
    the back pressure (forward flow, ``dp >= 0``) and the device acts as a
    subsonic constant-area jet pump (Keenan-Neumann-Lustwerk mixing via
    ``ejector_jetpump_discharge``).

The regime is selected smoothly (no branching discontinuity) by two
smootherstep weights: ``s_choke = smootherstep(mp / mdot_choked(Pt_p))``
(0 = unchoked jet pump, 1 = choked) and ``s_sub`` (0 = outlet below ``P_c*``,
critical; 1 = outlet above it, subcritical). See the residual design below.

Working fluid: combaero's real-fluid EOS covers combustion/air-relevant
species only (N2, O2, Ar, CO2, H2O, light hydrocarbons, CO, H2, NH3) -- no
halocarbon refrigerants. gamma and R are evaluated LIVE from that EOS at
each residual call, at the entrained-flow choking plane (same convention
validated in the reference model, there baked from CoolProp for R141b; here
computed on the fly via ``cb.cp_mass``/``cb.cv_mass``), using the secondary
stream's own composition -- Huang's model assumes a single working fluid, so
this is the representative gamma/R for the whole closed form.

Residual design. Ports are ``[primary, secondary, outlet]``, so
``_port_signs = [-1, -1, +1]`` and the base class's ``unknowns()``/
``n_equations()`` (one owned unknown ``{id}.P_jct``, N+1 = 4 rows) are
inherited unchanged -- but the owned unknown is REPURPOSED as the physical
mixing-plane static pressure ``P_py`` (the pressure common to the primary jet
and the entrained stream at the mixing-chamber entrance), and the row content
is ejector physics instead of impulse-CV physics. With ``s = s_choke``:

    R0 = mp - cd_nozzle(Pt_p, P_py, A_t, A_e)      (primary C-D nozzle)
    R1 = ms - [s*omega_ch*mp + (1-s)*ms_sec]       (blended entrainment)
    R2 = mdot_out - mp - ms                        (mass conservation)
    R3 = [w_pin*Pt_outlet + (1-w_pin)*P_py] - recovery*discharge

``omega_ch`` is the choked entrainment: omega_crit up to P_c*, then the
Tier-1 subcritical droop, linear in the outlet Pt down to 0 at the dead-head
P_b0 (#619; the sec 8.1 closure, which the element did not use before).

``P_py`` is a genuine Newton unknown (not derived): it is the quantity the two
nozzle relations R0/R1 share, so making it free keeps BOTH live and lets the
downstream back pressure select the operating regime. ``R0``'s C-D nozzle
smooth-mins the subsonic-flux and sonic-cap branches, so it is choked or
unchoked automatically as ``P_py`` moves. It is strictly monotone in ``P_py`` over the
whole line (#619): the exit Mach is held at 1 below the critical ratio (the
supersonic-exit flux used to dip under the cap at deep expansion), and above
the supply the nozzle passes REVERSE flow, so zero flow pins ``P_py`` to the
supply instead of leaving it free (zero flow at any ``P_py`` above both
supplies used to close every row: a spurious root at ~1e13 Pa). ``discharge = recovery * (s*P_c* +
(1-s)*P03_jetpump)`` blends the critical and jet-pump mixing closures; each is
evaluated ONLY where its weight is nonzero (the critical closures return NaN
for an unchoked primary, and ``0 * NaN = NaN`` would otherwise poison the
Jacobian -- safe because smootherstep is flat at ``s in {0, 1}``).

``R3`` closes the system differently per regime via
``w_pin = 1 - s*(1 - s_b0)`` (the weight on the outlet-pin form; ``s_b0`` is
the smootherstep of ``Pt_outlet / (recovery*P_b0)`` over [0.98, 1.02]):

  * jet pump (s=0): ``w_pin = 1`` -> ``Pt_outlet = recovery*discharge``. The
    outlet is what the ejector sets; the downstream boundary then forces
    ``P_py`` through this pin, and R0/R1 set the primary Pt and entrainment.
  * choked, outlet below the dead-head ``P_b0`` (critical AND subcritical,
    s=1, s_b0=0): ``w_pin = 0`` -> ``P_py = recovery*P_c*`` (a diagnostic) and
    the outlet floats FREE, set by the downstream network, the entrainment
    following it along the droop above ``P_c*``.
  * choked, outlet above ``P_b0`` (s=1, s_b0=1): ``w_pin = 1`` ->
    ``Pt_outlet = recovery*P_b0``. This pin excludes the spurious double-choked
    root (primary Pt below the back pressure), which can only occur above
    ``P_b0`` since ``P_b0 < Pt_primary``. Until #619 it switched in at ``P_c*``
    instead, which also pinned every subcritical outlet below its boundary:
    a choked ejector between ``P_c*`` and ``P_b0`` had no root.

A collapsed droop window (``P_b0`` within 0.1-1% of ``P_c*``, the unchoked
``M_py < 1`` corner) fades back to the plateau-and-pin-at-``P_c*`` form
continuously, since ``P_b0`` depends on Newton unknowns. Beyond ``P_b0`` (the
ejector cannot hold the back pressure: backflow) is outside the model.

``verify_solution_consistent`` now returns True unconditionally: every regime
is modeled, so there is no longer a subcritical design point to demote.

Assembly: the whole 4-row residual system AND its 9-column analytic
Jacobian are built in C++ -- a single ``_solver_tools.
ejector_element_residuals_and_jacobian`` call (``include/ejector.h`` /
``src/ejector.cpp``) that chains the scalar closures through a forward-mode
``DualN<9>`` across the regime blend, matching the whole-element (f, J)
practice of the base ``MultiPortChamberBase`` and ``TeeJunctionElement``.
``residuals()`` here is only the relabeling shim: physical-flow signs and
mapping the nine Jacobian seeds to solver column names.
(``combaero.network._ejector_huang1999`` remains the Python validation
reference the C++ is kept 1:1 with.) The C++ Jacobian is exact w.r.t. the
nine Newton unknowns holding ``gamma`` and ``r_gas`` FROZEN. Here
``gamma = choke_plane_gamma(secondary.Tt, ...)`` is a weak function of a
Newton unknown and ``r_gas`` depends only on composition (not a solved
unknown), so the element Jacobian is exact in the explicit (Pt, Tt)
dependence and drops only the tiny implicit ``d gamma / d Tt`` term -- a
standard frozen-property Jacobian. ``gamma`` is recomputed fresh every
residual call, so the converged root is still exact; only the Newton step
uses a slightly approximate slope. This matches the reference data's own
convention (``generate_reference.py`` treats gamma as a separate Jacobian
input, so its d/dt_e central-difference targets also hold gamma fixed).
``choke_plane_gamma`` therefore stays in Python.
"""

from __future__ import annotations

import contextlib
import math
from typing import Any

import combaero as cb
from combaero import _solver_tools
from combaero.network._ejector_huang1999 import ANCHOR_CLOSURES, ETA_P, ETA_S, EjectorGeometry
from combaero.network.components import (
    DEFAULT_REGIME,
    NEIGHBOUR_AREA_RTOL,
    CompressibilityLiteral,
    MomentumChamberNode,
    MultiPortChamberBase,
    NetworkMixtureState,
    _attached_pipe,
    _element_area_facing,
)

#: Mixing-section area A_3 [m^2] used when it is neither set nor inferable
#: from the outlet pipe (#559). Matches the GUI's historical default.
DEFAULT_EJECTOR_MIXING_AREA = 8.0e-4

#: Reverse flow at a port, relative to the largest port flow, beyond which a
#: converged root is rejected as outside the model (#619). Exactly at shut-off
#: the suction flow is zero to solver precision, about 1e-10 relative.
EJECTOR_REVERSE_FLOW_TOL = 1e-6


def _real_gamma(t: float, x: Any) -> float:
    """Ideal-gas-like gamma = cp/cv from combaero's own EOS at temperature
    t and composition x (mole fractions). cb.cp_mass/cv_mass are functions
    of (T, X) only, matching the pressure-independent NASA-polynomial-style
    cp0(T) convention the reference model's CoolProp-baked gamma also uses.
    """
    cp = cb.cp_mass(t, x)
    cv = cb.cv_mass(t, x)
    return cp / cv if cv > 0.0 else 1.4


def choke_plane_gamma(t_e: float, x: Any, *, iterations: int = 12) -> float:
    """gamma at the entrained-flow choking plane y-y (T_sy = T_e/((g+1)/2)),
    a self-consistent fixed point since gamma appears on both sides through
    T_sy -- same convention validated in the reference model (see
    ``validation/ejector/README.md``'s "gamma provenance" section), here
    evaluated live from combaero's EOS instead of baked from CoolProp.
    """
    gamma = _real_gamma(t_e, x)
    for _ in range(iterations):
        t_sy = t_e / ((gamma + 1.0) / 2.0)
        gamma = _real_gamma(t_sy, x)
    return gamma


class EjectorElement(MultiPortChamberBase):
    """Supersonic ejector across operating regimes: 2 inlets (primary, secondary),
    1 outlet. Critical (double-choked), subcritical, and unchoked-primary
    (subsonic jet-pump) in one C1 residual system."""

    # The ejector's port flows follow choking and entrainment, not a Bernoulli
    # share of the imposed pressure differences, and it carries its own
    # physics-based warm start for the port pressures and P_jct. The solver's
    # junction split seed would overwrite that with a split the ejector does
    # not obey, and the cold reference network then fails to converge.
    seeds_ports_by_pressure_split: bool = False

    # The incompressible validity report (#588) judges the primary jet: an
    # ejector's ports can be slow while its jet is not.
    validity_mach_keys = ("mach_jet",)

    # to_dict serialises the given mixing area, not the inherited one (#1099).
    _init_kwarg_attrs = {"mixing_area": lambda e: None if e._auto_mixing else e.mixing_area}

    def __init__(
        self,
        id: str,
        primary_node: str,
        secondary_node: str,
        outlet_node: str,
        throat_area: float,
        nozzle_exit_area: float,
        mixing_area: float | None = None,
        recovery_efficiency: float = 1.0,
        regime: CompressibilityLiteral = DEFAULT_REGIME,
        k_nozzle: float | None = None,
        k_suction: float | None = None,
        k_throat: float = 0.0,
        k_diffuser: float = 0.0,
        anchor_closure: str = "huang1999",
    ):
        """Build an ejector.

        Args:
            id: Element id.
            primary_node: Motive/high-pressure inlet port-MCN id.
            secondary_node: Suction/entrained inlet port-MCN id.
            outlet_node: Outlet port-MCN id.
            throat_area: Primary nozzle throat area A_t [m^2].
            nozzle_exit_area: Primary nozzle exit area A_p1 [m^2]. Must
                exceed throat_area (the nozzle must diverge to expand
                supersonically).
            mixing_area: Constant-area mixing-section area A_3 [m^2]. Must
                exceed nozzle_exit_area (room for the entrained flow). None =
                inherit it from the outlet pipe (the mixing duct and the pipe
                downstream are the same duct), else DEFAULT_EJECTOR_MIXING_AREA.
                The inlet faces never inherit: the nozzle exit and the suction
                annulus are internal geometry, not the feed pipes' areas.
            recovery_efficiency: Multiplies the lossless mixed stagnation
                pressure to give the critical back pressure P_c*, in (0, 1]:
                above 1 the ejector would gain pressure beyond lossless
                mixing. Defaults to 1.0 (no artificial loss); see this module's and
                ``_ejector_huang1999``'s docstrings for why that default is
                not fitted to any reference dataset. The incompressible
                closure applies it the same way, to the outlet's stagnation
                pressure.
            regime: ``"compressible"`` (the closures above) or
                ``"incompressible"``: the low-Mach limit, Sanger's 1-D jet
                pump (NASA TN D-4445, 1968), each stream at its own constant
                density, no choking and no critical regime (#588). Default
                ``DEFAULT_REGIME``, like every element.
            k_nozzle, k_suction, k_throat, k_diffuser: Sanger's component
                loss coefficients (K_p, K_s, K_t, K_d), on the nozzle-exit,
                suction-annulus and mixing-duct dynamic heads; incompressible
                regime only. ``None`` (the nozzle and suction defaults) takes
                the compressible closure's own nozzle efficiencies, Huang's
                ETA_P = 0.95 and ETA_S = 0.85, as K = 1/eta - 1 (an efficiency
                on kinetic energy IS a loss coefficient: Pt - p = rho V^2 /
                (2 eta)), so both regimes agree at low Mach out of the box.
                The throat and diffuser default to 0, as the compressible
                mixing has no friction and ``recovery_efficiency`` covers the
                outlet. Sanger calibrated 0.008-0.036 / 0.09-0.14 /
                0.076-0.0985 / 0.079-0.102 on his water jet pumps.
            anchor_closure: The (P_c*, P_b0) pair that ends the critical
                plateau and the subcritical droop; compressible regime only.
                ``"huang1999"`` (default): Huang 1999's printed chain (phi_m
                mixing, normal shock, isentropic diffuser) and the
                sudden-expansion dead head.
                ``"kracik_dvorak"``: Kracik & Dvorak's lossless mixing, at the
                critical entrainment and at zero entrainment (the only
                closure up to 0.6.0). Cross-source scores of both, and the
                printed chain's fidelity gap, are in validation/ejector/README.md.
        """
        # Validate geometry BEFORE super().__init__ so the port-face areas
        # passed below are guaranteed positive.
        # Written "not (x > y)" so NaN fails them too (#769); inf is refused
        # explicitly.
        if not (throat_area > 0.0) or not math.isfinite(throat_area):
            raise ValueError(f"EjectorElement '{id}': throat_area must be positive and finite.")
        if not (nozzle_exit_area > throat_area) or not math.isfinite(nozzle_exit_area):
            raise ValueError(
                f"EjectorElement '{id}': nozzle_exit_area ({nozzle_exit_area}) must "
                f"exceed throat_area ({throat_area}) -- the primary nozzle must "
                f"diverge to expand supersonically."
            )
        if mixing_area is not None and (
            not (mixing_area > nozzle_exit_area) or not math.isfinite(mixing_area)
        ):
            raise ValueError(
                f"EjectorElement '{id}': mixing_area ({mixing_area}) must exceed "
                f"nozzle_exit_area ({nozzle_exit_area}) -- the constant-area section "
                f"must leave room for the entrained flow."
            )
        if anchor_closure not in ANCHOR_CLOSURES:
            raise ValueError(
                f"EjectorElement '{id}': anchor_closure must be one of "
                f"{ANCHOR_CLOSURES}, got {anchor_closure!r}."
            )
        if not (0.0 < recovery_efficiency <= 1.0):
            raise ValueError(
                f"EjectorElement '{id}': recovery_efficiency must be positive and at most 1 "
                f"(above 1 is a pressure gain beyond lossless mixing), got {recovery_efficiency}."
            )
        if k_nozzle is None:
            k_nozzle = 1.0 / ETA_P - 1.0
        if k_suction is None:
            k_suction = 1.0 / ETA_S - 1.0
        for name, k in (
            ("k_nozzle", k_nozzle),
            ("k_suction", k_suction),
            ("k_throat", k_throat),
            ("k_diffuser", k_diffuser),
        ):
            if not (k >= 0.0) or not math.isfinite(k):
                raise ValueError(f"EjectorElement '{id}': {name} must be finite and >= 0, got {k}.")
        # Port-face areas from the ejector's OWN geometry (ordered
        # inlet_nodes + outlet_nodes = [primary, secondary, outlet]). Unlike a
        # generic tee, the ejector always knows these, so it supplies them
        # rather than relying on the base class inferring them from a
        # connecting ChannelElement -- an ejector can be fed by a bare
        # boundary + connection with no channel/diameter to infer from.
        # These set the port-MCN dynamic-head (Pt = P + 0.5*rho*v^2) only;
        # the ejector's own residuals use throat/nozzle_exit/mixing directly.
        #   primary  -> nozzle exit A_p1 (motive jet plane)
        #   secondary-> annular suction A_3 - A_p1
        #   outlet   -> mixing/diffuser duct A_3
        # An inherited A_3 is filled in by resolve_topology, which sets all
        # three (the base class's own inheritance would give the suction port
        # the suction pipe's area, which is not the annulus).
        port_areas: list[float | None] = (
            [
                float(nozzle_exit_area),
                float(mixing_area - nozzle_exit_area),
                float(mixing_area),
            ]
            if mixing_area is not None
            else [float(nozzle_exit_area), None, None]
        )
        super().__init__(
            id,
            inlet_nodes=[primary_node, secondary_node],
            outlet_nodes=[outlet_node],
            port_areas=port_areas,
        )
        self.primary_node = primary_node
        self.secondary_node = secondary_node
        self.outlet_node = outlet_node
        self.throat_area = float(throat_area)
        self.nozzle_exit_area = float(nozzle_exit_area)
        self._auto_mixing = mixing_area is None
        self.mixing_area: float | None = float(mixing_area) if mixing_area is not None else None
        #: Where mixing_area came from (#559): "user", "channel '<id>'" /
        #: "<Element> '<id>'" (the outlet pipe), or "default".
        self._mixing_source = "user" if mixing_area is not None else ""
        self.recovery_efficiency = float(recovery_efficiency)
        self.anchor_closure = anchor_closure
        self.regime = regime
        self.k_nozzle = float(k_nozzle)
        self.k_suction = float(k_suction)
        self.k_throat = float(k_throat)
        self.k_diffuser = float(k_diffuser)
        #: Warm-start values this element wrote (see seed_warm_start), keyed
        #: by (holder id, unknown); a guess the user gave is never in here.
        #: FlowNetwork.to_dict leaves these out of the saved guesses (#1101).
        self._auto_seeds: dict[tuple[str, str], float] = {}

    def area_facing(self, node_id: str) -> float | None:
        """The flow area this ejector presents at ``node_id``: the mixing duct
        A_3 at the outlet, so an outlet pipe can inherit it. None at the inlets
        (internal nozzle exit / suction annulus, not a feed pipe's area)."""
        if node_id == self.outlet_node and self.mixing_area is not None:
            return self.mixing_area
        return None

    def resolve_topology(self, graph: Any) -> None:
        if self._auto_mixing:
            area, source = self._outlet_pipe_area(graph)
            if area is not None:
                self.mixing_area, self._mixing_source = area, source
            elif self.mixing_area is None:
                self.mixing_area = DEFAULT_EJECTOR_MIXING_AREA
                self._mixing_source = "default"
            if not (self.mixing_area > self.nozzle_exit_area):
                raise ValueError(
                    f"EjectorElement '{self.id}': the mixing area inherited from "
                    f"the outlet pipe ({self.mixing_area:.4g} m^2, {self._mixing_source}) "
                    f"must exceed nozzle_exit_area ({self.nozzle_exit_area:.4g} m^2). "
                    f"Set mixing_area explicitly or widen the outlet pipe."
                )
            a3, ap1 = self.mixing_area, self.nozzle_exit_area
            self.port_areas = [ap1, a3 - ap1, a3]
            self._port_area_sources = ["user", self._mixing_source, self._mixing_source]
        else:
            self._check_mixing_against_outlet_pipe(graph)
        self._size_outlet_port_node(graph)
        super().resolve_topology(graph)
        self.seed_warm_start(graph)

    def seed_warm_start(self, graph: Any) -> None:
        """Seed the cold start from the ejector's own closed-form physics.

        The stiff compressible solve needs it: from the generic cold start
        (P_jct at the highest boundary pressure, i.e. the primary nozzle's
        zero-flow point) 2 of 22 stress draws converge, against 9 of 22
        seeded (#619). Each port is traced to its boundary through the
        network; the port chambers' P/Pt, ``P_jct`` (= P_c* critical, just
        below the back pressure in the jet-pump case) and the port-chain
        flows are seeded where the user gave no guess. Runs on every
        ``resolve_topology``: it first withdraws the seeds it wrote last
        time, so a re-solve with new boundary values is seeded afresh, and
        seeds nothing in the incompressible regime, whose ``P_jct`` means
        something else (the solver's incompressible proxy flips it there).
        Best-effort: a closed form that refuses leaves that seed unset.
        Moved from the GUI backend, so API-built networks get it too.
        """
        for (holder_id, key), value in self._auto_seeds.items():
            holder = graph.nodes.get(holder_id) or graph.elements.get(holder_id)
            guess = dict(getattr(holder, "initial_guess", {}) or {}) if holder else {}
            if guess.get(key) == value:
                del guess[key]
                holder.initial_guess = guess
        self._auto_seeds = {}
        if self.regime != "compressible":
            return
        with contextlib.suppress(Exception):  # a seed, never a failure
            self._seed_from_closed_forms(graph)

    def _seed(self, holder: Any, key: str, value: float) -> None:
        guess = dict(getattr(holder, "initial_guess", {}) or {})
        if key in guess:
            return
        guess[key] = float(value)
        holder.initial_guess = guess
        self._auto_seeds[(holder.id, key)] = float(value)

    def compressible_proxy_overrides(self) -> dict[str, Any]:
        """Attributes the solver sets while this incompressible element runs
        compressible in its warm-start proxy (``_compressible_seed``).

        The proxy only seeds the incompressible solve, which uses no anchors,
        so it takes the lossless pair: its wider droop window keeps a forward
        root where the Huang dead head has none. On the stress family's 191
        incompressible ejectors that proxy converges 24 times against 1, and
        the family's time falls 51 -> 40 s with the same 190 solved.
        """
        return {"anchor_closure": "kracik_dvorak"}

    def no_root_hint(self, graph: Any) -> str:
        """Why a failed solve may have no root here: the outlet boundary above
        the dead-head back pressure P_b0 the supply boundaries give. Above it a
        real ejector backflows (the suction turns outward), a regime this model
        does not cover (#619), so a network that needs it has no forward-flow
        root. Evaluated at the boundary states, so it needs no solved point.
        "" when the ports are not all pressure boundaries, or for the
        incompressible jet pump (no dead head)."""
        from combaero.network.components import PressureBoundary

        if self.regime != "compressible":
            return ""
        pb, sb, ob = (
            self._port_boundary(graph, n)
            for n in (self.primary_node, self.secondary_node, self.outlet_node)
        )
        if not all(isinstance(b, PressureBoundary) for b in (pb, sb, ob)):
            return ""
        a3 = self.mixing_area
        if a3 is None:
            a3 = self._outlet_pipe_area(graph)[0] or DEFAULT_EJECTOR_MIXING_AREA
        try:
            x = list(cb.mass_to_mole(list(sb.Y)))
            gamma = choke_plane_gamma(float(sb.Tt), x)
            r_gas = float(cb.specific_gas_constant(x))
            closure = (
                _solver_tools.ejector_sudden_expansion_dead_head_and_jacobian
                if self.anchor_closure == "huang1999"
                else _solver_tools.ejector_dead_head_back_pressure_and_jacobian
            )
            p_b0 = closure(
                float(pb.Pt),
                float(pb.Tt),
                float(sb.Pt),
                float(sb.Tt),
                self.area_ratio_nozzle,
                float(a3) / self.throat_area,
                gamma,
                r_gas,
                self.recovery_efficiency,
            ).value.p_c_pa
        except Exception:  # noqa: BLE001 -- a diagnostic never breaks a solve
            return ""
        if not float(ob.Pt) > p_b0:
            return ""
        return (
            f"EjectorElement '{self.id}': the outlet boundary ({float(ob.Pt) / 1e5:.4g} bar) is "
            f"above the dead-head back pressure P_b0 = {p_b0 / 1e5:.4g} bar that its supplies "
            f"give (anchor_closure='{self.anchor_closure}'). A real ejector backflows there, "
            "the suction turning outward, which this model does not cover: lower the outlet "
            "pressure or raise the primary supply."
        )

    def _port_boundary(self, graph: Any, port_id: str, max_hops: int = 8) -> Any:
        """The first flow boundary reached walking outward from a port, never
        through this ejector; None if there is none within ``max_hops``."""
        from combaero.network.components import MassFlowBoundary, PressureBoundary

        seen, frontier = {port_id}, [port_id]
        for _ in range(max_hops):
            nxt: list[str] = []
            for nid in frontier:
                for elem in graph.get_upstream_elements(nid) + graph.get_downstream_elements(nid):
                    if elem is self:
                        continue
                    for other in (getattr(elem, "from_node", None), getattr(elem, "to_node", None)):
                        if other is None or other in seen:
                            continue
                        seen.add(other)
                        node = graph.nodes.get(other)
                        if isinstance(node, (PressureBoundary, MassFlowBoundary)):
                            return node
                        nxt.append(other)
            frontier = nxt
        return None

    def _seed_port_chain_flows(self, graph: Any, port_id: str, flow: float, inflow: bool) -> None:
        """``<id>.m_dot`` on the linear chain of two-port elements running
        outward from a port to its boundary, signed by each element's
        orientation (towards the port on an inlet, away on an outlet). Stops
        at a branch, a multi-port element or a boundary."""
        from combaero.network.components import MassFlowBoundary, PressureBoundary

        node_id, prev = port_id, self
        for _ in range(32):
            around = graph.get_upstream_elements(node_id) + graph.get_downstream_elements(node_id)
            elems = [e for e in around if e is not prev and e is not self]
            if len(elems) != 1:
                return
            elem = elems[0]
            frm, to = getattr(elem, "from_node", None), getattr(elem, "to_node", None)
            if node_id not in (frm, to) or frm is None or to is None:
                return
            signed = flow if (to == node_id) == inflow else -flow
            unk = f"{elem.id}.m_dot"
            if unk in elem.unknowns():
                self._seed(elem, unk, signed)
            node_id = frm if to == node_id else to
            if isinstance(graph.nodes.get(node_id), (PressureBoundary, MassFlowBoundary)):
                return
            prev = elem

    def _seed_from_closed_forms(self, graph: Any) -> None:
        from combaero.network._ejector_huang1999 import (
            choked_mass_flow,
            critical_back_pressure,
            entrainment_ratio,
            huang_critical_back_pressure,
        )
        from combaero.network.components import MassFlowBoundary, PressureBoundary

        pb = self._port_boundary(graph, self.primary_node)
        sb = self._port_boundary(graph, self.secondary_node)
        ob = self._port_boundary(graph, self.outlet_node)

        # Nominal air properties for the SEED only (exact gamma/R not needed).
        X = cb.species.dry_air()
        t_g = float(getattr(pb, "Tt", 440.0))
        t_e = float(getattr(sb, "Tt", 280.0))
        try:
            gamma = choke_plane_gamma(t_e, X)
        except Exception:  # noqa: BLE001
            gamma = 1.4
        r_gas = float(cb.specific_gas_constant(X))

        # Primary stagnation: direct for a pressure boundary; inverse choked
        # flow for a mass boundary (m_dot is linear in p0).
        p_g = None
        if isinstance(pb, PressureBoundary):
            p_g = float(pb.Pt)
        elif isinstance(pb, MassFlowBoundary):
            unit = choked_mass_flow(1.0, t_g, self.throat_area, gamma, r_gas, eta=ETA_P)
            if unit > 0.0:
                p_g = float(pb.m_dot) / unit
        p_e = float(sb.Pt) if isinstance(sb, PressureBoundary) else None
        p_o = float(ob.Pt) if isinstance(ob, PressureBoundary) else None

        # A mass-driven primary whose choked-inverse Pt is below the back
        # pressure cannot choke: the physical branch is the subsonic jet pump,
        # primary Pt floating to ~the back pressure, P_py just below it. The
        # choked-inverse seed there parks Newton in the spurious double-choked
        # basin (Pt_primary < Pb), from which it stalls.
        jetpump = (
            isinstance(pb, MassFlowBoundary) and p_g is not None and p_o is not None and p_g < p_o
        )
        if jetpump:
            p_g = p_o

        for port_id, value in (
            (self.primary_node, p_g),
            (self.secondary_node, p_e),
            (self.outlet_node, p_o),
        ):
            node = graph.nodes.get(port_id)
            if node is None or value is None or value <= 0.0:
                continue
            existing = getattr(node, "initial_guess", {}) or {}
            if any(
                k.endswith((".P", ".Pt")) for k in existing if (port_id, k) not in self._auto_seeds
            ):
                continue  # the user seeded this port
            self._seed(node, f"{port_id}.P", value)
            self._seed(node, f"{port_id}.Pt", value)

        # An inherited A_3 resolves in resolve_topology, from the outlet pipe;
        # a seed taken before that (the GUI seeds at build) uses the same
        # value, or every closed form below would fail on area_ratio_mix.
        a3 = self.mixing_area
        if a3 is None:
            a3 = self._outlet_pipe_area(graph)[0] or DEFAULT_EJECTOR_MIXING_AREA
        geom = EjectorGeometry(self.area_ratio_nozzle, float(a3) / self.throat_area)
        p_jct = None
        if jetpump and p_o is not None:
            p_jct = 0.99 * p_o
        elif p_g and p_e:
            try:
                p_jct = float(
                    (
                        huang_critical_back_pressure
                        if self.anchor_closure == "huang1999"
                        else critical_back_pressure
                    )(p_g, t_g, p_e, t_e, geom, gamma, r_gas, self.recovery_efficiency).p_c_pa
                )
            except Exception:  # noqa: BLE001
                p_jct = None
        if p_jct is not None and p_jct > 0.0:
            self._seed(self, f"{self.id}.P_jct", p_jct)

        # Port-chain flows (#567): left to the generic defaults (tens of kg/s
        # for a ~0.1 kg/s ejector) a channel on a port sends the cold start to
        # Mach-20 port states. Critical regime only: at the jet-pump seed
        # (p_g = p_o, no drive) the closed form gives omega ~ 0 and a
        # zero-entrainment seed parks Newton there.
        if jetpump:
            return
        mp = None
        if isinstance(pb, MassFlowBoundary):
            mp = float(pb.m_dot)
        elif p_g:
            mp = float(choked_mass_flow(p_g, t_g, self.throat_area, gamma, r_gas, eta=ETA_P))
        omega = None
        if mp and p_g and p_e:
            try:
                omega = float(entrainment_ratio(p_g, t_g, p_e, t_e, geom, gamma).omega)
            except Exception:  # noqa: BLE001
                omega = None
        if omega is not None and not (0.0 < omega < 50.0):
            omega = None  # an out-of-envelope closed form is no seed
        if mp:
            self._seed_port_chain_flows(graph, self.primary_node, mp, inflow=True)
            if omega is not None:
                self._seed_port_chain_flows(graph, self.secondary_node, omega * mp, inflow=True)
                self._seed_port_chain_flows(
                    graph, self.outlet_node, (1.0 + omega) * mp, inflow=False
                )

    def _size_outlet_port_node(self, graph: Any) -> None:
        """An auto-sized outlet port chamber takes A_3 unless a pipe sized it.

        Nodes resolve before elements, so an outlet pipe whose D is itself
        inherited from this ejector is still unsized when the port chamber
        looks for an area. The chamber's search then hopped across the ejector
        to ANOTHER port's chamber -- in the GUI default, the 0.05 m^2 primary
        inlet chamber -- and the outlet ran as a near-plenum. The outlet port
        face is the mixing-duct exit, so A_3 is its area; a chamber sized by an
        explicit outlet pipe (source "channel ...") keeps the pipe's area.
        """
        node = graph.nodes.get(self.outlet_node)
        if not isinstance(node, MomentumChamberNode) or not getattr(node, "_auto_area", False):
            return
        if str(getattr(node, "_area_source", "")).startswith("channel"):
            return
        if self.mixing_area is None:
            return
        source = f"EjectorElement '{self.id}' (A_3)"
        node.area = float(self.mixing_area)
        node.Dh = 2.0 * math.sqrt(node.area / math.pi)
        node._area_source = source
        node._Dh_source = source

    def _check_mixing_against_outlet_pipe(self, graph: Any) -> None:
        """An explicit mixing area wider than the outlet pipe is refused
        (#1023): the mixed stream would contract into the pipe with no loss
        at the diffuser exit. A wider pipe is a plausible diffuser exit and
        stays accepted; equal (to NEIGHBOUR_AREA_RTOL on the equivalent
        diameter) is the duct continuing. The pipe is the single channel at
        the outlet node (_attached_pipe), through the port chamber."""
        hit = _attached_pipe(graph, self.outlet_node, self)
        if hit is None or self.mixing_area is None:
            return
        pipe, A_pipe = hit
        if math.sqrt(self.mixing_area / A_pipe) <= 1.0 + NEIGHBOUR_AREA_RTOL:
            return
        raise ValueError(
            f"EjectorElement '{self.id}': mixing_area {self.mixing_area:.6g} m^2 is larger "
            f"than its outlet pipe, channel '{pipe.id}' at '{self.outlet_node}' "
            f"(D = {pipe.diameter:.6g} m, A = {A_pipe:.6g} m^2): the mixed stream would "
            "contract into the pipe with no loss. Make the outlet pipe at least the "
            "mixing area (a wider pipe is a diffuser exit), or leave mixing_area unset "
            "to inherit the pipe's."
        )

    def _outlet_pipe_area(self, graph: Any) -> tuple[float | None, str]:
        """The area of the element on the outside of the outlet port, unless
        that element itself took its area from this ejector (which would turn a
        default into an "inherited" value on the next pass)."""
        own = f"'{self.id}'"
        for elem in graph.get_upstream_elements(self.outlet_node) + graph.get_downstream_elements(
            self.outlet_node
        ):
            if elem is self or own in str(getattr(elem, "_area_source", "")):
                continue
            area = _element_area_facing(elem, self.outlet_node)
            if area is not None:
                return (
                    area,
                    f"{'channel' if type(elem).__name__ == 'ChannelElement' else type(elem).__name__} '{elem.id}'",
                )
        return None, ""

    @property
    def area_ratio_nozzle(self) -> float:
        """A_p1 / A_t."""
        return self.nozzle_exit_area / self.throat_area

    @property
    def area_ratio_mix(self) -> float:
        """A_3 / A_t."""
        return self.mixing_area / self.throat_area

    @property
    def geometry(self) -> EjectorGeometry:
        return EjectorGeometry(self.area_ratio_nozzle, self.area_ratio_mix)

    def row_scale_kinds(self) -> list[str]:
        """Rows 0-2 (choked flow, entrainment closure, mass conservation)
        are mass-flow-valued; row 3 (P_jct - P_c*) is pressure-valued --
        the OPPOSITE pattern from the base class's own impulse-CV rows.
        See the module docstring's residual design section and
        RESIDUAL-ROW-SCALING (NetworkSolver._build_residual_scales)."""
        if self.regime == "incompressible":
            # Nozzle and suction rows are pressure balances there (#588).
            return ["p", "p", "mdot", "p"]
        return ["mdot", "mdot", "mdot", "p"]

    # smootherstep window for s_choke, on mp / choked_mass_flow(Pt_p): s = 0
    # when clearly unchoked (jet pump), 1 exactly at the primary-choke threshold.
    _S_CHOKE_LO = 0.90
    _S_CHOKE_HI = 0.999
    # Handover band, centred on 1, for two smootherstep windows on outlet.Pt:
    # over recovery * P_c* (s_sub: the plateau omega_crit -> the subcritical
    # droop) and over recovery * P_b0 (s_b0: outlet free -> pinned, #619).
    _S_SUB_LO = 0.98
    _S_SUB_HI = 1.02
    _EPS_FRAC = 1.0e-3  # C-D nozzle / secondary-flux smooth-min rounding

    def residuals(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        if self.regime == "incompressible":
            return self._incompressible_residuals(states, P_jct, port_mdots)
        primary, secondary, outlet = states

        # gamma/r_gas: frozen coefficients (see module docstring). Recomputed
        # each call so the converged root is exact; the closure Jacobians hold
        # both fixed (drops only the tiny implicit d(gamma)/d(Tt) term).
        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        arn = self.area_ratio_nozzle
        arm = self.area_ratio_mix
        A_t = self.throat_area
        A_e = self.nozzle_exit_area
        A_s = self.mixing_area - self.nozzle_exit_area
        rec = self.recovery_efficiency

        # Physical port flows. port_mdots are junction-sign convention
        # (positive = out); primary/secondary are inlets (sign -1), so their
        # physical flow is the negation. The owned unknown ``P_jct`` is
        # repurposed as the mixing-plane static pressure P_py (see the module
        # docstring's residual-design section).
        mp_v = -float(port_mdots[0])
        ms_v = -float(port_mdots[1])
        mo_v = float(port_mdots[2])

        # The 4-row residual system + its 9-column analytic Jacobian are
        # assembled in C++ (whole-element (f, J), matching the base
        # MultiPortChamberBase / TeeJunctionElement practice): the scalar
        # closures are chained through a forward-mode DualN<9> across the
        # regime blend. See src/ejector.cpp / OPERATING_REGIMES_DESIGN.md sec
        # 6c. This Python is only the relabeling shim.
        rj = _solver_tools.ejector_element_residuals_and_jacobian(
            mp_v,
            ms_v,
            mo_v,
            primary.Pt,
            primary.Tt,
            secondary.Pt,
            secondary.Tt,
            outlet.Pt,
            P_jct,
            arn,
            arm,
            A_t,
            A_e,
            A_s,
            gamma,
            r_gas,
            ETA_P,
            ETA_S,
            rec,
            self._EPS_FRAC,
            self._S_CHOKE_LO,
            self._S_CHOKE_HI,
            self._S_SUB_LO,
            self._S_SUB_HI,
            anchor_closure=self.anchor_closure,
        )

        # Map the 9 Jacobian seeds (in the C++ unknown order) to solver columns.
        # Inlet flows: the physical flow is -(port unknown), so d(phys)/d(unknown)
        # = -port_sign; the outlet unknown IS the physical flow (+port_sign).
        seed_cols: list[tuple[str, float]] = [
            (f"{self._port_element_ids[0]}.m_dot", -self._port_signs[0]),  # 0 mp
            (f"{self._port_element_ids[1]}.m_dot", -self._port_signs[1]),  # 1 ms
            (f"{self._port_element_ids[2]}.m_dot", self._port_signs[2]),  # 2 mdot_out
            (f"{self.primary_node}.Pt", 1.0),  # 3 p_g
            (f"{self.primary_node}.Tt", 1.0),  # 4 t_g
            (f"{self.secondary_node}.Pt", 1.0),  # 5 p_e
            (f"{self.secondary_node}.Tt", 1.0),  # 6 t_e
            (f"{self.outlet_node}.Pt", 1.0),  # 7 p_out
            (f"{self.id}.P_jct", 1.0),  # 8 P_py (owned)
        ]
        jac: dict[int, dict[str, float]] = {}
        for row in range(4):
            grad = rj.jacobian[row]
            cols: dict[str, float] = {}
            for seed, (name, coeff) in enumerate(seed_cols):
                g = grad[seed]
                if g != 0.0:
                    cols[name] = cols.get(name, 0.0) + coeff * g
            jac[row] = cols
        return (list(rj.residuals), jac)

    def _incompressible_residuals(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float],
    ) -> tuple[list[float], dict[int, dict[str, float]]]:
        """Sanger's jet pump (#588): C++ assembles the four rows and their
        Jacobian over (mp, ms, mo, Pt_p, v_p, Pt_s, v_s, Pt_o, p3); here the
        port densities (each port chamber's own closure, as every junction
        reads them) are chained into the v columns."""
        primary, secondary, outlet = states
        dens = self.port_densities(states)
        rho_p, drho_p = dens[0]
        rho_s, drho_s = dens[1]
        rj = _solver_tools.ejector_jetpump_incompressible_residuals_and_jacobian(
            -float(port_mdots[0]),
            -float(port_mdots[1]),
            float(port_mdots[2]),
            float(primary.Pt),
            1.0 / rho_p,
            float(secondary.Pt),
            1.0 / rho_s,
            float(outlet.Pt),
            float(P_jct),
            self.nozzle_exit_area,
            self.mixing_area - self.nozzle_exit_area,
            self.mixing_area,
            self.k_nozzle,
            self.k_suction,
            self.k_throat,
            self.k_diffuser,
            self.recovery_efficiency,
        )
        flow_cols: list[tuple[int, str, float]] = [
            (0, f"{self._port_element_ids[0]}.m_dot", -self._port_signs[0]),
            (1, f"{self._port_element_ids[1]}.m_dot", -self._port_signs[1]),
            (2, f"{self._port_element_ids[2]}.m_dot", self._port_signs[2]),
            (3, f"{self.primary_node}.Pt", 1.0),
            (5, f"{self.secondary_node}.Pt", 1.0),
            (7, f"{self.outlet_node}.Pt", 1.0),
            (8, f"{self.id}.P_jct", 1.0),
        ]
        # d(v)/d(name) = -d(rho)/d(name) / rho^2
        v_cols = [
            (4, {k: -d / (rho_p * rho_p) for k, d in drho_p.items()}),
            (6, {k: -d / (rho_s * rho_s) for k, d in drho_s.items()}),
        ]
        jac: dict[int, dict[str, float]] = {}
        for row in range(4):
            grad = rj.jacobian[row]
            cols: dict[str, float] = {}
            for seed, name, coeff in flow_cols:
                g = grad[seed]
                if g != 0.0:
                    cols[name] = cols.get(name, 0.0) + coeff * g
            for seed, chain in v_cols:
                g = grad[seed]
                if g == 0.0:
                    continue
                for name, dv in chain.items():
                    cols[name] = cols.get(name, 0.0) + g * dv
            jac[row] = cols
        return list(rj.residuals), jac

    def verify_solution_consistent(
        self,
        sol: dict[str, float],
        rel_tol: float = EJECTOR_REVERSE_FLOW_TOL,
    ) -> bool:
        """Reject a root the model was not derived for (#619).

        The element models the critical, subcritical and unchoked-primary
        (jet-pump) regimes in one C1 residual system; the old ``outlet.Pt <=
        P_c*`` check demoted valid subcritical/jet-pump points and stays
        dropped (OPERATING_REGIMES_DESIGN.md sec 6b). But every one of those
        regimes has its streams going forward: primary and suction INTO the
        mixing plane, the mixed stream OUT of the outlet. Past the jet pump's
        shut-off head the suction would reverse, and the rows' C1 continuation
        there has roots no ejector produces: seed 36 of the stress family
        "converged" at omega = -169 with 0.47 kg/s entering the 1.7 bar
        outlet, climbing to a 4 bar mixing plane and leaving through the 1 bar
        suction. So a root with any port flowing backwards, beyond ``rel_tol``
        of the largest port flow, is not reported as a solution. Missing keys
        mean True, as for the junction verifiers.
        """
        flows: list[float] = []
        for i in range(self.N):
            elem_id = self._port_element_ids[i]
            key = f"{elem_id}.m_dot"
            if not elem_id or key not in sol:
                return True
            # Junction convention: positive = OUT of the mixing plane.
            flows.append(float(self._port_signs[i]) * float(sol[key]))
        m_primary_in, m_suction_in, m_outlet_out = -flows[0], -flows[1], flows[2]
        tol = rel_tol * max(abs(f) for f in flows)
        return m_primary_in >= -tol and m_suction_in >= -tol and m_outlet_out >= -tol

    def diagnostics(
        self,
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None = None,
    ) -> dict[str, float]:
        """Parent per-port fields plus ejector-specific outputs: entrainment
        ratio, live gamma/R, and the critical-vs-actual back-pressure margin.

        Pressures follow the result-key convention (#773): ``P`` static,
        ``Pt`` total. ``P_jct`` is the owned unknown, the mixing-plane STATIC
        pressure P_py at which both streams arrive; ``Pt_out`` the outlet
        port's total; ``P_c_star`` / ``P_b0`` the critical and dead-head back
        pressures (outlet totals). There is no ``Pt_jct``: at the mixing plane
        the primary and secondary streams share the static P_py at different
        velocities (the primary supersonic, the secondary sonic in critical
        mode), so the plane has no single total pressure; the mixed stream's
        total exists only downstream of the mixing, which is not the plane
        P_jct describes. ``entrainment_ratio`` is m_dot_secondary /
        m_dot_primary (Huang's omega).
        """
        diag = super().diagnostics(states, P_jct, port_mdots)
        primary, secondary, outlet = states
        if self.regime == "incompressible":
            return self._incompressible_diagnostics(diag, states, P_jct, port_mdots)

        gamma = choke_plane_gamma(secondary.Tt, secondary.X)
        r_gas = float(cb.specific_gas_constant(secondary.X))
        diag["gamma"] = gamma
        diag["r_gas"] = r_gas
        diag["Pt_out"] = float(outlet.Pt)
        # Flag and clamp (#945): the nozzle-exit Mach M_p1 the closures use is
        # the supersonic root of A/A* = area_ratio_nozzle on the bracket
        # [1, 50]. Outside it the bisection returns the bracket end; +1 (area
        # ratio above A/A*(M = 50)) or -1 (below) says the result rests on a
        # clamped M_p1, 0 that it does not.
        diag["nozzle_exit_mach_out_of_range"] = float(
            _solver_tools.ejector_nozzle_exit_mach(self.area_ratio_nozzle, gamma).out_of_range
        )

        # Operating regime: s_choke = smootherstep(mp/cap), reported so the GUI
        # can label critical (1) / jet-pump (0) / transition, plus the actual
        # entrainment ratio and the choked-mode reference quantities. The
        # critical closures (entrainment_ratio_critical, P_c*) assume a choked primary and
        # return NaN in the jet-pump regime, so they are reported only when
        # meaningful (s_choke > 0).
        if port_mdots is not None and len(port_mdots) == self.N:
            mp = -float(port_mdots[0])
            ms = -float(port_mdots[1])
            cap = _solver_tools.ejector_choked_mass_flow_and_jacobian(
                primary.Pt, primary.Tt, self.throat_area, gamma, r_gas, ETA_P
            )[0]
            t = (mp / cap - self._S_CHOKE_LO) / (self._S_CHOKE_HI - self._S_CHOKE_LO)
            t = min(max(t, 0.0), 1.0)
            s_choke = t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
            diag["s_choke"] = s_choke
            diag["m_dot_primary"] = mp
            diag["m_dot_secondary"] = ms
            diag["m_dot_outlet"] = float(port_mdots[2])
            # ACTUAL entrainment ratio -- regime-independent, always valid.
            diag["entrainment_ratio"] = ms / mp if mp != 0.0 else float("nan")
            # Regime label. critical_mode = 1 only when the primary is choked
            # AND the outlet sits below the ejector's critical back pressure
            # (double-choked plateau); 0 in the subcritical-droop and unchoked
            # jet-pump regimes, where the outlet pressure sets the operating
            # point.
            critical_mode = 0.0
            # The choking anchors exist only for a choked primary: NaN, not
            # absent, otherwise (a result's keys never follow its regime).
            diag["entrainment_ratio_critical"] = math.nan
            diag["P_c_star"] = math.nan
            diag["P_b0"] = math.nan
            if s_choke > 0.0:
                omega_c = self._critical_omega(primary, secondary, gamma)
                p_c_star = self._critical_back_pressure(primary, secondary, gamma, r_gas)
                diag["entrainment_ratio_critical"] = omega_c
                diag["P_c_star"] = p_c_star
                # Upper anchor of the subcritical droop (#619): omega falls
                # linearly from entrainment_ratio_critical at P_c* to 0 here.
                diag["P_b0"] = self._dead_head_back_pressure(primary, secondary, gamma, r_gas)
                if s_choke >= 0.5 and math.isfinite(p_c_star) and float(outlet.Pt) <= p_c_star:
                    critical_mode = 1.0
            diag["critical_mode"] = critical_mode
        else:
            diag["entrainment_ratio"] = self._critical_omega(primary, secondary, gamma)
            diag["P_c_star"] = self._critical_back_pressure(primary, secondary, gamma, r_gas)
        return diag

    def _incompressible_diagnostics(
        self,
        diag: dict[str, float],
        states: list[NetworkMixtureState],
        P_jct: float,
        port_mdots: list[float] | None,
    ) -> dict[str, float]:
        """Jet-pump outputs (#588): flows, entrainment ratio, Sanger's head
        ratio N = (Pt_o - Pt_s) / (Pt_p - Pt_o), efficiency M N, and the primary
        jet's Mach for the validity report. No choking quantities: the
        incompressible jet pump has none."""
        primary, secondary, outlet = states
        diag["regime_incompressible"] = 1.0
        diag["critical_mode"] = 0.0
        diag["Pt_out"] = float(outlet.Pt)
        if port_mdots is None or len(port_mdots) != self.N:
            return diag
        mp, ms = -float(port_mdots[0]), -float(port_mdots[1])
        diag["m_dot_primary"] = mp
        diag["m_dot_secondary"] = ms
        diag["m_dot_outlet"] = float(port_mdots[2])
        diag["entrainment_ratio"] = ms / mp if mp != 0.0 else float("nan")
        rho_p = self.port_densities(states)[0][0]
        a_p = float(cb.speed_of_sound(primary.T, primary.X))
        diag["mach_jet"] = abs(mp) / (rho_p * self.nozzle_exit_area) / a_p if a_p > 0 else 0.0
        drive = float(primary.Pt) - float(outlet.Pt)
        # Undefined without a drive (or a primary flow): NaN, not absent.
        diag["head_ratio_N"] = math.nan
        diag["efficiency"] = math.nan
        if drive != 0.0:
            N = (float(outlet.Pt) - float(secondary.Pt)) / drive
            diag["head_ratio_N"] = N
            if mp != 0.0:
                diag["efficiency"] = (ms / mp) * N
        return diag

    def _critical_omega(
        self, primary: NetworkMixtureState, secondary: NetworkMixtureState, gamma: float
    ) -> float:
        return float(
            _solver_tools.ejector_entrainment_ratio_and_jacobian(
                primary.Pt,
                primary.Tt,
                secondary.Pt,
                secondary.Tt,
                self.area_ratio_nozzle,
                self.area_ratio_mix,
                gamma,
            ).value.omega
        )

    def _dead_head_back_pressure(
        self,
        primary: NetworkMixtureState,
        secondary: NetworkMixtureState,
        gamma: float,
        r_gas: float,
    ) -> float:
        closure = (
            _solver_tools.ejector_sudden_expansion_dead_head_and_jacobian
            if self.anchor_closure == "huang1999"
            else _solver_tools.ejector_dead_head_back_pressure_and_jacobian
        )
        return float(
            closure(
                primary.Pt,
                primary.Tt,
                secondary.Pt,
                secondary.Tt,
                self.area_ratio_nozzle,
                self.area_ratio_mix,
                gamma,
                r_gas,
                self.recovery_efficiency,
            ).value.p_c_pa
        )

    def _critical_back_pressure(
        self,
        primary: NetworkMixtureState,
        secondary: NetworkMixtureState,
        gamma: float,
        r_gas: float,
    ) -> float:
        closure = (
            _solver_tools.ejector_huang_critical_back_pressure_and_jacobian
            if self.anchor_closure == "huang1999"
            else _solver_tools.ejector_critical_back_pressure_and_jacobian
        )
        return float(
            closure(
                primary.Pt,
                primary.Tt,
                secondary.Pt,
                secondary.Tt,
                self.area_ratio_nozzle,
                self.area_ratio_mix,
                gamma,
                r_gas,
                self.recovery_efficiency,
            ).value.p_c_pa
        )
