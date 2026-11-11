"""Impingement array, channel exit: one configuration, built from rows (#465).

The configuration Florschuetz, Truman and Metzger (1981) measured: a jet
plate fed from a supply plenum, impinging on a target, the spent air leaving
in ONE direction down the gap between them. Optionally an initial (bypass)
crossflow enters the gap's upstream end (#467): Florschuetz, Metzger and Su
(1982, NASA CR-3630). Spent air leaving through the target (#468) is a
different configuration and gets its own element.

``ImpingementArray`` is the configuration, not a network element: it adds the
rows to a network, one ``ImpingementPlateElement`` and one
``ImpingementCrossflowElement`` per row, joined by internal crossflow nodes,
so the solver sees nothing new. Ids are namespaced under the array's own id::

    supply --+------------------+------------------+
             |                  |                  |
         {id}__p1           {id}__p2           {id}__pN
             |                  |                  |
         {id}__c1 --{id}__x1-- {id}__c2 -- ... -- {id}__cN --{id}__xN-- exit

With a bypass the closed end becomes an inlet, ``bypass --{id}__x0-- {id}__c1``
(see ``ImpingementCrossflowElement``'s ``bypass_inlet``): wired from another
array's exit the gap continues with its momentum; from any other node the
bypass is accelerated from that node's total pressure, so a feed duct's
dynamic head arrives with it.

``add_wall`` distributes one thermal wall to the hot side across the rows,
and ``summarize`` collects the per-row diagnostics under the array's id.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Literal

from .components import (
    CasingRadiation,
    GasRadiation,
    HoleCd,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    PlenumNode,
    ThermalWall,
    WallLayer,
    _max_over_min,
)
from .hole_cd import _hole_cd_spec

if TYPE_CHECKING:
    from .graph import FlowNetwork

# A bypass this much hotter or colder than the supply moves the target's
# reference temperature (CR-3936's eta_r, see ImpingementPlateElement); below
# it the reference is the supply's to solver noise (see summarize).
BYPASS_DT_TOL_K = 1.0

# A row whose plate total-pressure drop is below this fraction of the gap's
# dynamic head there is below what the crossflow chain resolves (#729): the
# station scheme is consistent to O(M^2 dm/m) of that head. On CR-3630 Table
# D.2 (B(5,4,1) at mc/mj 0.5 and 1.01, M 0.22) the dead rows sit at up to 0.4%
# of it, their flows alternating in sign from row to row; every other row is
# above 6.9%. 1% falls in that gap. It is a resolution estimate, not a
# measured band. The row's flow, Nu and eta_r are reported as unresolved
# (see summarize).
UNRESOLVED_ROW_DPT_FRACTION = 0.01


@dataclass
class ImpingementArray:
    """A jet-plate impingement array whose spent air exits down the channel.

    Parameters
    ----------
    id : str
        The array's id; every row element and node is namespaced under it.
    n_rows : int
        Spanwise rows, counted from the closed (upstream) end.
    d_jet, xn_d, yn_d, z_d, span, plate_thickness, pattern, Nu_multiplier :
        As ``ImpingementPlateElement``; every row shares them.
    hole_cd : HoleCd or None
        The holes' Cd, as ``ImpingementPlateElement(hole_cd=...)``; every row
        shares it. None: the row's defaults ('fixed' at Florschuetz's 0.79,
        the gap's crossflow on the Cd by Rogers & Hersh 1976, #800).
        ``HoleCd(exit_crossflow=...)`` picks the crossflow law:
        'rogers_hersh_1976' ("auto"), 'florschuetz_metzger_su_1982' (CR-3630
        Table 4.1, the plates it covers) or 'none' (the constant Cd).
    correlation, Cd, cd_crossflow :
        Deprecated since 0.12.0, removed in 0.13.0: ``hole_cd=HoleCd(
        correlation=..., Cd=..., exit_crossflow=...)``.
    roughness : float
        Of the target and jet plate in the gap, for the crossflow friction.
    bypass_length : float or None
        Lead-in from the bypass inlet to row 1 [m]. None: ``xn/2``, where
        CR-3630's array begins (half a pitch ahead of row 1).
    bypass_entry_K : float
        Entry loss into the gap, on its dynamic head: 0 for a bypass arriving
        along a duct, 0.5 for an abrupt sharp-edged entry from a large
        plenum. Unused when the gap continues an upstream array's.
    """

    id: str
    n_rows: int
    d_jet: float
    xn_d: float
    yn_d: float
    z_d: float
    span: float
    plate_thickness: float
    pattern: Literal["inline", "staggered"] = "inline"
    correlation: str | None = None
    Cd: float | None = None
    Nu_multiplier: float = 1.0
    roughness: float = 0.0
    bypass_length: float | None = None
    bypass_entry_K: float = 0.0
    cd_crossflow: Literal["rogers_hersh_1976", "florschuetz_metzger_su_1982", "none"] | None = None
    hole_cd: HoleCd | None = None
    _rows_hole_cd: HoleCd = field(default_factory=HoleCd, init=False, repr=False)
    _bypass: str | None = field(default=None, init=False, repr=False)
    _rows: list[tuple[str, str, str]] = field(default_factory=list, init=False, repr=False)

    def __post_init__(self) -> None:
        if int(self.n_rows) != self.n_rows or self.n_rows < 1:
            raise ValueError(f"ImpingementArray {self.id!r}: n_rows must be a positive integer")
        if self.bypass_length is not None and self.bypass_length <= 0.0:
            raise ValueError(f"ImpingementArray {self.id!r}: bypass_length must be positive")
        if self.bypass_entry_K < 0.0:
            raise ValueError(f"ImpingementArray {self.id!r}: bypass_entry_K must be non-negative")
        # The rows' HoleCd: given, or built from the deprecated fields (one
        # DeprecationWarning, here, not per row). A Cd with a correlation
        # (#872) is refused here, at construction, by HoleCd. Every field
        # keeps what was passed (so dataclasses.replace rebuilds the same
        # array); the resolved spec the rows take is _rows_hole_cd.
        self._rows_hole_cd = _hole_cd_spec(
            f"ImpingementArray {self.id!r}",
            self.hole_cd,
            {"correlation": self.correlation, "Cd": self.Cd, "cd_crossflow": self.cd_crossflow},
            stacklevel=4,
        )

    def plate_id(self, row: int) -> str:
        return f"{self.id}__p{row}"

    def segment_id(self, row: int) -> str:
        return f"{self.id}__x{row}"

    def node_id(self, row: int) -> str:
        return f"{self.id}__c{row}"

    def add_to(self, net: FlowNetwork, supply: str, exit: str, bypass: str | None = None) -> None:
        """Add the rows between ``supply`` (the jet plenum) and ``exit``.

        ``bypass``, when given, is the node an initial crossflow enters the
        gap from; None is the closed upstream end of Florschuetz et al. (1981).
        """
        n = int(self.n_rows)
        for i in range(1, n + 1):
            net.add_node(PlenumNode(self.node_id(i)))
        self._rows = []
        self._bypass = bypass
        if bypass is not None:
            net.add_element(
                ImpingementCrossflowElement(
                    self.segment_id(0),
                    bypass,
                    self.node_id(1),
                    length=(
                        self.bypass_length
                        if self.bypass_length is not None
                        else 0.5 * self.xn_d * self.d_jet
                    ),
                    height=self.z_d * self.d_jet,
                    span=self.span,
                    roughness=self.roughness,
                    bypass_inlet=True,
                    entry_K=self.bypass_entry_K,
                )
            )
        for i in range(1, n + 1):
            to = self.node_id(i + 1) if i < n else exit
            net.add_element(
                ImpingementPlateElement(
                    self.plate_id(i),
                    supply,
                    self.node_id(i),
                    d_jet=self.d_jet,
                    xn_d=self.xn_d,
                    yn_d=self.yn_d,
                    z_d=self.z_d,
                    span=self.span,
                    plate_thickness=self.plate_thickness,
                    pattern=self.pattern,
                    row=i,
                    Nu_multiplier=self.Nu_multiplier,
                    hole_cd=self._rows_hole_cd,
                )
            )
            net.add_element(
                ImpingementCrossflowElement(
                    self.segment_id(i),
                    self.node_id(i),
                    to,
                    length=self.xn_d * self.d_jet,
                    height=self.z_d * self.d_jet,
                    span=self.span,
                    roughness=self.roughness,
                )
            )
            self._rows.append((self.plate_id(i), self.segment_id(i), self.node_id(i)))

    def add_wall(
        self,
        net: FlowNetwork,
        wall_id: str,
        hot_element: str,
        layers: list[WallLayer],
        R_fouling: float = 0.0,
        array_is_side_a: bool = False,
        radiation_hot: GasRadiation | CasingRadiation | None = None,
        radiation_target: GasRadiation | CasingRadiation | None = None,
    ) -> list[str]:
        """One wall between ``hot_element`` and the array's target, as one
        ``ThermalWall`` per row, each on that row's target footprint.

        The hot side is evaluated once per row at its own state, so every row
        sees the same hot-gas h and T_aw: a hot side that changes along the
        array needs its own segments. Returns the per-row wall ids.

        Radiation (#793, #1071) goes on each row's wall as ThermalWall's
        ``radiation_a``/``radiation_b``: ``radiation_hot`` on the hot
        element's face (a GasRadiation from a flame, or a CasingRadiation),
        ``radiation_target`` on the target face (typically a CasingRadiation
        toward the jet plate). Every row's wall gets the same specs.
        """
        ids = []
        for i in range(1, int(self.n_rows) + 1):
            plate = self.plate_id(i)
            a, b = (plate, hot_element) if array_is_side_a else (hot_element, plate)
            rad_a, rad_b = (
                (radiation_target, radiation_hot)
                if array_is_side_a
                else (radiation_hot, radiation_target)
            )
            wid = f"{wall_id}__r{i}"
            net.add_wall(
                ThermalWall(
                    id=wid,
                    element_a=a,
                    element_b=b,
                    layers=list(layers),
                    contact_area=net.elements[plate].surface.area,
                    R_fouling=R_fouling,
                    radiation_a=rad_a,
                    radiation_b=rad_b,
                )
            )
            ids.append(wid)
        net.thermal_coupling_enabled = True
        return ids

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The array's own diagnostics from its rows' (``__element_diag__``).

        Scalars for the array as a whole, plus a per-row list for each row
        quantity. The pressure block is supply-to-exit.

        ``rows_unresolved`` marks (1.0) a starved row: its plate drop is below
        ``UNRESOLVED_ROW_DPT_FRACTION`` of the gap's dynamic head, so its jet
        flow -- sign included -- and the Nu and eta_r that follow from it are
        below the crossflow model's resolution (#729). Such rows carry a tiny
        share of the flow; their heat transfer is the gap's channel flow,
        which the jet correlation does not model (#467).
        """
        n = int(self.n_rows)
        plates = [element_diag[self.plate_id(i)] for i in range(1, n + 1)]
        segments = [element_diag[self.segment_id(i)] for i in range(1, n + 1)]
        m_total = sum(p["m_dot"] for p in plates)
        out: dict[str, Any] = {
            "m_dot": float(m_total),
            "n_rows": float(n),
            "n_holes": float(sum(p["n_holes"] for p in plates)),
            "Pt_in": float(plates[0]["Pt_in"]),
            "P_out": float(segments[-1]["P_out"]),
            "dP": float(plates[0]["Pt_in"] - segments[-1]["P_out"]),
            "Re_j_min": float(min(p["Re_j"] for p in plates)),
            "Re_j_max": float(max(p["Re_j"] for p in plates)),
            "Gc_Gj_max": float(max(p["Gc_Gj"] for p in plates)),
            "Nu_mean": float(sum(p["Nu"] for p in plates) / n),
            "htc_mean": float(sum(p["htc"] for p in plates) / n),
            "Cd": float(plates[0]["Cd"]),
            "Cd_in_range": float(all(p["Cd_in_range"] for p in plates)),
            # Which Cd input is outside its source's range on any row, by
            # variable: information, never a warning (the plates' own keys).
            **{
                k: float(any(p[k] for p in plates))
                for k in (
                    "Cd_extrapolated_Re",
                    "Cd_extrapolated_l_d",
                    "Cd_extrapolated_edge",
                    "Cd_crossflow_extrapolated_velocity_ratio",
                    "Cd_crossflow_extrapolated_l_d",
                    "Cd_crossflow_extrapolated_grazing_velocity",
                )
            },
            "surface_extrapolated": float(any(p["surface_extrapolated"] for p in plates)),
            "jet_flow_nonuniformity": _max_over_min([p["m_dot"] for p in plates]),
        }
        for key in ("m_dot", "Re_j", "Gc_Gj", "Nu", "htc", "T_aw", "eta_r", "bypass_fraction"):
            out[f"rows_{key}"] = [float(p[key]) for p in plates]
        # The gap's dynamic head at each row, on the segment leaving it.
        unresolved = [
            float(
                abs(p["dPt"])
                < UNRESOLVED_ROW_DPT_FRACTION * 0.5 * s["rho_ref"] * s["velocity"] ** 2
            )
            for p, s in zip(plates, segments, strict=True)
        ]
        out["rows_unresolved"] = unresolved
        out["n_rows_unresolved"] = float(sum(unresolved))
        if self._bypass is not None:
            lead = element_diag[self.segment_id(0)]
            m_c = float(lead["m_dot"])
            dT = float(lead["Tt_in"]) - float(plates[0]["Tt_in"])
            out["m_dot_bypass"] = m_c
            out["bypass_to_jet_flow"] = m_c / m_total if m_total > 0.0 else float("inf")
            out["bypass_Tt_minus_supply"] = dT
            # A bypass at another temperature moves each row's reference
            # temperature through CR-3936's eta_r: a closure fitted to its
            # INLINE data. Outside that data -- a staggered array, or a row
            # whose geometry or Gc/Gj the data do not cover -- it is
            # extrapolated. 1 K is above solver noise, not a physical band.
            outside = self.pattern != "inline" or any(p["eta_r_extrapolated"] for p in plates)
            out["heat_transfer_outside_source"] = float(abs(dT) > BYPASS_DT_TOL_K and outside)
            # The gap draining back out of the inlet is outside CR-3630: the
            # lead-in then discharges at static continuity (channel_entry_drop
            # is inflow only), a dump nothing in the data checks. Reported,
            # not hidden.
            out["bypass_reversed"] = float(m_c < 0.0)
        return out
