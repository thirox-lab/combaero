"""Effusion liner: a duct-fed effusion wall as one configuration (#471).

A backside duct (an annulus between liner and casing) carries coolant from
``coolant_in`` to ``coolant_out``; its inner wall is an effusion plate that
bleeds coolant into the hot gas. Both sides have crossflow: the coolant runs
past the holes, the gas past their exits.

``EffusionLiner`` is the configuration, not a network element. It adds N
stations along the duct, each a ``PlenumNode`` (Pt = P: the node pressure is
the duct's static pressure, the drive a hole sees), joined by
``CrossflowSegmentElement`` s, and one ``EffusionPlateElement`` panel per
station into the gas node::

    coolant_in --s0-- b1 --s1-- b2 -- ... -- bN --sN-- coolant_out
                      |         |               |
                    panel1    panel2    ...   panelN     (into gas)

* Each station is a BLEED: the segments carry Bassett's separating
  straight-run momentum (kappa = 0.75, half per side; ``station_half_drop``).
* ``s0`` enters the duct from ``coolant_in`` as a reservoir (``entry_K``);
  ``sN`` dumps into ``coolant_out`` (static continuity).
* Each panel is DUCT-FED: its Cd is McGreehan-Schotsch with the supply-side
  U1/Vi from the duct's mean velocity at its station. Its wall is the #475
  plate wall; the coolant-side heat transfer is still Andrews' still-plenum
  correlation, flagged ``coolant_ht_plenum_assumed`` (crossflow-supply heat
  transfer is the next step, pending sources).

Ids are namespaced under the liner's id: ``{id}__b{i}``, ``{id}__s{i}``,
``{id}__p{i}``.

**A film provider (#695).** The liner registers its panels' film on the
network under its own id, so a ThermalWall downstream on the gas path can
name it (``film=<liner id>``). Each panel is one group of rows (Baldauf per
row, ``pitch_x`` apart, as ``EffusionPlateElement.film_span_average``) at
its own blowing ratio, density ratio and coolant temperature, layered row by
row along the gas path (``cb.film_span_average_baldauf_2002_groups``).
``film_direction`` says which way the gas runs past the panels: with the
coolant (station 1 upstream on the gas side) or against it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Literal

import combaero as cb

from .components import (
    CrossflowSegmentElement,
    EffusionPlateElement,
    FilmSource,
    GasRadiation,
    HoleCd,
    PlenumNode,
    _max_over_min,
)
from .hole_cd import _with_correlation

if TYPE_CHECKING:
    from .graph import FlowNetwork

# Bassett et al. (2001) measured the separating straight run at main/branch
# area ratios psi = F_C/F_B of 1 to 3. An effusion station bleeds through
# holes far smaller than the duct; K5 is psi-independent by derivation, so
# the station is used there, but flagged.
BASSETT_PSI_MEASURED_MAX = 3.0


@dataclass
class EffusionLiner:
    """A duct-fed effusion liner, N stations along the backside duct.

    Parameters
    ----------
    id : str
    n_segments : int
        Stations (panels) along the duct. More resolve the bleed's decline
        along the liner -- the first panel passes more than the last.
    length, width : float
        Liner panel length (along the backside flow) and width [m]. Each
        station's panel is ``length/n_segments`` by ``width``.
    duct_height : float
        Backside duct height [m]; its cross-section is ``duct_height*width``.
    hole_diameter, wall_thickness, pitch_x, pitch_y, angle_deg :
        As ``EffusionPlateElement``; every panel shares them.
    entry_K : float
        Loss on entering the duct from ``coolant_in`` (0.5 sharp).
    wall_conductivity, gas_film, gas_augmentation, turbulence_intensity,
    gas_heat_flux, internal_Nu_multiplier, edge_radius, material :
        As ``EffusionPlateElement``.
    film_direction : {"with_coolant", "against_coolant"}
        Which way the hot gas runs past the panels, for the film the liner
        lays on a wall downstream (#695): with the backside coolant (station
        1 is upstream on the gas side too) or against it.
    radiation : GasRadiation | None
        Gas radiation on every panel's hot face, as
        ``EffusionPlateElement(radiation=...)`` (#1071); None (default): none.
    hole_cd : HoleCd | None
        Every panel's hole Cd, as ``EffusionPlateElement(hole_cd=...)``. Its
        correlation left at None is 'McGreehanSchotsch', the liner's default
        (the duct-fed panels read its supply crossflow U1/Vi); None
        (default): ``HoleCd()``, i.e. that and every term at "auto".
    """

    id: str
    n_segments: int
    length: float
    width: float
    duct_height: float
    hole_diameter: float
    wall_thickness: float
    pitch_x: float
    pitch_y: float
    angle_deg: float = 90.0
    entry_K: float = 0.5
    roughness: float = 0.0
    edge_radius: float = 0.0
    wall_conductivity: float = 20.0
    gas_film: Literal["none", "baldauf_sellers"] = "none"
    gas_augmentation: float = 1.0
    turbulence_intensity: float = 0.05
    gas_heat_flux: float = 0.0
    internal_Nu_multiplier: float = 1.0
    film_direction: Literal["with_coolant", "against_coolant"] = "with_coolant"
    material: str | None = None
    radiation: GasRadiation | None = None
    hole_cd: HoleCd | None = None

    def __post_init__(self) -> None:
        if self.hole_cd is not None and not isinstance(self.hole_cd, HoleCd):
            raise ValueError(
                f"EffusionLiner {self.id!r}: hole_cd is a HoleCd, not {type(self.hole_cd).__name__}"
            )
        if int(self.n_segments) != self.n_segments or self.n_segments < 1:
            raise ValueError(f"EffusionLiner {self.id!r}: n_segments must be a positive integer")
        for name in ("length", "width", "duct_height"):
            if getattr(self, name) <= 0.0:
                raise ValueError(f"EffusionLiner {self.id!r}: {name} must be positive")
        if self.film_direction not in ("with_coolant", "against_coolant"):
            raise ValueError(
                f"EffusionLiner {self.id!r}: film_direction must be "
                "'with_coolant' or 'against_coolant'"
            )

    @property
    def duct_area(self) -> float:
        return self.duct_height * self.width

    @property
    def duct_Dh(self) -> float:
        return 2.0 * self.duct_height * self.width / (self.duct_height + self.width)

    def node_id(self, i: int) -> str:
        return f"{self.id}__b{i}"

    def segment_id(self, i: int) -> str:
        return f"{self.id}__s{i}"

    def panel_id(self, i: int) -> str:
        return f"{self.id}__p{i}"

    def add_to(
        self, net: FlowNetwork, coolant_in: str, coolant_out: str, gas: str | list[str]
    ) -> None:
        """Add the stations between ``coolant_in`` and ``coolant_out``, each
        panel discharging into ``gas`` (one node, or one per station for an
        axial gas-side gradient)."""
        n = int(self.n_segments)
        gas_nodes = [gas] * n if isinstance(gas, str) else list(gas)
        if len(gas_nodes) != n:
            raise ValueError(f"EffusionLiner {self.id!r}: give one gas node or {n}")
        dx = self.length / n
        bleed = cb.STATION_KAPPA_BLEED_BASSETT
        A, Dh = self.duct_area, self.duct_Dh
        hole_cd = _with_correlation(
            self.hole_cd or HoleCd(), "McGreehanSchotsch", f"EffusionLiner {self.id!r}"
        )
        for i in range(1, n + 1):
            net.add_node(PlenumNode(self.node_id(i)))
        # s0: reservoir entry into station 1.
        net.add_element(
            CrossflowSegmentElement(
                self.segment_id(0),
                coolant_in,
                self.node_id(1),
                length=0.5 * dx,
                area=A,
                Dh=Dh,
                roughness=self.roughness,
                entry_K=self.entry_K,
                to_kappa=bleed,
                next_seg=self.segment_id(1),
            )
        )
        for i in range(1, n + 1):
            last = i == n
            net.add_element(
                CrossflowSegmentElement(
                    self.segment_id(i),
                    self.node_id(i),
                    coolant_out if last else self.node_id(i + 1),
                    length=(0.5 if last else 1.0) * dx,
                    area=A,
                    Dh=Dh,
                    roughness=self.roughness,
                    from_kappa=bleed,
                    prev_seg=self.segment_id(i - 1),
                    to_kappa=None if last else bleed,
                    next_seg=None if last else self.segment_id(i + 1),
                )
            )
            net.add_element(
                EffusionPlateElement(
                    self.panel_id(i),
                    self.node_id(i),
                    gas_nodes[i - 1],
                    hole_diameter=self.hole_diameter,
                    wall_thickness=self.wall_thickness,
                    pitch_x=self.pitch_x,
                    pitch_y=self.pitch_y,
                    panel_length=dx,
                    panel_width=self.width,
                    angle_deg=self.angle_deg,
                    hole_cd=hole_cd,
                    edge_radius=self.edge_radius,
                    internal_Nu_multiplier=self.internal_Nu_multiplier,
                    wall_conductivity=self.wall_conductivity,
                    gas_film=self.gas_film,
                    gas_augmentation=self.gas_augmentation,
                    turbulence_intensity=self.turbulence_intensity,
                    gas_heat_flux=self.gas_heat_flux,
                    crossflow_segments=(self.segment_id(i - 1), self.segment_id(i)),
                    crossflow_area=A,
                    material=self.material,
                    radiation=self.radiation,
                )
            )
        panels = [net.elements[self.panel_id(i)] for i in range(1, n + 1)]
        net.film_providers[self.id] = LinerFilm(self, panels)

    def summarize(self, element_diag: dict[str, dict[str, Any]]) -> dict[str, Any]:
        """The liner's own diagnostics from its stations' (``__element_diag__``)."""
        n = int(self.n_segments)
        panels = [element_diag[self.panel_id(i)] for i in range(1, n + 1)]
        m_in = float(element_diag[self.segment_id(0)]["m_dot"])
        m_out = float(element_diag[self.segment_id(n)]["m_dot"])
        bleed = [float(p["m_dot"]) for p in panels]
        hole_area = math.pi * self.hole_diameter**2 / 4.0
        psi = [self.duct_area / (p["n_holes"] * hole_area) for p in panels]
        out: dict[str, Any] = {
            "m_dot": m_in,
            "m_dot_out": m_out,
            "m_dot_bleed": float(sum(bleed)),
            "bleed_fraction": float(sum(bleed) / m_in) if m_in else 0.0,
            "n_segments": float(n),
            "n_holes": float(sum(p["n_holes"] for p in panels)),
            "U1_over_Vi_max": float(max(p.get("U1_over_Vi", 0.0) for p in panels)),
            "crossflow_cd_beyond_8pct": float(
                any(p.get("crossflow_cd_beyond_8pct", 0.0) for p in panels)
            ),
            "crossflow_cd_degraded": float(
                any(p.get("crossflow_cd_degraded", 0.0) for p in panels)
            ),
            # McGreehan-Schotsch outside its source's range on any panel, by
            # variable (Re below 1e4: Cd held; l/d beyond 10): information,
            # never a warning.
            "Cd_extrapolated_Re": float(any(p.get("Cd_extrapolated_Re", 0.0) for p in panels)),
            "Cd_extrapolated_l_d": float(any(p.get("Cd_extrapolated_l_d", 0.0) for p in panels)),
            "bleed_nonuniformity": _max_over_min(bleed),
            "station_psi_min": float(min(psi)),
            "station_psi_beyond_bassett": float(min(psi) > BASSETT_PSI_MEASURED_MAX),
            "coolant_ht_plenum_assumed": 1.0,
            "is_ingesting": float(any(p.get("is_ingesting", 0.0) for p in panels)),
            # #998: any panel with no cooling path (the solver refuses such a
            # root), and any panel face above the material's limit (NaN when
            # a panel's wall temperature is undefined).
            "no_cooling_path": float(any(p.get("no_cooling_path", 0.0) for p in panels)),
        }
        if all("wall_above_limit" in p for p in panels):
            above = [float(p["wall_above_limit"]) for p in panels]
            out["T_wall_limit"] = float(panels[0]["T_wall_limit"])
            out["wall_above_limit"] = (
                1.0
                if any(a == 1.0 for a in above)
                else (0.0 if all(math.isfinite(a) for a in above) else math.nan)
            )
        if all("T_wall_hot" in p for p in panels):
            out["T_wall_hot_max"] = float(max(p["T_wall_hot"] for p in panels))
            out["rows_T_wall_hot"] = [float(p["T_wall_hot"]) for p in panels]
        out["rows_m_dot_bleed"] = bleed
        out["rows_U1_over_Vi"] = [float(p.get("U1_over_Vi", 0.0)) for p in panels]
        out["rows_Cd"] = [float(p["Cd"]) for p in panels]
        out["rows_P_backside"] = [float(p["P_in"]) for p in panels]
        return out


class LinerFilm:
    """An EffusionLiner's panels as ONE film provider (#695): what a
    ThermalWall's ``film=<liner id>`` reads. Built by ``EffusionLiner.add_to``.

    Sources and groups run upstream first ON THE GAS SIDE (``film_direction``).
    Each panel is a group of ``round(panel_length / pitch_x)`` rows, ``pitch_x``
    apart across panel boundaries too; x is measured downstream of the LAST
    row of the most downstream panel.
    """

    def __init__(self, liner: EffusionLiner, panels: list[EffusionPlateElement]) -> None:
        self.id = liner.id
        self._liner = liner
        self._panels = panels if liner.film_direction == "with_coolant" else panels[::-1]

    def film_sources(self) -> list[FilmSource]:
        return [FilmSource(p.id, p.from_node, float(p.total_hole_area)) for p in self._panels]

    def film_discharge_node(self) -> str | None:
        """The gas node, when every panel discharges into the same one."""
        nodes = {p.to_node for p in self._panels}
        return nodes.pop() if len(nodes) == 1 else None

    def film_provider_error(self) -> str | None:
        return None

    def rows_per_panel(self) -> int:
        lv = self._liner
        return max(1, round(lv.length / lv.n_segments / lv.pitch_x))

    def film_span_average_groups(
        self, x_start: float, x_end: float, M: list[float], density_ratio: list[float]
    ) -> Any:
        """The film over [x_start, x_end] [m] downstream of the last row."""
        lv = self._liner
        D = lv.hole_diameter
        return cb.film_span_average_baldauf_2002_groups(
            [self.rows_per_panel()] * len(self._panels),
            lv.pitch_x / D,
            x_start / D,
            x_end / D,
            list(M),
            list(density_ratio),
            lv.angle_deg,
            lv.pitch_y / D,
            lv.turbulence_intensity,
        )
