import FlowEdge from "./FlowEdge";
import AreaChangeNode from "./nodes/AreaChangeNode";
import ChannelNode from "./nodes/ChannelNode";
import CombustorNode from "./nodes/CombustorNode.tsx";
import DiscreteLossNode from "./nodes/DiscreteLossNode";
import EffusionLinerNode from "./nodes/EffusionLinerNode";
import EffusionPlateNode from "./nodes/EffusionPlateNode";
import EjectorNode from "./nodes/EjectorNode";
import FilmCoolingNode from "./nodes/FilmCoolingNode";
import FlowSplitterNode from "./nodes/FlowSplitterNode";
import FuelBoundaryNode from "./nodes/FuelBoundaryNode";
import HeatSourceNode from "./nodes/HeatSourceNode";
import ImpingementArrayNode from "./nodes/ImpingementArrayNode";
import LosslessNode from "./nodes/LosslessNode";
import MassBoundaryNode from "./nodes/MassBoundaryNode";
import MomentumChamberNode from "./nodes/MomentumChamberNode.tsx";
import MPCETeeNode from "./nodes/MPCETeeNode";
import OrificeNode from "./nodes/OrificeNode";
import PlenumNode from "./nodes/PlenumNode";
import PressureBoundaryNode from "./nodes/PressureBoundaryNode";
import ProbeNode from "./nodes/ProbeNode";
import SlotFilmNode from "./nodes/SlotFilmNode";
import TurbomachineNode from "./nodes/TurbomachineNode";
import VortexNode from "./nodes/VortexNode";
import WallNode from "./nodes/WallNode";
import ProbeLinkEdge from "./ProbeLinkEdge";
import ThermalEdge from "./ThermalEdge";

export const nodeTypes = {
	plenum: PlenumNode,
	mass_boundary: MassBoundaryNode,
	fuel_boundary: FuelBoundaryNode,
	pressure_boundary: PressureBoundaryNode,
	channel: ChannelNode,
	orifice: OrificeNode,
	compressor: TurbomachineNode,
	turbine: TurbomachineNode,
	impingement_array: ImpingementArrayNode,
	effusion_plate: EffusionPlateNode,
	film_cooling: FilmCoolingNode,
	slot_film: SlotFilmNode,
	effusion_liner: EffusionLinerNode,
	combustor: CombustorNode,
	momentum_chamber: MomentumChamberNode,
	discrete_loss: DiscreteLossNode,
	lossless_connection: LosslessNode,
	probe: ProbeNode,
	area_change: AreaChangeNode,
	mpce_tee: MPCETeeNode,
	flow_splitter: FlowSplitterNode,
	ejector: EjectorNode,
	vortex: VortexNode,
	wall: WallNode,
	heat_source: HeatSourceNode,
};

export const edgeTypes = {
	thermal: ThermalEdge,
	default: FlowEdge,
	probe_link: ProbeLinkEdge,
};

// A hole node's discharge coefficient as ONE object, `data.hole_cd`: the
// backend's HoleCdData, i.e. the library's HoleCd (#1098). A field left at
// null or "auto" is the node type's own default; utils/holeCd.ts holds which
// choices each node type takes and what its defaults resolve to.
export type HoleCdCorrelation =
	| "fixed"
	| "ReaderHarrisGallagher"
	| "IdelchikThick"
	| "IdelchikBeveled"
	| "IdelchikRounded"
	| "McGreehanSchotsch"
	| "Lichtarowicz";
export type HoleCdInclination = "auto" | "none" | "idelchik_3_2";
export type HoleCdInletCrossflow = "auto" | "none" | "mcgreehan_schotsch";
export type HoleCdExitCrossflow =
	| "auto"
	| "none"
	| "rowbury_2001"
	| "rogers_hersh_1976"
	| "florschuetz_metzger_su_1982";
// The duct-fed hole's drive (#1087): "auto" = "total" (upstream total to
// downstream static); "static" the pre-0.12 drive. Carried through the file,
// no editor select yet (#1124).
export type HoleCdSupplyDrive = "auto" | "static" | "total";

export interface HoleCdData {
	correlation?: HoleCdCorrelation | null;
	// The 'fixed' correlation's value and nothing else (#872).
	Cd?: number | null;
	inclination?: HoleCdInclination;
	inlet_crossflow?: HoleCdInletCrossflow;
	exit_crossflow?: HoleCdExitCrossflow;
	supply_drive?: HoleCdSupplyDrive;
}

// The node types built as a hole element, each with a `hole_cd`.
export type HoleCdNodeType =
	| "orifice"
	| "effusion_plate"
	| "film_cooling"
	| "slot_film"
	| "impingement_array"
	| "effusion_liner";
