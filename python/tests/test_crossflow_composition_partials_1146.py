"""A crossflow segment's terms carry the composition partials of their gas (#1146).

Every momentum term of a ``CrossflowSegmentElement`` forms a density from a
gas state: the reservoir entry (``entry``, also the lossless ``_feed_K``
entry), the station halves (``from``, ``to``), the reversed reservoir exit
(``exit_in``), the exit head into a moving chamber (``exit_head``) and the
taper's Bernoulli term (``area``); a joining jet's axial momentum
(``injection_momentum``) is its supply's isentropic velocity. The C++ kernels
(``station_half_drop``, ``channel_entry_drop``, ``segment_area_change_drop``,
``jet_impulse``) returned no partial in the composition, so a gas whose Y
moves with the flows (a burning combustor, a mixing plenum, a station a mixed
jet joins) lacked d/dY in every one of those rows: 8.8e-4 scaled on the
segment row fed by a burning combustor, 7.2e-4 reversed out of one, 2.7e-5 on
the reversed chain's legs and 3.5e-5 with a jet from a mixing plenum. Each
partial is keyed to the node that supplies the gas (``gas_node``, #1137).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from _fd_jacobian import fd_jacobian, scaled_error
from test_crossflow_entry_gas_1137 import _merge

import combaero as cb
from combaero.network import (
    ChannelElement,
    CombustorNode,
    CrossflowSegmentElement,
    FlowNetwork,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    NetworkMixtureState,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))
CH4 = list(cb.species.to_mass(cb.species.pure_species("CH4")))
#: A non-burning mixed stream: 3 % CH4 by mass in air.
MIXED = [0.97 * a + 0.03 * f for a, f in zip(Y_AIR, CH4, strict=True)]
A_DUCT = 2.0e-4
P0 = 1.0e5
N_SPECIES = len(Y_AIR)
REF = int(cb.species_index_from_name("N2"))


def _reservoir(g: FlowNetwork, kind: str, p_sup: float) -> str:
    """'comb': a CombustorNode burning 2e-4 kg/s CH4 in the supply's air;
    'mixpl': a PlenumNode mixing the air with a 400 K MIXED stream. Both
    relay T0 and Y from the flows."""
    g.add_node(PressureBoundary("sup", Pt=p_sup, Tt=300.0, Y=Y_AIR))
    if kind == "comb":
        g.add_node(CombustorNode("res"))
        g.add_node(MassFlowBoundary("fuel", m_dot=2e-4, Tt=300.0, Y=CH4))
        g.add_element(LosslessConnectionElement("inj", "fuel", "res"))
    else:
        g.add_node(PlenumNode("res"))
        g.add_node(PressureBoundary("mix", Pt=p_sup, Tt=400.0, Y=MIXED))
        g.add_element(
            OrificeElement("w", "mix", "res", diameter=0.012, Cd=0.7, correlation="fixed",
                           in_pipe=False)
        )  # fmt: skip
    g.add_element(
        OrificeElement("o", "sup", "res", diameter=0.03, Cd=0.8, correlation="fixed",
                       in_pipe=False)
    )  # fmt: skip
    return "res"


def _network(regime: str, res: str, topo: str, rev: bool, dp: float = 8000.0) -> FlowNetwork:
    """A relayed-gas reservoir (``_reservoir``) feeding, by ``topo``:

    * 'entry': one segment (entry_K 0.5) into a plenum taking a 700 K MIXED
      side stream, then a channel to the exit;
    * 'chamber': the same into a merge chamber (the segment its main inlet);
    * 'chain': a stationed chain res -u0-> st -u1-> exit, u0 fed at rest (the
      lossless entry, ``_feed_K``), a 30 deg MIXED jet joining at st, u1
      tapered (the ``area`` term);
    * 'chainmc': the same chain into a MomentumChamberNode (``exit_head``);
    * 'chainjet': the chain with the jet's supply a plenum mixing air and
      MIXED (the jet's own Y relayed).

    ``rev`` swaps the supply and exit pressures."""
    p_sup, p_exit = (P0, P0 + dp) if rev else (P0 + dp, P0)
    g = FlowNetwork()
    g.add_node(PressureBoundary("exit", Pt=p_exit, Tt=300.0, Y=Y_AIR))
    r = _reservoir(g, res, p_sup)
    D = math.sqrt(4.0 * A_DUCT / math.pi)
    if topo in ("entry", "chamber"):
        if topo == "entry":
            g.add_node(PlenumNode("mc"))
        else:
            g.add_node(MomentumChamberNode("mc", area=A_DUCT, regime=regime, main_inlet="s1x"))
        g.add_element(
            CrossflowSegmentElement("s1x", r, "mc", length=0.05, area=A_DUCT, entry_K=0.5,
                                    regime=regime)
        )  # fmt: skip
        g.add_node(PressureBoundary("hot", Pt=p_sup, Tt=700.0, Y=MIXED))
        g.add_element(OrificeElement("h", "hot", "mc", diameter=0.008, Cd=0.6, correlation="fixed"))
        g.add_element(
            ChannelElement("tail", "mc", "exit", length=0.05, diameter=1.3 * D, regime=regime)
        )
        return g
    g.add_node(PlenumNode("st"))
    end = "exit"
    if topo == "chainmc":
        g.add_node(MomentumChamberNode("mc", area=1.5 * A_DUCT, regime=regime))
        g.add_element(
            ChannelElement("tail", "mc", "exit", length=0.05, diameter=1.3 * D, regime=regime)
        )
        end = "mc"
    g.add_element(
        CrossflowSegmentElement("u0", r, "st", length=0.05, area=A_DUCT, to_kappa=0.0,
                                next_seg="u1", regime=regime)
    )  # fmt: skip
    g.add_element(
        CrossflowSegmentElement("u1", "st", end, length=0.05, area=A_DUCT, area_to=1.4 * A_DUCT,
                                from_kappa=0.0, prev_seg="u0", regime=regime)
    )  # fmt: skip
    p_jet = (p_exit if rev else p_sup) + 0.3 * dp
    if topo == "chainjet":
        g.add_node(PressureBoundary("jA", Pt=p_jet, Tt=300.0, Y=Y_AIR))
        g.add_node(PressureBoundary("jM", Pt=p_jet, Tt=600.0, Y=MIXED))
        g.add_node(PlenumNode("jet"))
        for eid, src in (("ja", "jA"), ("jm", "jM")):
            g.add_element(
                OrificeElement(eid, src, "jet", diameter=0.01, Cd=0.8, correlation="fixed",
                               in_pipe=False)
            )  # fmt: skip
    else:
        g.add_node(PressureBoundary("jet", Pt=p_jet, Tt=600.0, Y=MIXED))
    g.add_element(
        OrificeElement("j", "jet", "st", diameter=0.005, Cd=0.7, correlation="fixed",
                       injection_angle_deg=30.0)
    )  # fmt: skip
    return g


# ---- element level -----------------------------------------------------------


def _mixed_Y(seed: int) -> list[float]:
    """A composition with every species present (no FD step crosses Y = 0,
    where the state's clamp and s() drop a species)."""
    rng = np.random.default_rng(seed)
    y = 0.8 * np.array(MIXED) + 0.2 * rng.uniform(0.2, 1.0, N_SPECIES) / N_SPECIES
    return list(y / y.sum())


def _state(P: float, Tt: float, m: float, Y: list[float], inc: bool) -> NetworkMixtureState:
    return NetworkMixtureState(P, P, Tt=Tt, m_dot=m, Y=Y, incompressible=inc)


def _resolved(regime: str, topo: str, seg: str) -> CrossflowSegmentElement:
    g = _network(regime, "mixpl", topo, rev=False)
    for el in g.elements.values():
        el.resolve_topology(g)
    return g.elements[seg]


def _term_dY_errors(
    seg: CrossflowSegmentElement,
    kind: str,
    states: dict[str, NetworkMixtureState],
    flows: dict[str, float],
) -> tuple[float, float]:
    """(error, scale) of one momentum term's Y partials of BOTH end nodes
    against central differences of its dP, along e_k - e_ref (sum Y = 1)."""
    ids = (seg.from_node, seg.to_node)

    def term(st: dict[str, NetworkMixtureState]) -> tuple:
        terms = seg._momentum_terms(st[ids[0]], st[ids[1]], flows)
        found = [t for t in terms if t[5] == kind]
        assert len(found) == 1, (kind, [t[5] for t in terms])
        return found[0]

    partials = seg._term_partials(term(states))
    h = 1e-7
    err = scale = 0.0
    for nid in ids:
        base = states[nid]
        for k in range(N_SPECIES):
            if k == REF:
                continue
            vals = []
            for sgn in (1.0, -1.0):
                y = list(base.Y)
                y[k] += sgn * h
                y[REF] -= sgn * h
                st = dict(states)
                st[nid] = _state(base.P, base.Tt, base.m_dot, y, base.incompressible)
                vals.append(term(st)[0].dP)
            fd = (vals[0] - vals[1]) / (2.0 * h)
            an = partials.get(f"{nid}.Y[{k}]", 0.0) - partials.get(f"{nid}.Y[{REF}]", 0.0)
            err = max(err, abs(an - fd))
            scale = max(scale, abs(fd))
    return err, scale


# (regime, topology, segment, term kind, flow sign)
TERM_CASES = [
    ("incompressible", "entry", "s1x", "entry", 1.0),
    ("compressible", "entry", "s1x", "entry", 1.0),
    ("incompressible", "entry", "s1x", "exit_in", -1.0),
    ("incompressible", "chamber", "s1x", "exit_head", 1.0),
    ("incompressible", "chamber", "s1x", "exit_head", -1.0),
    ("incompressible", "chain", "u0", "entry", 1.0),
    ("incompressible", "chain", "u0", "to", 1.0),
    ("incompressible", "chain", "u0", "to", -1.0),
    ("compressible", "chain", "u0", "to", 1.0),
    ("compressible", "chain", "u0", "to", -1.0),
    ("incompressible", "chain", "u1", "from", 1.0),
    ("incompressible", "chain", "u1", "from", -1.0),
    ("compressible", "chain", "u1", "from", 1.0),
    ("compressible", "chain", "u1", "from", -1.0),
    ("incompressible", "chain", "u1", "area", 1.0),
    ("incompressible", "chain", "u1", "area", -1.0),
]


@pytest.mark.parametrize(("regime", "topo", "seg_id", "kind", "sign"), TERM_CASES)
def test_term_composition_partials_match_fd(
    regime: str, topo: str, seg_id: str, kind: str, sign: float
) -> None:
    """Each term's d/dY, keyed to the node whose gas it reads and to no other,
    against central differences in both end nodes' compositions. Every end
    node holds its own gas, so a partial on the wrong node fails both."""
    seg = _resolved(regime, topo, seg_id)
    inc = regime == "incompressible"
    m = sign * 0.012
    states = {
        seg.from_node: _state(1.06e5, 650.0, m, _mixed_Y(1), inc),
        seg.to_node: _state(1.0e5, 1400.0, m, _mixed_Y(2), inc),
    }
    flows = {"u0": 0.9 * m, "u1": 1.1 * m, "s1x": m}
    err, scale = _term_dY_errors(seg, kind, states, flows)
    assert scale > 0.0
    assert err < 1e-7 * scale, (err, scale)


@pytest.mark.parametrize("Pt", [1.2e5, 1.9e5, 3.0e5])
def test_jet_momentum_composition_partials_match_fd(Pt: float) -> None:
    """A 30 deg jet's axial momentum carries its supply's d/dY, unchoked and
    choked (``jet_impulse``)."""
    g = _network("compressible", "mixpl", "chainjet", rev=False)
    for el in g.elements.values():
        el.resolve_topology(g)
    jet = g.elements["j"]
    st_out = _state(1.0e5, 600.0, 0.0, _mixed_Y(3), False)

    def J(y: list[float]) -> tuple[float, dict[str, float]]:
        st_in = NetworkMixtureState(Pt, Pt, Tt=700.0, m_dot=0.02, Y=y)
        return jet.injection_momentum(st_in, st_out)

    y0 = _mixed_Y(4)
    _, d = J(y0)
    h = 1e-7
    err = scale = 0.0
    for k in range(N_SPECIES):
        if k == REF:
            continue
        yp, ym = list(y0), list(y0)
        yp[k] += h
        yp[REF] -= h
        ym[k] -= h
        ym[REF] += h
        fd = (J(yp)[0] - J(ym)[0]) / (2.0 * h)
        an = d.get(f"jet.Y[{k}]", 0.0) - d.get(f"jet.Y[{REF}]", 0.0)
        err = max(err, abs(an - fd))
        scale = max(scale, abs(fd))
    assert scale > 0.0
    assert err < 1e-6 * scale, (err, scale)


# ---- network -----------------------------------------------------------------


def _all_rows_fd_error(build, own: str, point: str) -> float:
    """Worst scaled FD error (J_ij |x_j|) over EVERY row, fresh solver at the
    solved root, perturbed, or with ``own`` at m = +-1e-7."""
    s = NetworkSolver(build())
    r = s.solve()
    assert r["__success__"], r.get("__message__")
    x = np.array(r["__x_solution__"], dtype=float)
    if point == "perturbed":
        x = x * (1.0 + 1e-3 * np.random.default_rng(7).uniform(-1.0, 1.0, x.size))
    elif point in ("m+", "m-"):
        x[s.unknown_names.index(own)] = 1e-7 if point == "m+" else -1e-7
    f = NetworkSolver(build())
    f._build_x0()
    J = f._residuals_and_jacobian(x)[1]
    fd = fd_jacobian(f, x, floor=1e-10)
    return float(np.max(scaled_error(J, fd, x)))


NET_CASES = [
    (res, topo) for res in ("comb", "mixpl") for topo in ("entry", "chamber", "chain", "chainmc", "chainjet")
]  # fmt: skip


@pytest.mark.parametrize("point", ["root", "perturbed", "m+", "m-"])
@pytest.mark.parametrize("rev", [False, True])
@pytest.mark.parametrize(("res", "topo"), NET_CASES)
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_relayed_gas_reservoir_jacobian(
    regime: str, res: str, topo: str, rev: bool, point: str
) -> None:
    """A burning combustor or a mixing plenum feeding a segment, a merge
    chamber or a chain (and reversed: the chain's gas from the exit and the
    mixed jet): every row to 1e-7 scaled. Measured <= 2.8e-8 (main: up to
    1.0e-3 on the segment rows)."""
    own = "s1x.m_dot" if topo in ("entry", "chamber") else "u0.m_dot"
    err = _all_rows_fd_error(lambda: _network(regime, res, topo, rev), own, point)
    assert err < 1e-7, err


@pytest.mark.parametrize("point", ["root", "perturbed", "m+", "m-"])
@pytest.mark.parametrize("reservoir", ["fed", "mix", "chain", "mcres"])
@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_reversed_out_of_a_burning_combustor_jacobian(
    regime: str, reservoir: str, point: str
) -> None:
    """Reversed flow out of a burning CombustorNode through an entry segment
    (#1137's networks): incompressible the reversed exit term ``exit_in``
    reads the burnt gas; behind the 'chain' reservoir the burnt gas runs back
    through the station, both legs' halves at its density. Every row to 1e-7
    scaled. Measured <= 1.3e-8 (main: 7.2e-4 incompressible, 2.7e-5 on the
    compressible chain's legs)."""
    err = _all_rows_fd_error(
        lambda: _merge(regime, 8000.0, "combustor", rev=True, reservoir=reservoir, fuel=2e-4),
        "s1x.m_dot",
        point,
    )
    assert err < 1e-7, err
