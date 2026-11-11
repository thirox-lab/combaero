"""The hole-Cd chain's C++ composition (hole_cd_compose) against the Python
composition it replaced.

_ORACLE_* below is that Python composition, kept verbatim here as the oracle
only (python/combaero/network/hole_cd.py before the port). Every Cd and dCd
the network asks of a hole element over a grid of networks and states -- the
base correlations (fixed, ISO metering, discharge-hole with inclination and a
duct-fed U1/Vi), the exit-crossflow term, and the impingement laws suspended
or not -- is evaluated both ways. Values are identical. The partials are
identical too, except where the exit-crossflow term acts: there the chain
rule's association moved from Python into C++ (dCd/dCd0 dCd0/dRe is formed
before dRe/dm is applied), which changes the last bits only.
"""

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero import _solver_tools
from combaero.network import NetworkSolver
from combaero.network.components import (
    FilmCoolingElement,
    LosslessConnectionElement,
    MassFlowBoundary,
    MomentumChamberNode,
    OrificeElement,
    PressureBoundary,
    SlotFilmElement,
)
from combaero.network.graph import FlowNetwork
from combaero.network.hole_cd import (
    _DISCHARGE_SELECTORS,
    _METERING_SELECTORS,
    _HoleCdChain,
)
from combaero.network.impingement_array import ImpingementArray

sys.path.insert(0, str(Path(__file__).resolve().parent))
import test_effusion_liner as tl  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))

# Rounding of the reassociated chain rule where the exit term acts.
_EXIT_REL = 1e-13


# --- the oracle: the Python composition before the port --------------------------


def _oracle_base_cd(ch, state_in, state_out, flows):  # noqa: ANN001, ANN202
    el = ch._el
    if not el.use_correlation or el._orifice_geom is None:
        return max(1e-4, min(el._seed_Cd, 1.0))
    D_up = el._orifice_geom.D
    if D_up > 0.0:
        ts = cb.transport_state(state_in.Tt, state_in.P, state_in.X)
        mu = ts.mu if ts.mu > 0.0 else 1.8e-5
        m_dot_abs = abs(state_in.m_dot)
        Re_D = (4.0 * m_dot_abs) / (math.pi * D_up * mu) if m_dot_abs > 1e-12 else 1e4
    else:
        Re_D = 1e5
    Re_hole = el._hole_reynolds(state_in)
    metering = _METERING_SELECTORS.get(el.correlation)
    if metering is not None:
        return float(cb.metering_cd_and_dRe(metering, el._orifice_geom, Re_D)[0])
    u1_vi, _ = ch.supply_crossflow(state_in, state_out, flows)
    return float(
        cb.discharge_cd(
            _DISCHARGE_SELECTORS[el.correlation],
            ch.discharge_hole(),
            cb.DischargeHoleState(Re=Re_hole, U1_over_Vi=u1_vi),
        )
    )


def _oracle_metering_dcd(ch, state_in):  # noqa: ANN001, ANN202
    el = ch._el
    D_up = el._orifice_geom.D
    m = float(state_in.m_dot)
    if not (D_up > 0.0) or abs(m) <= 1e-12:
        return {}
    mu = cb.transport_state(state_in.Tt, state_in.P, state_in.X).mu
    mu_fallback = not mu > 0.0
    mu = 1.8e-5 if mu_fallback else mu
    Re = 4.0 * abs(m) / (math.pi * D_up * mu)
    metering = _METERING_SELECTORS.get(el.correlation)
    dCd_dRe = float(cb.metering_cd_and_dRe(metering, el._orifice_geom, Re)[1])
    out = {f"{el.id}.m_dot": dCd_dRe * Re / abs(m) * (1.0 if m > 0 else -1.0)}
    if not mu_fallback:
        for name, d in ch._re_viscosity_dnames(state_in, Re).items():
            out[name] = out.get(name, 0.0) + dCd_dRe * d
    return out


def _oracle_base_dcd(ch, state_in, state_out, flows, continue_re_floor):  # noqa: ANN001, ANN202
    el = ch._el
    if el._orifice_geom is None or not el.use_correlation:
        return {}
    selector = _DISCHARGE_SELECTORS.get(el.correlation)
    if selector is None:
        return _oracle_metering_dcd(ch, state_in)
    u1_vi, du = ch.supply_crossflow(state_in, state_out, flows)
    Re = el._hole_reynolds(state_in)
    _, dCd_dRe, dCd_dU = cb.discharge_cd_and_derivatives(
        selector, ch.discharge_hole(), cb.DischargeHoleState(Re=Re, U1_over_Vi=u1_vi)
    )
    if not continue_re_floor and Re < cb.discharge_cd_re_floor(selector):
        dCd_dRe = 0.0
    out: dict[str, float] = {}
    m = float(state_in.m_dot)
    if abs(m) > 1e-12:
        out[f"{el.id}.m_dot"] = float(dCd_dRe) * Re / abs(m) * (1.0 if m > 0 else -1.0)
        if el._orifice_geom.d > 0.0:
            for name, d in ch._re_viscosity_dnames(state_in, Re).items():
                out[name] = out.get(name, 0.0) + float(dCd_dRe) * d
    for name, d in du.items():
        out[name] = out.get(name, 0.0) + float(dCd_dU) * d
    return out


def _oracle_exit_crossflow(ch, cd0, state_in, state_out):  # noqa: ANN001, ANN202
    el = ch._el
    spec = ch.exit_xf
    q, dq = ch.mainstream_momentum_flux(state_out)
    if spec is None or not q > 0.0:
        return cd0, None, {}
    ex = cb.rowbury_2001_exit_density(state_in.Pt, state_in.Tt, state_out.P, state_in.X)
    rho_h = float(ex.rho)
    m = float(state_in.m_dot)
    a_exit = el.area * spec.exit_area_ratio
    kp = cb.rowbury_2001_k_perp(abs(m) / a_exit, rho_h, q, spec.theta_deg, spec.alpha_deg)
    r = cb.rowbury_2001_exit_crossflow_cd(
        cd0, kp.k, spec.coefficient_set, spec.theta_deg, spec.alpha_deg
    )
    dk: dict[str, float] = {}
    if m != 0.0:
        dk[f"{el.id}.m_dot"] = kp.dk_dG * (1.0 if m > 0.0 else -1.0) / a_exit
    dk_drho = float(kp.dk_drho_h)
    for key, d in (
        (f"{el.from_node}.Pt", ex.drho_dPt),
        (f"{el.from_node}.Tt", ex.drho_dTt),
        (f"{el.to_node}.P", ex.drho_dP),
    ):
        dk[key] = dk.get(key, 0.0) + dk_drho * float(d)
    for k, d in enumerate(ex.drho_dY):
        if d != 0.0:
            key = f"{el.from_node}.Y[{k}]"
            dk[key] = dk.get(key, 0.0) + dk_drho * float(d)
    for name, d in dq.items():
        dk[name] = dk.get(name, 0.0) + float(kp.dk_dq_m) * d
    return float(r.Cd), r, dk


def _oracle_hole_cd(ch, state_in, state_out, flows):  # noqa: ANN001, ANN202
    cd0 = _oracle_base_cd(ch, state_in, state_out, flows)
    if ch.exit_xf is None or ch._el._orifice_geom is None:
        return cd0
    return _oracle_exit_crossflow(ch, cd0, state_in, state_out)[0]


def _oracle_hole_dcd(ch, state_in, state_out, flows, continue_re_floor):  # noqa: ANN001, ANN202
    d0 = _oracle_base_dcd(ch, state_in, state_out, flows, continue_re_floor)
    if ch.exit_xf is None or ch._el._orifice_geom is None:
        return d0
    cd0 = _oracle_base_cd(ch, state_in, state_out, flows)
    _, r, dk = _oracle_exit_crossflow(ch, cd0, state_in, state_out)
    if r is None:
        return d0
    out = {name: float(r.dCd_dCd0) * v for name, v in d0.items()}
    if r.dCd_dk != 0.0:
        for name, v in dk.items():
            out[name] = out.get(name, 0.0) + float(r.dCd_dk) * v
    return out


def _oracle_jet_crossflow_cd(ch, state_in, state_out, flows, base_cd):  # noqa: ANN001, ANN202
    if ch.jet_crossflow == "none" or ch.jet_crossflow_suspended:
        return base_cd, None
    el = ch._el
    forward = state_in.Pt >= state_out.P
    up, down = (state_in, state_out) if forward else (state_out, state_in)
    up_id, down_id = (el.from_node, el.to_node) if forward else (el.to_node, el.from_node)
    dP = float(up.Pt) - float(down.P)
    rho, drho_dT, drho_dP = _solver_tools.density_and_jacobians(
        float(up.Tt), float(up.P), list(up.X)
    )
    G_i = math.sqrt(2.0 * rho * dP) if dP > 0.0 else 0.0
    flows = flows or {}
    m_c = float(sum(flows.get(eid, 0.0) for eid in el._crossflow_sources))
    G_c = abs(m_c) / el.channel_area
    if ch.jet_crossflow == "rogers_hersh_1976":
        r = cb.rogers_hersh_1976_cd(base_cd, G_i, G_c, not forward)
    else:
        r = cb.florschuetz_metzger_su_1982_cd(el.xn_d, el.yn_d, el.z_d, G_i, G_c)
    if r.capped and ch.jet_crossflow == "rogers_hersh_1976":
        return base_cd, None
    out: dict[str, float] = {}
    if r.capped or not G_i > 0.0:
        return float(r.Cd), out

    def add(name: str, v: float) -> None:
        if v != 0.0:
            out[name] = out.get(name, 0.0) + v

    dGi_ddP = rho / G_i
    dGi_drho = dP / G_i
    dcd_dGi = float(r.dCd_dG_ideal)
    add(f"{up_id}.Pt", dcd_dGi * dGi_ddP)
    add(f"{down_id}.P", -dcd_dGi * dGi_ddP)
    add(f"{up_id}.Tt", dcd_dGi * dGi_drho * float(drho_dT))
    add(f"{up_id}.P", dcd_dGi * dGi_drho * float(drho_dP))
    W = float(cb.mwmix(up.X))
    c_rho = dcd_dGi * dGi_drho * rho
    for k in range(len(up.Y)):
        add(f"{up_id}.Y[{k}]", c_rho * (1.0 - W / float(cb.species_molar_mass(k))))
    if m_c != 0.0:
        dcd_dmc = float(r.dCd_dG_c) * (1.0 if m_c > 0.0 else -1.0) / el.channel_area
        for eid in el._crossflow_sources:
            src = el._source_elements[eid]
            names = src.unknowns()
            for local, coeff in src.flow_jac_at_node(el.to_node, list(range(len(names)))).items():
                add(names[local], dcd_dmc * coeff)
    return float(r.Cd), out


def _oracle_cd(ch, state_in, state_out, flows):  # noqa: ANN001, ANN202
    base = _oracle_hole_cd(ch, state_in, state_out, flows)
    if ch.jet_crossflow == "none":
        return base
    return _oracle_jet_crossflow_cd(ch, state_in, state_out, flows, base)[0]


def _oracle_dcd(ch, state_in, state_out, flows, continue_re_floor):  # noqa: ANN001, ANN202
    if ch.jet_crossflow != "none" and not ch.jet_crossflow_suspended:
        base = _oracle_hole_cd(ch, state_in, state_out, flows)
        _, d = _oracle_jet_crossflow_cd(ch, state_in, state_out, flows, base)
        if d is not None:
            return d
    return _oracle_hole_dcd(ch, state_in, state_out, flows, continue_re_floor)


# --- the grid ------------------------------------------------------------------------


def _chamber(element, dp: float = 0.03, m_gas: float = 2.0) -> FlowNetwork:  # noqa: ANN001
    g = FlowNetwork()
    g.add_node(MassFlowBoundary("gas_in", m_dot=m_gas, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("exit", Pt=5e5, Tt=1500.0, Y=Y_AIR))
    g.add_node(PressureBoundary("cool", Pt=5e5 * (1.0 + dp), Tt=600.0, Y=Y_AIR))
    g.add_node(MomentumChamberNode("ch", area=0.01, main_inlet="feed"))
    g.add_element(LosslessConnectionElement("feed", "gas_in", "ch"))
    g.add_element(LosslessConnectionElement("tail", "ch", "exit"))
    g.add_element(element)
    return g


def _networks() -> dict[str, FlowNetwork]:
    nets: dict[str, FlowNetwork] = {}
    for corr in ("IdelchikThick", "McGreehanSchotsch", "fixed"):
        cdkw = {"Cd": 0.7} if corr == "fixed" else {}
        for xf in ("none", "rowbury_2001"):
            for ang in (30.0, 90.0):
                for dp, m_gas in ((0.03, 2.0), (0.3, 6.0)):
                    film = FilmCoolingElement(
                        "h",
                        "cool",
                        "ch",
                        hole_diameter=1e-3,
                        holes_per_row=10,
                        lateral_pitch=3e-3,
                        wall_thickness=3e-3,
                        angle_deg=ang,
                        correlation=corr,
                        cd_exit_crossflow=xf,
                        **cdkw,
                    )
                    nets[f"film|{corr}|{xf}|{ang}|{dp}"] = _chamber(film, dp, m_gas)
    for corr in ("Lichtarowicz", "IdelchikThick"):
        nets[f"slot|{corr}"] = _chamber(
            SlotFilmElement("h", "cool", "ch", 1e-3, 0.05, 4e-3, correlation=corr)
        )
    for corr, d in (("ReaderHarrisGallagher", 0.02), ("IdelchikBeveled", 0.01)):
        g = FlowNetwork()
        g.add_node(PressureBoundary("hi", Pt=1.2e5, Tt=1600.0))
        g.add_node(PressureBoundary("ho", Pt=1.0e5, Tt=300.0))
        g.add_element(OrificeElement("o", "hi", "ho", diameter=d, correlation=corr))
        nets[f"orifice|{corr}"] = g
    for law in ("rogers_hersh_1976", "florschuetz_metzger_su_1982", "none"):
        for corr in ("fixed", "McGreehanSchotsch"):
            for susp in (False, True):
                g = FlowNetwork()
                D = 0.00254
                g.add_node(PressureBoundary("sup", Pt=3e5 * 1.03, Tt=600.0, Y=Y_AIR))
                g.add_node(PressureBoundary("ex", Pt=3e5, Tt=600.0, Y=Y_AIR))
                ImpingementArray(
                    "ia",
                    n_rows=4,
                    d_jet=D,
                    xn_d=5.0,
                    yn_d=4.0,
                    z_d=1.0,
                    span=12 * 4.0 * D,
                    plate_thickness=D,
                    correlation=corr,
                    cd_crossflow=law,
                ).add_to(g, "sup", "ex")
                if susp:
                    for el in g.elements.values():
                        if hasattr(el, "cd_crossflow"):
                            el._cd_crossflow_suspended = True
                nets[f"imp|{law}|{corr}|{susp}"] = g
    nets["liner"] = tl._net(tl._liner(4), p_in=1.5e5, p_out=1.1e5)
    return nets


_NETS = _networks()


@pytest.mark.parametrize("name", sorted(_NETS))
def test_the_composed_cd_and_partials_are_the_python_compositions(
    name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[tuple] = []
    new_cd, new_dcd = _HoleCdChain.cd, _HoleCdChain.dcd

    def cd(self, state_in, state_out, flows=None):  # noqa: ANN001, ANN202
        v = new_cd(self, state_in, state_out, flows)
        calls.append(("cd", self, v, _oracle_cd(self, state_in, state_out, flows)))
        return v

    def dcd(self, state_in, state_out, flows=None, continue_re_floor=False):  # noqa: ANN001, ANN202
        v = new_dcd(self, state_in, state_out, flows, continue_re_floor)
        ref = _oracle_dcd(self, state_in, state_out, flows, continue_re_floor)
        calls.append(("dcd", self, v, ref))
        return v

    monkeypatch.setattr(_HoleCdChain, "cd", cd)
    monkeypatch.setattr(_HoleCdChain, "dcd", dcd)
    g = _NETS[name]
    g.resolve_all_topology()
    s = NetworkSolver(g)
    x0 = np.array(s._build_x0(), dtype=float)
    rng = np.random.default_rng(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for k in range(3):
            x = x0 if k == 0 else x0 * (1.0 + 0.02 * rng.standard_normal(x0.size))
            s._residuals_and_jacobian(x)
    assert any(c[0] == "dcd" for c in calls)
    for kind, chain, got, ref in calls:
        exit_on = chain.exit_xf is not None
        if kind == "cd":
            assert got == ref
            continue
        assert set(got) == set(ref)
        for key, v in ref.items():
            if exit_on:
                assert got[key] == pytest.approx(v, rel=_EXIT_REL, abs=1e-300), key
            else:
                assert got[key] == v, key
