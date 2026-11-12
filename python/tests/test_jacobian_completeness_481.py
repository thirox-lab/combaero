"""Jacobian completeness (#481 C1, C4, C5): the analytic network Jacobian
equals differences where it used to drop a dependence.

C1  A theta-sourced PressureLossElement read the combustor's UNBURNED
    temperature (theta = T_b/T_u - 1, and the correlation's reference T) and,
    for a head loss, the burned composition (its density) -- neither relayed:
    2.6% off in d/d(air flow), 5% in d/d(fuel flow).
C4  On the [200, 5000] K clamp a node's T is flat, but its sensitivities were
    relayed in full (5% off). The clamp is gone (#626): a node heated past
    5000 K relays its own sensitivities.
C5  The wall relay took dT/dQ = 1/(cp sum m); the mixer applies Q phi(m)
    (cb.mixer_heat_per_unit_flow: 1/|m| above MIXER_HEAT_MDOT_FLOOR, fading
    below it): 2x off below the floor and of the wrong sign for a negative
    total.
#760/#801  An EffusionPlateElement's self-consistent row m - m*(x) divided
    every sensitivity by 1 - K dCd/dm with C++'s Re-floor CONTINUATION slope,
    although below McGreehan-Schotsch's floor (Re 1e4) the Cd value is held
    and its true slope is 0: each duct-fed panel row was uniformly 0.987x
    (1.26e-2) to 0.97x (2.98e-2, duct Mach 0.5) the differences.
"""

from __future__ import annotations

import numpy as np
import pytest
import test_combustor_regression as tcr
from _fd_jacobian import fd_jacobian

import combaero as cb
from combaero import _solver_tools
from combaero.network import (
    EnergyBoundary,
    FlowNetwork,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
)


def _jac_err(s: NetworkSolver, x: np.ndarray, owner: str | None = None) -> float:
    """RESIDUAL-ROW-SCALING: compare J_ij |x_j| (see _build_residual_scales).
    ``owner`` restricts it to the rows whose owner label starts with it."""
    _, J = s._residuals_and_jacobian(x)
    fd = fd_jacobian(s, x)
    cols = np.maximum(np.abs(x), 1e-12)[None, :]
    Js, fds = J.toarray() * cols, fd * cols
    scale = np.maximum(np.max(np.abs(fds), axis=1, keepdims=True), 1e-300)
    err = np.max(np.abs(Js - fds) / scale, axis=1)
    if owner is not None:
        rows = [i for i, (o, _) in enumerate(s._row_kinds(x)) if o.startswith(owner)]
        assert rows, owner
        err = err[rows]
    return float(np.max(err))


@pytest.mark.parametrize(
    "corr",
    [LinearThetaFractionLoss(k=0.01, xi0=0.02), LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05)],
    ids=["fraction", "head"],
)
def test_a_theta_sourced_loss_carries_the_unburned_state(corr: object) -> None:
    """Errors were 4.4e-4 (fraction) and 6.2e-5 (head) in these scaled units."""
    net, _ = tcr._make_effective_area_inlet_network(loss_correlation=corr)
    s = NetworkSolver(net)
    r = s.solve()
    assert r["__success__"]
    assert net.elements["loss"]._theta_source_resolved == "comb"
    assert _jac_err(s, np.array(r["__x_solution__"])) < 2e-5


def test_a_node_heated_past_the_old_clamp_relays_its_temperature_sensitivity() -> None:
    """20 kW into ~1 g/s: 14 000 K. A [200, 5000] K clamp held such a plenum
    flat (its sensitivities zeroed) and discarded the energy; now it is simply
    hot, energy conserved, and its Jacobian is complete."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    g = FlowNetwork()
    for name, Pt in (("a", 1.1e5), ("b", 1.0e5)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, 300.0, Y
        g.add_node(bnd)
    p = PlenumNode("p")
    g.add_node(p)
    g.add_node(PlenumNode("q"))
    for eid, a, b in (("o1", "a", "p"), ("o2", "p", "q"), ("o3", "q", "b")):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=0.01, correlation="fixed"))
    p.add_energy_boundary(EnergyBoundary(id="heat", Q=2e4))
    s = NetworkSolver(g)
    r = s.solve(auto_retry=False)
    assert r["__success__"]
    assert r["__T_clamped__"] == [] and r["__heat_faded__"] == []
    assert r["p.Tt"] > 1e4
    assert r["__ledger__"]["energy_closure_rel"] < 1e-12
    assert _jac_err(s, np.array(r["__x_solution__"])) < 1e-6


@pytest.mark.parametrize("m", [2.5e-7, 1e-6, 1e-3, -1e-4])
def test_the_mixer_reports_its_own_heat_sensitivity(m: float) -> None:
    """dT_mix/dQ at the heat per unit flow the mixer actually applies."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))

    def mix(Q: float) -> object:
        return _solver_tools.mixer_from_streams_and_jacobians(
            [cb.MassStream(m, 400.0, 1e5, Y)], Q=Q, fraction=0.0
        )

    h = 1e-7
    fd = (mix(1e-4 + h).Tt_mix - mix(1e-4 - h).Tt_mix) / (2 * h)
    assert mix(1e-4).dTt_mix_dQ == pytest.approx(fd, rel=1e-6)
    phi, _ = cb.mixer_heat_per_unit_flow(m)
    assert mix(1e-4).dTt_mix_dQ == pytest.approx(mix(1e-4).dTt_mix_d_delta_h * phi, rel=1e-12)


def _premixed_and_dilution_network(corr: object) -> FlowNetwork:
    """A premix plenum (air + CH4, its Y moving with the flow split) and a hot
    dilution stream both feed the combustor, so T_unburned depends on a
    stream's COMPOSITION, not only its flow and Tt (#718)."""
    from combaero.network import CombustorNode, PressureLossElement

    g = FlowNetwork()
    for name, Pt, Tt, Y in (
        ("air", 1.6e5, 600.0, cb.species.dry_air_mass()),
        ("fuel", 1.6e5, 300.0, cb.species.pure_species("CH4")),
        ("dil", 1.6e5, 800.0, cb.species.dry_air_mass()),
        ("out", 1.0e5, 300.0, cb.species.dry_air_mass()),
    ):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, Tt, list(Y)
        g.add_node(bnd)
    g.add_node(PlenumNode("premix"))
    g.add_node(CombustorNode("comb", method="complete"))
    g.add_node(PlenumNode("post"))
    for eid, a, b, d in (
        ("o_air", "air", "premix", 0.03),
        ("o_fuel", "fuel", "premix", 0.004),
        ("o_pre", "premix", "comb", 0.04),
        ("o_dil", "dil", "comb", 0.02),
        ("nozzle", "post", "out", 0.05),
    ):
        g.add_element(OrificeElement(eid, a, b, Cd=0.8, diameter=d, correlation="fixed"))
    g.add_element(PressureLossElement("loss", "comb", "post", correlation=corr))
    return g


@pytest.mark.parametrize(
    "corr",
    [LinearThetaFractionLoss(k=0.01, xi0=0.02), LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05)],
    ids=["fraction", "head"],
)
def test_the_unburned_temperature_is_the_enthalpy_mix(corr: object) -> None:
    """T_unburned is the enthalpy-mixed reactant temperature, the one the C++
    combustor context defines theta with (#718); it was a mass-weighted Tt.
    Its Jacobian carries the inflows' composition through "{node}.Y[k]"."""
    g = _premixed_and_dilution_network(corr)
    s = NetworkSolver(g)
    r = s.solve()
    assert r["__success__"]
    comb = g.nodes["comb"]
    # The enthalpy balance, independently of the mixer: the inflows' total
    # enthalpy is the mixture's at T_unburned.
    n = cb.num_species()
    m_pre, m_dil = r["o_pre.m_dot"], r["o_dil.m_dot"]
    Y_pre = [r[f"premix.Y[{k}]"] for k in range(n)]
    Y_dil = list(cb.species.dry_air_mass())
    P = r["comb.P"]

    def h(T: float, Y: list[float]) -> float:
        return cb.thermo_state(T, P, list(cb.mass_to_mole(Y))).h

    M = m_pre + m_dil
    Y_mix = [(m_pre * a + m_dil * b) / M for a, b in zip(Y_pre, Y_dil, strict=True)]
    h_in = (m_pre * h(r["premix.Tt"], Y_pre) + m_dil * h(800.0, Y_dil)) / M
    assert h(comb._T_unburned, Y_mix) == pytest.approx(h_in, rel=1e-9, abs=1e-3)
    # The mass-weighted Tt it replaced misses that balance.
    T_mass = (m_pre * r["premix.Tt"] + m_dil * 800.0) / M
    assert abs(T_mass - comb._T_unburned) > 0.3
    assert _jac_err(s, np.array(r["__x_solution__"])) < 2e-5


def test_a_duct_fed_liner_at_duct_mach_half_differentiates_the_held_cd() -> None:
    """#801: a duct-fed EffusionLiner at duct Mach 0.5, its panels straddling
    McGreehan-Schotsch's Re floor (1e4). Below it the Cd value is held, so
    the panel row's implicit-function denominator 1 - K dCd/dm is 1, not the
    C++ continuation's; with the continuation every panel row was 2.98e-2 off
    (uniformly 0.97x the differences), at every FD step.

    Measured with the fix: 2e-7 for the whole network at a 1e-6 step (6.6e-9
    on the panel rows; the rest is FD truncation in the duct rows)."""
    import test_effusion_liner as tel

    s, r = tel._solve(tel._net(tel._liner(3, duct_height=0.005), p_in=1.4e5, p_out=1.2e5))
    assert r["__element_diag__"]["ln__s0"]["M_in_duct"] == pytest.approx(0.5, abs=0.02)
    x = np.array(r["__x_solution__"])
    re = [r["__element_diag__"][f"ln__p{i}"]["Re_hole"] for i in (1, 2, 3)]
    assert min(re) < 1e4 < max(re), re
    assert _jac_err(s, x, owner="EffusionPlateElement") < 1e-6
    assert _jac_err(s, x) < 1e-6
