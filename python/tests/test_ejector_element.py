"""EjectorElement: 3-port supersonic ejector on the momentum-CV junction
topology.

Physics is Huang (1999) Eqs. 1-8 (entrainment ratio) + Kracik & Dvorak
(2016) Eqs. 7-13 (critical back pressure) -- see
combaero.network._ejector_huang1999 and validation/ejector/README.md for
the derivation, paper references, and accuracy comparison against
alternatives. This file tests the ELEMENT (residuals/diagnostics/topology/
solver integration), not the underlying physics functions, which are
already covered by python/tests/test_ejector_huang.py against Huang's
digitized Table 3/4 data.

References: sources/ejector/Huang_1d_analysis_ejector.pdf,
sources/ejector/Development_of_an_Analytical_Method_for.pdf
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
import combaero.species as sp
from combaero.network import ChannelElement, FlowNetwork, MomentumChamberNode, PressureBoundary
from combaero.network.components import NetworkMixtureState
from combaero.network.ejector_element import EjectorElement, choke_plane_gamma
from combaero.network.solver import NetworkSolver

_DRY_AIR_Y = list(cb.mole_to_mass(sp.dry_air()))

_THROAT_AREA = 3.14e-5  # ~6.3 mm diameter
_NOZZLE_EXIT_AREA = 1.0e-4
_MIXING_AREA = 8.0e-4
_CHANNEL_D = 0.05


def _states() -> list[NetworkMixtureState]:
    primary = NetworkMixtureState(P=700000.0, Pt=700000.0, Tt=440.0, m_dot=0.0, Y=_DRY_AIR_Y)
    secondary = NetworkMixtureState(P=22850.0, Pt=22850.0, Tt=280.0, m_dot=0.0, Y=_DRY_AIR_Y)
    outlet = NetworkMixtureState(P=40000.0, Pt=40000.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
    return [primary, secondary, outlet]


def _element(**overrides) -> EjectorElement:
    kwargs = {
        "id": "ej1",
        "primary_node": "p",
        "secondary_node": "s",
        "outlet_node": "o",
        "throat_area": _THROAT_AREA,
        "nozzle_exit_area": _NOZZLE_EXIT_AREA,
        "mixing_area": _MIXING_AREA,
        # The compressible closures are the subject (#588: the library
        # default is incompressible until #585).
        "regime": "compressible",
        # K&D's anchors: these fixtures place the outlet inside its critical window;
        # huang1999 (the default) has its own tests (test_ejector_huang_anchors.py,
        # test_ejector_anchor_closure_element.py, test_ejector_subcritical_element.py).
        "anchor_closure": "kracik_dvorak",
    }
    kwargs.update(overrides)
    return EjectorElement(**kwargs)


# ---------------------------------------------------------------------------
# Construction / validation
# ---------------------------------------------------------------------------


def test_construction_sets_area_ratios():
    e = _element()
    assert e.area_ratio_nozzle == pytest.approx(_NOZZLE_EXIT_AREA / _THROAT_AREA)
    assert e.area_ratio_mix == pytest.approx(_MIXING_AREA / _THROAT_AREA)
    assert e.geometry.area_ratio_nozzle == e.area_ratio_nozzle
    assert e.geometry.area_ratio_mix == e.area_ratio_mix


def test_construction_rejects_nonpositive_throat_area():
    with pytest.raises(ValueError, match="throat_area"):
        _element(throat_area=0.0)


def test_construction_rejects_nozzle_exit_not_exceeding_throat():
    with pytest.raises(ValueError, match="nozzle_exit_area"):
        _element(nozzle_exit_area=_THROAT_AREA * 0.5)


def test_construction_rejects_mixing_area_not_exceeding_nozzle_exit():
    with pytest.raises(ValueError, match="mixing_area"):
        _element(mixing_area=_NOZZLE_EXIT_AREA * 0.5)


def test_construction_rejects_nonpositive_recovery_efficiency():
    with pytest.raises(ValueError, match="recovery_efficiency"):
        _element(recovery_efficiency=0.0)


def test_ports_are_primary_then_secondary_then_outlet():
    e = _element()
    assert e.inlet_nodes == ["p", "s"]
    assert e.outlet_nodes == ["o"]
    assert e.port_nodes == ["p", "s", "o"]
    assert e._port_signs == [-1.0, -1.0, 1.0]


def test_port_areas_derived_from_geometry():
    # The ejector supplies its own port-face areas from geometry (ordered
    # primary, secondary, outlet) so the base class never needs to infer them
    # from a connecting channel -- an ejector can be fed by a bare boundary.
    e = _element()
    assert e.port_areas == [
        _NOZZLE_EXIT_AREA,
        _MIXING_AREA - _NOZZLE_EXIT_AREA,
        _MIXING_AREA,
    ]


def test_row_scale_kinds_is_mdot_mdot_mdot_p():
    # Opposite of the base class's own N-pressure-rows + 1-mass-row pattern
    # -- see solver.py's scaling block and the module docstring.
    assert _element().row_scale_kinds() == ["mdot", "mdot", "mdot", "p"]


# ---------------------------------------------------------------------------
# Residuals: unit-level, against a hand-guessed (not necessarily converged) state
# ---------------------------------------------------------------------------


def test_residuals_cross_check_against_reference_physics():
    """residuals() must reproduce the SAME numbers as calling the reference
    physics functions directly with the live-evaluated gamma/R -- this is
    the cross-check against combaero.network._ejector_huang1999 (which is
    itself validated against Huang's Table 3/4 in test_ejector_huang.py)."""
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()

    mp_guess, ms_guess = 0.5, 0.3
    mdot_out_guess = mp_guess + ms_guess
    port_mdots = [-mp_guess, -ms_guess, mdot_out_guess]
    P_jct_guess = 100000.0

    res, jac = e.residuals([primary, secondary, outlet], P_jct_guess, port_mdots)
    assert len(res) == 4
    # All four rows are now analytic (full C++ (f,J) path); see
    # test_residuals_jacobian_matches_finite_difference for the accuracy check.
    assert set(jac.keys()) == {0, 1, 2, 3}

    gamma = choke_plane_gamma(secondary.Tt, secondary.X)
    r_gas = cb.specific_gas_constant(secondary.X)
    geom = e.geometry
    mp_choked = choked_mass_flow(primary.Pt, primary.Tt, e.throat_area, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma).omega
    pc = critical_back_pressure(
        primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma, r_gas
    ).p_c_pa

    assert res[0] == pytest.approx(mp_guess - mp_choked)
    assert res[1] == pytest.approx(ms_guess - omega * mp_guess)
    assert res[2] == pytest.approx(0.0, abs=1e-9)  # mdot_out_guess == mp+ms by construction
    assert res[3] == pytest.approx(P_jct_guess - pc)


def test_residuals_all_zero_at_the_converged_point():
    """Solve the 4 equations directly (no full network) for the point where
    all residuals vanish, then verify residuals() actually returns zero
    there -- a tighter check than the cross-check above, which only proves
    internal consistency, not that R=0 has the physically expected root."""
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()

    gamma = choke_plane_gamma(secondary.Tt, secondary.X)
    r_gas = cb.specific_gas_constant(secondary.X)
    geom = e.geometry
    mp = choked_mass_flow(primary.Pt, primary.Tt, e.throat_area, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma).omega
    ms = omega * mp
    mdot_out = mp + ms
    pc = critical_back_pressure(
        primary.Pt, primary.Tt, secondary.Pt, secondary.Tt, geom, gamma, r_gas
    ).p_c_pa

    port_mdots = [-mp, -ms, mdot_out]
    res, _jac = e.residuals([primary, secondary, outlet], pc, port_mdots)
    for r in res:
        assert r == pytest.approx(0.0, abs=1e-6)


def test_residuals_jacobian_matches_finite_difference():
    """Central-difference-check every analytic jac-dict entry against
    perturbed residuals() calls, at a non-converged operating point. This
    directly proves the column-name / sign / relay mapping in residuals()
    is correct (the error-prone part of the C++ (f,J) rewire).

    gamma is live: the secondary-temperature column carries its
    d(gamma)/d(Tt) term (#1246; it used to be held frozen here). The
    composition keys are checked along a direction in
    test_gas_partials_chain_onto_the_secondary_composition_1246.
    """
    import dataclasses

    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()

    # A deliberately non-converged base point (residuals nonzero here, so the
    # Jacobian is exercised away from R=0).
    P_jct0 = 90000.0
    port_mdots0 = [-0.070, -0.030, 0.100]
    states0 = [primary, secondary, outlet]

    _res0, jac = e.residuals(states0, P_jct0, port_mdots0)

    def res_row(row_idx, *, P_jct=P_jct0, port_mdots=None, states=None):
        pm = list(port_mdots0) if port_mdots is None else port_mdots
        st = states0 if states is None else states
        return e.residuals(st, P_jct, pm)[0][row_idx]

    # Perturbation drivers keyed by jac column name. Each returns res_row as a
    # function of the scalar unknown behind that column.
    def perturb_mdot(port_i, sign):
        # column "chan_x.m_dot" = raw unknown u; port_mdots[port_i] = sign*u.
        def f(row_idx, u):
            pm = list(port_mdots0)
            pm[port_i] = sign * u
            return res_row(row_idx, port_mdots=pm)

        # base value of the raw unknown u = port_mdots[port_i]/sign
        return f, port_mdots0[port_i] / sign

    def perturb_state(idx, field):
        def f(row_idx, v):
            st = list(states0)
            st[idx] = dataclasses.replace(states0[idx], **{field: v})
            return res_row(row_idx, states=st)

        return f, getattr(states0[idx], field)

    def perturb_pjct(row_idx, v):
        return res_row(row_idx, P_jct=v)

    drivers = {
        "chan_p.m_dot": perturb_mdot(0, e._port_signs[0]),
        "chan_s.m_dot": perturb_mdot(1, e._port_signs[1]),
        "chan_o.m_dot": perturb_mdot(2, e._port_signs[2]),
        "p.Pt": perturb_state(0, "Pt"),
        "p.Tt": perturb_state(0, "Tt"),
        "s.Pt": perturb_state(1, "Pt"),
        "s.Tt": perturb_state(1, "Tt"),
        "ej1.P_jct": (perturb_pjct, P_jct0),
    }

    def central_diff(fn, row_idx, x0):
        eps = 1e-6 * abs(x0) if x0 != 0.0 else 1e-6
        return (fn(row_idx, x0 + eps) - fn(row_idx, x0 - eps)) / (2.0 * eps)

    checked = 0
    for row_idx, row in jac.items():
        for col, analytic in row.items():
            if ".Y[" in col:
                continue
            fn, x0 = drivers[col]
            fd = central_diff(fn, row_idx, x0)
            assert analytic == pytest.approx(fd, rel=1e-5, abs=1e-9), (
                f"row {row_idx} col {col}: analytic {analytic} vs FD {fd}"
            )
            checked += 1
    # Guard against a silently-empty jac dict masking a regression.
    assert checked >= 12


@pytest.mark.parametrize("anchor", ["kracik_dvorak", "huang1999"])
@pytest.mark.parametrize(
    "port_mdots",
    [
        [-0.44, -0.13, 0.57],  # forward, choked
        [-0.44, 0.07, 0.37],  # secondary reversed: the port gas is the ejector's
        [-0.0051, 0.010, -0.0049],  # jet pump, secondary and outlet reversed
    ],
)
def test_gas_partials_chain_onto_the_secondary_composition_1246(port_mdots, anchor):
    """#1246: gamma and R are the secondary port gas's, which moves with the
    flows when a leg reverses (and with a side stream). The rows' partials
    in them are chained onto the secondary node's T0 and composition; check
    both against residuals() differenced in that T0 and along a composition
    direction (CO2 for N2, so the mass fractions still sum to 1)."""
    import dataclasses

    e = _element(anchor_closure=anchor)
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, _, outlet = _states()
    if abs(port_mdots[0]) < 0.1:  # the jet pump: near-atmospheric supplies
        primary = NetworkMixtureState(P=105000.0, Pt=105000.0, Tt=305.0, m_dot=0.0, Y=_DRY_AIR_Y)
        outlet = NetworkMixtureState(P=100500.0, Pt=100500.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
        sec_P, P_jct = 100000.0, 99000.0
    else:
        sec_P, P_jct = 22850.0, 60000.0
    i_co2 = cb.species_index_from_name("CO2")
    i_n2 = cb.species_index_from_name("N2")
    Y = list(_DRY_AIR_Y)
    Y[i_co2] += 0.2
    Y[i_n2] -= 0.2
    secondary = NetworkMixtureState(P=sec_P, Pt=sec_P, Tt=320.0, m_dot=0.0, Y=Y)
    dY = [0.0] * len(Y)
    dY[i_co2], dY[i_n2] = 1.0, -1.0

    def rows(sec: NetworkMixtureState) -> list[float]:
        return e.residuals([primary, sec, outlet], P_jct, list(port_mdots))[0]

    res, jac = e.residuals([primary, secondary, outlet], P_jct, list(port_mdots))
    h_t, h_y = 1e-3, 1e-6
    t_p = rows(dataclasses.replace(secondary, Tt=secondary.Tt + h_t))
    t_m = rows(dataclasses.replace(secondary, Tt=secondary.Tt - h_t))
    y_p = rows(dataclasses.replace(secondary, Y=[y + h_y * d for y, d in zip(Y, dY, strict=True)]))
    y_m = rows(dataclasses.replace(secondary, Y=[y - h_y * d for y, d in zip(Y, dY, strict=True)]))
    moved = 0
    for row in range(4):
        fd_t = (t_p[row] - t_m[row]) / (2 * h_t)
        fd_y = (y_p[row] - y_m[row]) / (2 * h_y)
        an_t = jac[row].get("s.Tt", 0.0)
        an_y = sum(jac[row].get(f"s.Y[{k}]", 0.0) * d for k, d in enumerate(dY))
        # On the row's own scale: its value and these terms.
        scale = abs(res[row]) + abs(fd_t) * secondary.Tt + abs(fd_y) + 1e-12
        assert an_t * secondary.Tt == pytest.approx(fd_t * secondary.Tt, abs=1e-6 * scale), row
        assert an_y == pytest.approx(fd_y, abs=1e-6 * scale), row
        moved += abs(fd_y) > 1e-6 * scale
    # The composition moves the rows (R and gamma both), so this is not vacuous.
    assert moved >= 1


def test_diagnostics_reports_actual_omega_and_regime():
    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary, secondary, outlet = _states()
    # Actual entrainment ratio is now reported directly from the port flows
    # (regime-independent), NOT the choked-mode model value.
    port_mdots = [-0.5, -0.3, 0.8]
    diag = e.diagnostics([primary, secondary, outlet], 100000.0, port_mdots)
    assert diag["entrainment_ratio"] == pytest.approx(0.3 / 0.5)  # ms / mp
    assert diag["gamma"] == pytest.approx(1.4, abs=0.05)  # air
    assert diag["r_gas"] == pytest.approx(287.0, abs=1.0)  # air
    assert diag["m_dot_primary"] == pytest.approx(0.5)
    assert diag["m_dot_secondary"] == pytest.approx(0.3)
    assert diag["m_dot_outlet"] == pytest.approx(0.8)
    # Primary at 700 kPa is choked (s_choke -> 1) and the outlet (40 kPa) is
    # below the ejector's P_c*: the double-choked critical plateau.
    assert diag["s_choke"] == pytest.approx(1.0)
    assert diag["critical_mode"] == 1.0
    assert diag["P_c_star"] > float(outlet.Pt)


def test_diagnostics_reports_jet_pump_regime():
    # A near-atmospheric primary that cannot choke: the unchoked jet-pump
    # regime. s_choke -> 0 and the critical-only fields are suppressed.
    e = _element()
    e._port_element_ids = ["chan_p", "chan_s", "chan_o"]
    primary = NetworkMixtureState(P=101325.0, Pt=101325.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
    secondary = NetworkMixtureState(P=101325.0, Pt=101325.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
    outlet = NetworkMixtureState(P=101325.0, Pt=101325.0, Tt=300.0, m_dot=0.0, Y=_DRY_AIR_Y)
    port_mdots = [-0.0051, -0.0338, 0.0389]
    diag = e.diagnostics([primary, secondary, outlet], 100142.0, port_mdots)
    assert diag["s_choke"] == pytest.approx(0.0)
    assert diag["critical_mode"] == 0.0
    assert diag["entrainment_ratio"] == pytest.approx(0.0338 / 0.0051)
    # Critical-only diagnostics have no value here: present, NaN (a result's
    # keys never follow its regime).
    for key in ("P_c_star", "P_b0", "entrainment_ratio_critical"):
        assert math.isnan(diag[key])


def test_verify_solution_consistent_always_true_all_regimes_modeled():
    # The subcritical/jet-pump regimes are now MODELED (not demoted): the old
    # "outlet exceeds P_c* -> inconsistent" check is gone, so verification is
    # unconditionally satisfied once the residual system has converged.
    e = _element()
    assert e.verify_solution_consistent({"ej1.P_jct": 100000.0, "o.Pt": 40000.0}) is True
    assert e.verify_solution_consistent({"ej1.P_jct": 40000.0, "o.Pt": 100000.0}) is True
    assert e.verify_solution_consistent({}) is True


# ---------------------------------------------------------------------------
# choke_plane_gamma: sanity against known air properties
# ---------------------------------------------------------------------------


def test_choke_plane_gamma_matches_known_air_value():
    X = sp.dry_air()
    gamma = choke_plane_gamma(281.15, X)
    assert gamma == pytest.approx(1.4, abs=0.02)


# ---------------------------------------------------------------------------
# Full network integration: build a real graph, solve it, check the result.
# ---------------------------------------------------------------------------


def _expected_operating_point() -> dict[str, float]:
    """mp / ms / mdot_out / P_c* for the network's boundary conditions, from
    the same reference physics the element itself uses -- a good Newton
    starting point (the element now has a full analytic Jacobian, but a
    realistic seed still helps this stiff double-choked system converge fast).
    Port static P/Pt are seeded at the boundary Pt (negligible drop expected
    over the short, large-diameter test channels).
    """
    from combaero.network._ejector_huang1999 import (
        ETA_P,
        EjectorGeometry,
        choked_mass_flow,
        critical_back_pressure,
        entrainment_ratio,
    )

    geom = EjectorGeometry(
        area_ratio_nozzle=_NOZZLE_EXIT_AREA / _THROAT_AREA,
        area_ratio_mix=_MIXING_AREA / _THROAT_AREA,
    )
    gamma = choke_plane_gamma(280.0, sp.dry_air())
    r_gas = cb.specific_gas_constant(sp.dry_air())
    mp = choked_mass_flow(700000.0, 440.0, _THROAT_AREA, gamma, r_gas, eta=ETA_P)
    omega = entrainment_ratio(700000.0, 440.0, 22850.0, 280.0, geom, gamma).omega
    ms = omega * mp
    pc = critical_back_pressure(700000.0, 440.0, 22850.0, 280.0, geom, gamma, r_gas).p_c_pa
    return {"mp": mp, "ms": ms, "mout": mp + ms, "pc": pc}


def _build_network() -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(PressureBoundary("pb_primary", Pt=700000.0, Tt=440.0, Y=_DRY_AIR_Y))
    net.add_node(PressureBoundary("pb_secondary", Pt=22850.0, Tt=280.0, Y=_DRY_AIR_Y))
    net.add_node(PressureBoundary("pb_outlet", Pt=40000.0, Tt=300.0, Y=_DRY_AIR_Y))

    op = _expected_operating_point()

    mc_primary = MomentumChamberNode("mc_primary", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_primary.initial_guess = {"mc_primary.P": 700000.0, "mc_primary.Pt": 700000.0}
    net.add_node(mc_primary)
    mc_secondary = MomentumChamberNode("mc_secondary", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_secondary.initial_guess = {"mc_secondary.P": 22850.0, "mc_secondary.Pt": 22850.0}
    net.add_node(mc_secondary)
    mc_outlet = MomentumChamberNode("mc_outlet", area=math.pi * (_CHANNEL_D / 2.0) ** 2)
    mc_outlet.initial_guess = {"mc_outlet.P": 40000.0, "mc_outlet.Pt": 40000.0}
    net.add_node(mc_outlet)

    ch_primary = ChannelElement(
        "ch_primary",
        "pb_primary",
        "mc_primary",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_primary.initial_guess = {"ch_primary.m_dot": op["mp"]}
    net.add_element(ch_primary)

    ch_secondary = ChannelElement(
        "ch_secondary",
        "pb_secondary",
        "mc_secondary",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_secondary.initial_guess = {"ch_secondary.m_dot": op["ms"]}
    net.add_element(ch_secondary)

    ch_outlet = ChannelElement(
        "ch_outlet",
        "mc_outlet",
        "pb_outlet",
        length=0.2,
        diameter=_CHANNEL_D,
        regime="incompressible",
    )
    ch_outlet.initial_guess = {"ch_outlet.m_dot": op["mout"]}
    net.add_element(ch_outlet)

    ej1 = EjectorElement(
        id="ej1",
        primary_node="mc_primary",
        secondary_node="mc_secondary",
        outlet_node="mc_outlet",
        throat_area=_THROAT_AREA,
        nozzle_exit_area=_NOZZLE_EXIT_AREA,
        mixing_area=_MIXING_AREA,
        regime="compressible",
        # K&D's anchors: these fixtures place the outlet inside its critical window;
        # huang1999 (the default) has its own tests (test_ejector_huang_anchors.py,
        # test_ejector_anchor_closure_element.py, test_ejector_subcritical_element.py).
        anchor_closure="kracik_dvorak",
    )
    ej1.initial_guess = {"ej1.P_jct": op["pc"]}
    net.add_element(ej1)
    return net


def test_full_network_solves_and_conserves_mass():
    net = _build_network()
    solver = NetworkSolver(net)
    sol = solver.solve(timeout=60.0)
    assert sol["__success__"], sol.get("__message__")

    mp = sol["ch_primary.m_dot"]
    ms = sol["ch_secondary.m_dot"]
    mout = sol["ch_outlet.m_dot"]
    assert mp > 0.0
    assert ms > 0.0
    assert mout == pytest.approx(mp + ms, rel=1e-4)

    # Entrainment ratio should be a plausible value, not a degenerate root.
    omega = ms / mp
    assert 0.0 < omega < 3.0

    # Critical mode: the actual outlet back pressure (imposed by pb_outlet,
    # 40 kPa) must not exceed the ejector's own P_c* (ej1.P_jct).
    assert sol["mc_outlet.Pt"] <= sol["ej1.P_jct"] * 1.01


@pytest.mark.parametrize(
    ("Pt", "gap"),
    [
        (2.0e5, 1e-6),  # above kJetDriveFloor (drive = 0.2 Pa): raw closure
        (1.0e7, 1e-9),  # above kJetDriveFloor (drive = 0.01 Pa): raw closure at 1-r = 1e-9
        (1.0e10, 1e-12),  # above kJetDriveFloor (drive = 0.01 Pa): raw closure at 1-r = 1e-12
        (2.0e5, 1e-9),  # below kJetDriveFloor (drive = 2e-4 Pa): #1269 sub-floor cubic
        (2.0e5, 1e-12),  # below kJetDriveFloor (drive = 2e-7 Pa): #1269 sub-floor cubic
    ],
)
def test_cd_nozzle_flow_and_slope_have_no_cancellation_near_zero_drive(
    Pt: float, gap: float
) -> None:
    """The C-D nozzle's subsonic exit Mach, M^2 = 2/(g-1) (r^(-(g-1)/g) - 1),
    cancelled as r = P_py / Pt -> 1: a relative error ~ eps / (1 - r) in the
    flow and its P_py slope (1.5e-7 at 1 - r = 1e-9, 1e-4 at 1e-12). Written
    with expm1/log1p it stays at rounding above kJetDriveFloor (1e-3 Pa), and
    below kJetDriveFloor takes the #1259/#1269 odd-cubic sub-floor law; both
    are checked against a 40-digit evaluation of the same closure."""
    mpmath = pytest.importorskip("mpmath")
    from combaero import _solver_tools

    mpmath.mp.dps = 40
    g, R, eta, Tt, At, Ae = 1.4, 287.0, 0.95, 300.0, 3.14e-5, 1.0e-4
    floor = mpmath.mpf("1e-3")
    p_py = Pt * (1.0 - gap)

    def flux_raw(p):
        gm, r = mpmath.mpf(g), mpmath.mpf(p) / mpmath.mpf(Pt)
        m = mpmath.sqrt(2 / (gm - 1) * (r ** (-(gm - 1) / gm) - 1))
        f = m * (1 + (gm - 1) / 2 * m * m) ** (-(gm + 1) / (2 * (gm - 1)))
        ex = Ae * mpmath.sqrt(gm / R) * mpmath.sqrt(eta) * mpmath.mpf(Pt) / mpmath.sqrt(Tt) * f
        cc = mpmath.sqrt(gm / R * (2 / (gm + 1)) ** ((gm + 1) / (gm - 1)))
        cap = cc * At * mpmath.sqrt(eta) * mpmath.mpf(Pt) / mpmath.sqrt(Tt)
        e = mpmath.mpf("1e-3") * cap
        return (ex + cap - mpmath.sqrt((ex - cap) ** 2 + e * e)) / 2

    def flux(p):
        drive = mpmath.mpf(Pt) - p
        if drive >= floor:
            return flux_raw(p)
        t = drive / floor
        h = mpmath.mpf("1.25") * t - mpmath.mpf("0.25") * t * t * t
        return flux_raw(mpmath.mpf(Pt) - floor) * h

    j = _solver_tools.ejector_cd_nozzle_mass_flow_and_jacobian(
        Pt, Tt, p_py, At, Ae, g, R, eta, 1e-3
    )
    ref = flux(mpmath.mpf(p_py))
    dref = mpmath.diff(flux, mpmath.mpf(p_py))
    assert j.mdot == pytest.approx(float(ref), rel=1e-8)
    assert j.dmdot_dp_py == pytest.approx(float(dref), rel=1e-8)


def test_python_cd_nozzle_mirror_matches_cpp_across_sub_floor_band() -> None:
    """Pin _ejector_huang1999.cd_nozzle_mass_flow to the C++ kernel across +-floor.

    Catches silent drift in the Python mirror's jet_floor (1e-3 Pa) or odd cubic
    (1.25 t - 0.25 t^3) relative to C++ (#1269).
    """
    from combaero.network._ejector_huang1999 import cd_nozzle_mass_flow

    Pt, Tt, At, Ae, g, R, eta, eps = 1.0e5, 300.0, 1.0e-4, 1.0e-4, 1.4, 287.058, 0.95, 1e-3
    floor = 1.0e-3
    for frac in (-2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0):
        p_py = Pt - frac * floor
        py_val = cd_nozzle_mass_flow(Pt, Tt, p_py, At, Ae, g, R, eta, eps)
        cpp_val = cb._core.ejector_cd_nozzle_mass_flow_and_jacobian(
            Pt, Tt, p_py, At, Ae, g, R, eta, eps
        ).mdot
        assert py_val == pytest.approx(cpp_val, rel=1e-11, abs=1e-15), f"frac={frac}"
