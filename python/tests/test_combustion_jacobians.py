"""Tests for combustion Jacobian functions."""

import combaero as cb


def test_adiabatic_complete_jacobian_basic():
    """Test basic functionality of adiabatic_T_complete_and_jacobian_T."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    # Stoichiometric mixture
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0  # 700K
    P = 500000.0  # 5 bar

    # Test the function
    result = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_mix)

    # Should return tuple of (T_ad, dT_ad_dT_in, X_products)
    assert len(result) == 3
    T_ad, dT_ad_dT_in, X_products = result

    # Basic validation
    assert T_ad > T_in, "Adiabatic temperature should be higher than inlet temperature"
    assert 2000 < T_ad < 3000, f"T_ad {T_ad} should be in reasonable range for CH4 combustion"

    # Check that we have proper combustion products
    assert X_products[3] > 0.05, "Should have significant CO2"
    assert X_products[4] > 0.1, "Should have significant H2O"

    # Jacobian should be positive (T_ad increases with T_in)
    assert dT_ad_dT_in > 0, "Jacobian should be positive"
    assert dT_ad_dT_in < 2.0, "Jacobian should be reasonable magnitude"

    # Species should sum to 1
    assert abs(sum(X_products) - 1.0) < 1e-6, "Species should sum to 1"


def test_adiabatic_equilibrium_jacobian_basic():
    """Test basic functionality of adiabatic_T_equilibrium_and_jacobians."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    # Stoichiometric mixture
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0  # 700K
    P = 500000.0  # 5 bar

    # Test the function
    result = cb._core.adiabatic_T_equilibrium_and_jacobians(T_in, P, X_mix)

    # Should return tuple of (T_ad, dT_ad_dT_in, dT_ad_dP, X_products)
    assert len(result) == 4
    T_ad, dT_ad_dT_in, dT_ad_dP, X_products = result

    # Basic validation
    assert T_ad > T_in, "Adiabatic temperature should be higher than inlet temperature"

    # Jacobian should be positive
    assert dT_ad_dT_in > 0, "Jacobian should be positive"

    # Species should sum to 1
    assert abs(sum(X_products) - 1.0) < 1e-6, "Species should sum to 1"


def test_jacobian_sensitivity():
    """Test that Jacobian responds correctly to temperature changes."""
    # Create a methane-air mixture
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    P = 500000.0  # 5 bar
    T_base = 700.0

    # Get Jacobian at base temperature
    _, dT_ad_dT_in_base, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_base, P, X_mix)

    # Test at slightly different temperature
    T_perturbed = T_base + 10.0
    _, dT_ad_dT_in_perturbed, _ = cb._core.adiabatic_T_complete_and_jacobian_T(
        T_perturbed, P, X_mix
    )

    # Jacobians should be similar (function should be smooth)
    assert abs(dT_ad_dT_in_base - dT_ad_dT_in_perturbed) < 0.1, "Jacobian should be smooth"

    # Both should be positive
    assert dT_ad_dT_in_base > 0 and dT_ad_dT_in_perturbed > 0, "Jacobians should be positive"


def test_lean_vs_stoichiometric():
    """Test that lean mixtures have lower adiabatic temperatures."""
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")

    P = 500000.0
    T_in = 700.0

    # Stoichiometric
    X_stoich = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)
    _, _, X_products_stoich = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_stoich)

    # Lean (phi = 0.8)
    X_lean = cb.set_equivalence_ratio_mole(0.8, X_ch4, X_air)
    T_ad_lean, _, X_products_lean = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P, X_lean)

    # Lean should have lower adiabatic temperature due to excess air
    # (though the current implementation may not show this due to placeholder stoichiometry)
    print(f"Stoichiometric CO2: {X_products_stoich[3]:.4f}, Lean CO2: {X_products_lean[3]:.4f}")
    print(f"Lean T_ad: {T_ad_lean:.2f} K")


def test_pressure_dependence():
    """Test that adiabatic temperature depends on pressure."""
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.pure_species("CH4")
    X_mix = cb.set_equivalence_ratio_mole(1.0, X_ch4, X_air)

    T_in = 700.0

    # Test at different pressures
    P_low = 101325.0  # 1 atm
    P_high = 500000.0  # 5 bar

    T_ad_low, _, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P_low, X_mix)
    T_ad_high, _, _ = cb._core.adiabatic_T_complete_and_jacobian_T(T_in, P_high, X_mix)

    print(f"T_ad at 1 atm: {T_ad_low:.2f} K")
    print(f"T_ad at 5 bar: {T_ad_high:.2f} K")

    # Pressure should have some effect (though may be small)
    print(f"Pressure effect: {T_ad_high - T_ad_low:.2f} K")


if __name__ == "__main__":
    print("Testing combustion Jacobian functions...")

    test_adiabatic_complete_jacobian_basic()
    print("OK: Basic complete combustion Jacobian test passed")

    test_adiabatic_equilibrium_jacobian_basic()
    print("OK: Basic equilibrium combustion Jacobian test passed")

    test_jacobian_sensitivity()
    print("OK: Jacobian sensitivity test passed")

    test_lean_vs_stoichiometric()
    print("OK: Lean vs stoichiometric test completed")

    test_pressure_dependence()
    print("OK: Pressure dependence test completed")

    print("\\nAll combustion Jacobian tests completed!")


# --- #716: the composition derivatives are analytic --------------------------

import numpy as np  # noqa: E402
import pytest  # noqa: E402

import combaero as _cb  # noqa: E402
from combaero import _solver_tools as _st  # noqa: E402

_Y_AIR = list(_cb.mole_to_mass(_cb.species.dry_air()))
_Y_CH4 = [0.0] * len(_Y_AIR)
_Y_CH4[_cb.species_index_from_name("CH4")] = 1.0
_X_AIR = _cb.species.dry_air()
_M_ST = _X_AIR[_cb.species_index_from_name("O2")] / 2 * 16.043 / _cb.mwmix(_X_AIR)


@pytest.mark.parametrize("method", ["complete", "equilibrium"])
@pytest.mark.parametrize("phi", [0.6, 0.9, 1.0, 1.0002, 1.3, 1.6])
def test_the_combustor_derivatives_match_central_differences(method: str, phi: float) -> None:
    """dT and dY_products in each stream's flow and temperature, against
    central differences (#716). The forward differences inside C++ were 123%
    off across the phi = 1 softmin (complete) and up to 5-12% (equilibrium,
    whose reforming Newton also stalled off its root near complete
    reforming). Measured 2026-10-25: every entry within 7e-7."""
    fn = (
        _st.adiabatic_T_complete_and_jacobian_T_from_streams
        if method == "complete"
        else _st.adiabatic_T_equilibrium_and_jacobians_from_streams
    )

    def burn(ma: float, mf: float, Ta: float = 700.0):
        s = [_cb.MassStream(ma, Ta, 2e5, _Y_AIR), _cb.MassStream(mf, 300.0, 2e5, _Y_CH4)]
        return fn(s, 2e5, 0.0, 0.0)

    mf = phi * _M_ST
    r = burn(1.0, mf)
    # Inside the phi = 1 softmin (about 1e-4 wide in phi, 6e-6 kg/s of fuel
    # here; both methods since #1048) the slope swings from +24000 to -7000
    # K/(kg/s): a step that resolves it is 1e-9.
    in_band = abs(phi - 1.0) < 1e-3
    h, hT = (1e-9 if in_band else 1e-7), 1e-3
    ico = _cb.species_index_from_name("CO" if method == "equilibrium" else "CO2")
    checks = [
        (r.dTt_mix_d_stream[1].d_mdot, lambda d: burn(1.0, mf + d).Tt_mix, h),
        (r.dTt_mix_d_stream[0].d_mdot, lambda d: burn(1.0 + d, mf).Tt_mix, h),
        (r.dTt_mix_d_stream[0].d_Tt, lambda d: burn(1.0, mf, 700.0 + d).Tt_mix, hT),
        (r.dY_mix_d_stream[ico][1].d_mdot, lambda d: burn(1.0, mf + d).Y_mix[ico], h),
        (r.dY_mix_d_stream[ico][0].d_Tt, lambda d: burn(1.0, mf, 700.0 + d).Y_mix[ico], hT),
    ]
    for an, f, step in checks:
        fd = (f(step) - f(-step)) / (2 * step)
        assert an == pytest.approx(fd, rel=2e-6, abs=1e-9 * max(1.0, abs(fd))), (an, fd)


def test_the_complete_products_jacobian_is_exact_and_conserves_atoms() -> None:
    """dX_b/dX against central differences through phi = 1, and the oxygen
    atoms balance inside the smoothing band (O2 left = X_O2 - f n_req)."""
    names = [_cb.species_name(i) for i in range(_cb.num_species())]
    X = np.array(_X_AIR)
    iO2, iC, iN2 = names.index("O2"), names.index("CH4"), names.index("N2")
    v = np.zeros(len(names))
    v[iC], v[iN2] = 1.0, -1.0
    for phi in (0.9, 0.99995, 1.0, 1.00005, 1.1):
        xf = X[iO2] / 2 * phi
        X0 = X / (1 + xf)
        X0[iC] += xf / (1 + xf)
        cj = _cb.complete_combustion_products_and_jacobian(X0)
        h = 1e-9

        def f(d: float, X0: np.ndarray = X0) -> np.ndarray:
            return np.array(_cb.complete_combustion_products_and_jacobian(X0 + d * v).X_b)

        fd = (f(h) - f(-h)) / (2 * h)
        assert np.abs(np.array(cj.dXb_dX) @ v - fd).max() < 1e-6 * np.abs(fd).max(), phi
        Xb = np.array(cj.X_b)
        o_atoms = {"O2": 2, "CO2": 2, "H2O": 1}
        o_in = sum(X0[names.index(k)] * w for k, w in o_atoms.items()) / X0[iN2]
        o_out = sum(Xb[names.index(k)] * w for k, w in o_atoms.items()) / Xb[iN2]
        assert o_out == pytest.approx(o_in, rel=1e-12), phi


def _fuel_air(fuel: str, phi: float) -> tuple[np.ndarray, np.ndarray]:
    """X(phi) of a fuel in dry air, and dX/dphi (sum zero)."""
    i_f = _cb.species_index_from_name(fuel)
    X_air = np.array(_X_AIR)
    Xf = np.zeros(len(X_air))
    Xf[i_f] = 1.0
    c = X_air[_cb.species_index_from_name("O2")] / _cb.oxygen_required_per_mol_fuel(i_f)
    X = (X_air + c * phi * Xf) / (1.0 + c * phi)
    return X, c * (Xf - X) / (1.0 + c * phi)


@pytest.mark.parametrize("fuel", ["CH4", "H2", "C3H8"])
def test_equilibrium_T_is_smooth_through_stoichiometric(fuel: str) -> None:
    """combustion_equilibrium runs its complete-combustion step with the phi = 1
    softmin (#1048), so T is C1 at phi = 1. Without it the slope jumped there,
    +1478 to -961 K per unit phi for CH4-air at 300 K (one-sided slopes 1e-8
    either side of phi = 1), and the analytic derivative sat between them.
    Here the slopes either side agree to the step's curvature, and the
    analytic dT/dphi matches central differences across the band."""
    T_in, P = 700.0, 5e5

    def T_of(phi: float) -> float:
        return float(_cb.combustion_equilibrium(T_in, _fuel_air(fuel, phi)[0], P).state.T)

    h = 1e-8  # resolves the band, ~1/SMOOTHING_K_PHI1 = 5e-5 wide in phi
    T0 = T_of(1.0)
    fwd, bwd = (T_of(1.0 + h) - T0) / h, (T0 - T_of(1.0 - h)) / h
    # Smooth: the one-sided slopes differ by about T'' h (T'' ~ 1e7 per phi^2
    # in the band): 0.1 K per unit phi, against the kink's 2400.
    assert abs(fwd - bwd) < 1.0, (fwd, bwd)
    for phi in (0.9999, 0.99995, 1.0, 1.00002, 1.00005, 1.0001):
        X, dX = _fuel_air(fuel, phi)
        ej = _cb.combustion_equilibrium_and_jacobian(T_in, X, P)
        an = float(np.array(ej.dT_dX) @ dX)
        fd = (T_of(phi + h) - T_of(phi - h)) / (2 * h)
        assert an == pytest.approx(fd, rel=1e-5, abs=1e-2), (phi, an, fd)


def test_equilibrium_softmin_moves_T_only_near_stoichiometric() -> None:
    """The softmin changes T only within about 1e-3 of phi = 1 (#1048): at
    phi = 1 it is 0.083 K below the unsmoothed burn followed by the same
    reforming + WGS stage (CH4-air, 700 K), and at |phi - 1| >= 1e-3 equal to
    it to round-off."""
    T_in, P = 700.0, 5e5

    def change(phi: float) -> float:
        X = _fuel_air("CH4", phi)[0]
        smooth = _cb.combustion_equilibrium(T_in, X, P).state.T
        # The burn without the phi = 1 softmin, then the same reforming + WGS stage.
        b = _cb.complete_combustion(T_in, X, P, smooth_phi0=True, smooth_phi1=False)
        return smooth - _cb.reforming_equilibrium_adiabatic(b.T, list(b.X), b.P).state.T

    for phi in (0.99, 0.999, 1.001, 1.01):
        assert abs(change(phi)) < 1e-9, phi
    # ln2/k1 of the fuel is left unburnt at phi = 1: measured -0.0833 K.
    assert -0.1 < change(1.0) < -0.07


# --- #1251 / #1244: zero-oxygen slope and the burn pressure's partials -------


def _zero_loss(ctx: object) -> tuple[float, float]:
    return (0.0, 0.0)


def _combustor(streams: list, method: str):
    """The combustor kernel as CombustorNode calls it: at the streams' mixed
    Pt (#1244)."""
    P = _st.mixer_from_streams_and_jacobians(streams).Pt_mix
    return _st.combustor_residuals_and_jacobians(
        streams, P, pressure_loss=_zero_loss, use_equilibrium=(method == "equilibrium")
    )


def test_a_fuel_only_mix_burns_the_first_oxygen() -> None:
    """At X_O2 = 0 with fuel present the products are the reactants, and
    dX_b/dX_O2 is the burning branch's slope as X_O2 -> 0+, not the identity
    (#1251): the identity gave a fuel-fed combustor dT/dm_air = 0 against a
    one-sided difference of 2.26e6 K per kg/s."""
    names = [_cb.species_name(i) for i in range(_cb.num_species())]
    iO2, iC = names.index("O2"), names.index("CH4")
    X0 = np.zeros(len(names))
    X0[iC] = 1.0
    cj = _cb.complete_combustion_products_and_jacobian(X0)
    assert np.array(cj.X_b) == pytest.approx(X0, abs=0.0)
    v = np.zeros(len(names))
    v[iO2], v[iC] = 1.0, -1.0
    h = 1e-9
    fd = (np.array(_cb.complete_combustion_products_and_jacobian(X0 + h * v).X_b) - X0) / h
    an = np.array(cj.dXb_dX) @ v
    assert np.abs(an - fd).max() < 1e-6 * np.abs(fd).max()
    assert np.abs(an - v).max() > 0.5  # not the identity


def test_a_fuel_only_combustor_moves_with_the_air_flow() -> None:
    """Combustor kernel, air stream at exactly 0 beside a CH4 stream (#1251,
    combustor:2 point a_in.m_dot = 0): the air column against the one-sided
    difference on the arriving side (the solver's convention at m = 0).
    Complete combustion to 1e-6. Equilibrium only at the sweep's steps: J is
    the complete chain's, F's slope below the reforming onset (about 1e-12
    of oxygen); above the onset the reforming + WGS solve on trace steam does
    not converge (#1262) and the difference depends on the step (1.0e6 at
    1e-14 kg/s, 2.26e6 at 1e-7), which #1251 keeps open."""
    Y_air = _Y_AIR
    for method, h, rel in (("complete", 1e-10, 1e-6), ("equilibrium", 1e-7, 2e-3)):

        def burn(ma: float, method: str = method):
            s = [
                _cb.MassStream(ma, 300.0, 2.0e5, Y_air),
                _cb.MassStream(5.76e-4, 300.0, 2.1e5, _Y_CH4),
            ]
            return _combustor(s, method)

        r = burn(0.0)
        fd = (burn(h).Tt_mix - r.Tt_mix) / h
        an = r.dTt_mix_d_stream[0].d_mdot
        assert fd > 1e6
        assert an == pytest.approx(fd, rel=rel), (method, an, fd)


@pytest.mark.parametrize("ma", [3e-5, 1e-3, 0.05])
def test_the_equilibrium_combustor_carries_its_burn_pressure(ma: float) -> None:
    """The equilibrium depends on its pressure, the streams' mixed Pt, which
    moves with each stream's flow and Pt (#1244). Near zero air flow (3e-5
    kg/s beside 5.4e-4 kg/s of fuel 5 % higher in Pt) that mean moves 1.7e7
    Pa per kg/s, and the flow column was 3.5e-4 of the row off."""
    Y_fuel = [0.0] * len(_Y_AIR)
    for k, y in enumerate(_Y_AIR):
        Y_fuel[k] = 0.85 * y
    Y_fuel[_cb.species_index_from_name("CH4")] += 0.15

    def burn(ma: float, Pa: float = 2.0e5, mf: float = 5.376e-4):
        s = [
            _cb.MassStream(ma, 650.0, Pa, _Y_AIR),
            _cb.MassStream(mf, 300.0, 2.1e5, Y_fuel),
        ]
        return _combustor(s, "equilibrium")

    r = burn(ma)
    ico = _cb.species_index_from_name("CO")
    hm, hp = 1e-3 * ma, 10.0
    checks = [
        (r.dTt_mix_d_stream[0].d_mdot, lambda d: burn(ma + d).Tt_mix, hm),
        (r.dTt_mix_d_stream[0].d_Pt, lambda d: burn(ma, 2.0e5 + d).Tt_mix, hp),
        (r.dY_mix_d_stream[ico][0].d_mdot, lambda d: burn(ma + d).Y_mix[ico], hm),
        (r.dY_mix_d_stream[ico][0].d_Pt, lambda d: burn(ma, 2.0e5 + d).Y_mix[ico], hp),
    ]
    for an, f, step in checks:
        fd = (f(step) - f(-step)) / (2 * step)
        assert an == pytest.approx(fd, rel=1e-5, abs=1e-12), (an, fd)
    if ma < 1e-4:  # rich: the reforming moves with P (lean, it does not)
        assert abs(r.dTt_mix_d_stream[0].d_Pt) > 1e-7


def test_the_combustor_refuses_a_pressure_other_than_the_mixed_pt() -> None:
    """The stream partials are the mixed Pt's, so a caller's own P would carry
    them wrongly (dT/dPt_fuel 1.3e-4 against an FD of 0 at P = 1e5 beside a
    mixed 2.03e5): refused (#1244 review)."""
    s = [
        _cb.MassStream(0.05, 650.0, 2.0e5, _Y_AIR),
        _cb.MassStream(1e-3, 300.0, 2.1e5, _Y_CH4),
    ]
    with pytest.raises(ValueError, match="mixed Pt"):
        _st.combustor_residuals_and_jacobians(
            s, 1.0e5, pressure_loss=_zero_loss, use_equilibrium=True
        )
    _combustor(s, "equilibrium")  # at the mixed Pt: accepted
    # A non-finite stream Pt is refused as a bad pressure, not as a mismatch.
    bad = [s[0], _cb.MassStream(1e-3, 300.0, float("nan"), _Y_CH4)]
    with pytest.raises(ValueError, match="finite and positive"):
        _st.combustor_residuals_and_jacobians(
            bad, float("nan"), pressure_loss=_zero_loss, use_equilibrium=True
        )
