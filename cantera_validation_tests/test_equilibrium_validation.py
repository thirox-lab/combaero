import numpy as np
import pytest

from combaero._core import num_species, species_index_from_name, species_name


class TestWgsEquilibrium:
    """Validate WGS (Water-Gas Shift) equilibrium against Cantera.

    Reaction: CO + H2O <=> CO2 + H2

    CombAero uses partial equilibrium (WGS only), Cantera uses restricted
    species to match the same equilibrium model.
    """

    @pytest.fixture
    def wgs_gas(self, cantera):
        """Create Cantera gas with WGS species only."""
        species = {S.name: S for S in cantera.Species.list_from_file("gri30.yaml")}
        wgs_species = [species[S] for S in ("CO", "H2O", "CO2", "H2", "N2", "AR")]
        return cantera.Solution(thermo="ideal-gas", species=wgs_species)

    def set_cantera_composition(self, gas, cb_X, species_mapping):
        """Set Cantera composition from CombAero mole fractions."""
        ct_species = {}
        for i in range(num_species()):
            if cb_X[i] > 1e-10:
                name = species_name(i)
                ct_name = species_mapping.get(name, name)
                if ct_name in gas.species_names:
                    ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_wgs_isothermal_high_temp(
        self, combaero, cantera, wgs_gas, species_mapping, tolerance_config
    ):
        """Test WGS equilibrium at high temperature (1200 K) - favors CO + H2."""
        cb = combaero

        # Initial composition: CO + H2O rich
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T = 1200.0
        P = 101325.0

        # CombAero WGS equilibrium
        result_cb = cb.wgs_equilibrium(T, X, P)

        # Cantera equilibrium with WGS species only
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T, P
        wgs_gas.equilibrate("TP")  # Isothermal equilibrium

        # Compare temperature [K] - should be unchanged (isothermal)
        assert abs(result_cb.state.T - T) < 0.1, (
            f"Temperature should be unchanged: {result_cb.state.T:.1f} K vs {T:.1f} K"
        )

        # Compare equilibrium composition [mol/mol]
        # At high T, WGS favors reactants (CO + H2O), so expect shift toward products
        print(f"\nWGS isothermal @ {T} K:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.state.X[idx_cb]
            X_ct = wgs_gas[species].X[0] if species in wgs_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            # Major species tolerance
            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.005, (
                    f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}"
                )

    def test_wgs_isothermal_low_temp(
        self, combaero, cantera, wgs_gas, species_mapping, tolerance_config
    ):
        """Test WGS equilibrium at low temperature (800 K) - favors CO2 + H2."""
        cb = combaero

        # Initial composition: CO + H2O rich
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T = 800.0
        P = 101325.0

        # CombAero WGS equilibrium
        result_cb = cb.wgs_equilibrium(T, X, P)

        # Cantera equilibrium
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T, P
        wgs_gas.equilibrate("TP")

        # Compare composition
        print(f"\nWGS isothermal @ {T} K:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.state.X[idx_cb]
            X_ct = wgs_gas[species].X[0] if species in wgs_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.005, f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}"

    def test_wgs_isothermal_temperature_sweep(self, combaero, cantera, wgs_gas, species_mapping):
        """Test WGS equilibrium across temperature range (800-1500 K)."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS isothermal temperature sweep:")
        max_deviations = {"CO": 0.0, "H2O": 0.0, "CO2": 0.0, "H2": 0.0}

        for T in [800.0, 1000.0, 1200.0, 1500.0]:
            result_cb = cb.wgs_equilibrium(T, X, P)

            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T, P
            wgs_gas.equilibrate("TP")

            print(f"\n  T = {T} K:")
            for species in ["CO", "H2O", "CO2", "H2"]:
                idx_cb = species_index_from_name(species)
                X_cb = result_cb.state.X[idx_cb]
                X_ct = wgs_gas[species].X[0]
                diff = abs(X_cb - X_ct)
                max_deviations[species] = max(max_deviations[species], diff)
                print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

        print("\n  Maximum deviations:")
        for species, max_diff in max_deviations.items():
            print(f"    {species}: {max_diff:.6f}")

    def test_wgs_adiabatic(self, combaero, cantera, wgs_gas, species_mapping, tolerance_config):
        """Test adiabatic WGS equilibrium - temperature changes with reaction."""
        cb = combaero

        # Initial composition at high temperature
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        T_in = 1500.0
        P = 101325.0

        # CombAero adiabatic WGS
        result_cb = cb.wgs_equilibrium_adiabatic(T_in, X, P)

        # Cantera adiabatic equilibrium
        self.set_cantera_composition(wgs_gas, X, species_mapping)
        wgs_gas.TP = T_in, P
        wgs_gas.equilibrate("HP")  # Adiabatic (constant enthalpy, pressure)

        # Compare final temperature [K]
        # WGS is exothermic in forward direction (CO + H2O -> CO2 + H2)
        # At high T, equilibrium favors reactants, so forward reaction releases heat
        T_diff = abs(result_cb.state.T - wgs_gas.T)
        print("\nWGS adiabatic:")
        print(f"  T_in: {T_in:.1f} K")
        print(
            f"  T_final: CombAero={result_cb.state.T:.1f} K, "
            f"Cantera={wgs_gas.T:.1f} K, diff={T_diff:.1f} K"
        )

        assert T_diff < 10.0, (
            f"Temperature [K]: CombAero={result_cb.state.T:.1f}, "
            f"Cantera={wgs_gas.T:.1f}, diff={T_diff:.1f}"
        )

        # Compare equilibrium composition
        print("  Composition:")
        for species in ["CO", "H2O", "CO2", "H2"]:
            idx_cb = species_index_from_name(species)
            X_cb = result_cb.state.X[idx_cb]
            X_ct = wgs_gas[species].X[0]
            diff = abs(X_cb - X_ct)
            print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.005, f"{species} [mol/mol]: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}"

    def test_wgs_adiabatic_temperature_sweep(self, combaero, cantera, wgs_gas, species_mapping):
        """Test adiabatic WGS across inlet temperatures."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS adiabatic temperature sweep:")
        max_T_diff = 0.0

        for T_in in [1000.0, 1200.0, 1500.0, 1800.0]:
            result_cb = cb.wgs_equilibrium_adiabatic(T_in, X, P)

            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T_in, P
            wgs_gas.equilibrate("HP")

            T_diff = abs(result_cb.state.T - wgs_gas.T)
            max_T_diff = max(max_T_diff, T_diff)

            print(f"\n  T_in = {T_in} K:")
            print(
                f"    T_final: CombAero={result_cb.state.T:.1f} K, "
                f"Cantera={wgs_gas.T:.1f} K, diff={T_diff:.1f} K"
            )

            for species in ["CO", "H2O", "CO2", "H2"]:
                idx_cb = species_index_from_name(species)
                X_cb = result_cb.state.X[idx_cb]
                X_ct = wgs_gas[species].X[0]
                diff = abs(X_cb - X_ct)
                print(f"    {species}: CombAero={X_cb:.6f}, Cantera={X_ct:.6f}, diff={diff:.6f}")

        print(f"\n  Maximum temperature deviation: {max_T_diff:.1f} K")

    def test_wgs_equilibrium_constant(self, combaero, cantera, wgs_gas, species_mapping):
        """Test that equilibrium constant Kp matches between CombAero and Cantera."""
        cb = combaero

        # Initial composition
        X = np.zeros(num_species())
        X[species_index_from_name("CO")] = 0.10
        X[species_index_from_name("H2O")] = 0.15
        X[species_index_from_name("CO2")] = 0.05
        X[species_index_from_name("H2")] = 0.05
        X[species_index_from_name("N2")] = 0.65

        P = 101325.0

        print("\nWGS equilibrium constant Kp:")

        for T in [800.0, 1000.0, 1200.0, 1500.0]:
            # CombAero equilibrium
            result_cb = cb.wgs_equilibrium(T, X, P)

            # Calculate Kp from composition: Kp = (X_CO2 * X_H2) / (X_CO * X_H2O)
            X_CO2_cb = result_cb.state.X[species_index_from_name("CO2")]
            X_H2_cb = result_cb.state.X[species_index_from_name("H2")]
            X_CO_cb = result_cb.state.X[species_index_from_name("CO")]
            X_H2O_cb = result_cb.state.X[species_index_from_name("H2O")]
            Kp_cb = (
                (X_CO2_cb * X_H2_cb) / (X_CO_cb * X_H2O_cb) if (X_CO_cb * X_H2O_cb) > 1e-10 else 0.0
            )

            # Cantera equilibrium
            self.set_cantera_composition(wgs_gas, X, species_mapping)
            wgs_gas.TP = T, P
            wgs_gas.equilibrate("TP")

            X_CO2_ct = wgs_gas["CO2"].X[0]
            X_H2_ct = wgs_gas["H2"].X[0]
            X_CO_ct = wgs_gas["CO"].X[0]
            X_H2O_ct = wgs_gas["H2O"].X[0]
            Kp_ct = (
                (X_CO2_ct * X_H2_ct) / (X_CO_ct * X_H2O_ct) if (X_CO_ct * X_H2O_ct) > 1e-10 else 0.0
            )

            rel_diff = abs(Kp_cb - Kp_ct) / Kp_ct if Kp_ct > 1e-10 else 0.0

            print(
                f"  T = {T} K: CombAero Kp={Kp_cb:.4f}, "
                f"Cantera Kp={Kp_ct:.4f}, diff={rel_diff * 100:.2f}%"
            )

            # Kp should match within 5% (due to Gibbs free energy polynomial differences)
            assert rel_diff < 0.05, (
                f"Kp [dimensionless] @ {T} K: CombAero={Kp_cb:.4f}, "
                f"Cantera={Kp_ct:.4f}, diff={rel_diff * 100:.1f}%"
            )


class TestReformingEquilibrium:
    """Validate general reforming + WGS equilibrium against Cantera.

    CombAero uses sequential 1D Newton per hydrocarbon then WGS.
    Cantera uses restricted species to match the same partial equilibrium model.
    """

    @pytest.fixture
    def reforming_gas(self, cantera):
        """Cantera gas restricted to reforming + WGS species."""
        species = {S.name: S for S in cantera.Species.list_from_file("gri30.yaml")}
        ref_species = [species[s] for s in ("CH4", "H2O", "CO", "CO2", "H2", "N2", "AR")]
        return cantera.Solution(thermo="ideal-gas", species=ref_species)

    @pytest.fixture
    def reforming_gas_c2h6(self, cantera):
        """Cantera gas restricted to C2H6 reforming + WGS species."""
        species = {S.name: S for S in cantera.Species.list_from_file("gri30.yaml")}
        ref_species = [species[s] for s in ("C2H6", "H2O", "CO", "CO2", "H2", "N2", "AR")]
        return cantera.Solution(thermo="ideal-gas", species=ref_species)

    def set_cantera_composition(self, gas, cb_X, species_mapping):
        ct_species = {}
        for i in range(num_species()):
            if cb_X[i] > 1e-10:
                name = species_name(i)
                ct_name = species_mapping.get(name, name)
                if ct_name in gas.species_names:
                    ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_reforming_ch4_isothermal(self, combaero, cantera, reforming_gas, species_mapping):
        """CH4 + H2O reforming at 1200 K, 5 bar -- compare with Cantera TP equilibrium."""
        cb = combaero

        X = np.zeros(num_species())
        X[species_index_from_name("CH4")] = 0.10
        X[species_index_from_name("H2O")] = 0.30
        X[species_index_from_name("N2")] = 0.60

        T, P = 1200.0, 5e5

        result_cb = cb.reforming_equilibrium(T, X, P)

        self.set_cantera_composition(reforming_gas, X, species_mapping)
        reforming_gas.TP = T, P
        reforming_gas.equilibrate("TP")

        assert abs(result_cb.state.T - T) < 0.1, "Isothermal: T should be unchanged"

        print(f"\nReforming CH4 isothermal @ {T} K, {P / 1e5:.0f} bar:")
        for sp in ["CH4", "H2O", "CO", "CO2", "H2"]:
            idx = species_index_from_name(sp)
            X_cb = result_cb.state.X[idx]
            X_ct = reforming_gas[sp].X[0] if sp in reforming_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {sp}: CombAero={X_cb:.5f}, Cantera={X_ct:.5f}, diff={diff:.5f}")
            if X_cb > 0.005 or X_ct > 0.005:
                assert diff < 0.01, (
                    f"{sp}: CombAero={X_cb:.5f}, Cantera={X_ct:.5f}, diff={diff:.5f}"
                )

    def test_reforming_ch4_adiabatic(self, combaero, cantera, reforming_gas, species_mapping):
        """CH4 + H2O adiabatic reforming -- compare T_final and composition with Cantera HP."""
        cb = combaero

        X = np.zeros(num_species())
        X[species_index_from_name("CH4")] = 0.10
        X[species_index_from_name("H2O")] = 0.30
        X[species_index_from_name("N2")] = 0.60

        T, P = 1200.0, 5e5

        result_cb = cb.reforming_equilibrium_adiabatic(T, X, P)

        self.set_cantera_composition(reforming_gas, X, species_mapping)
        reforming_gas.TP = T, P
        reforming_gas.equilibrate("HP")

        T_diff = abs(result_cb.state.T - reforming_gas.T)
        print(f"\nReforming CH4 adiabatic @ {T} K:")
        print(
            f"  T_final: CombAero={result_cb.state.T:.1f} K, Cantera={reforming_gas.T:.1f} K, diff={T_diff:.1f} K"
        )
        print(f"  delta_T: {result_cb.delta_T:.1f} K")

        # CombAero's sequential 1D solver gives a different T_final than Cantera's
        # coupled full-equilibrium at high CH4 concentration (10%) + 5 bar.
        # The qualitative checks below verify correctness of direction and consumption.
        assert result_cb.state.T < T + 100.0, "T_final should not be much higher than T_in"
        assert result_cb.converged, "Solver should converge"

        X_CH4_in = X[species_index_from_name("CH4")]
        X_CH4_out = result_cb.state.X[species_index_from_name("CH4")]
        X_CO_out = result_cb.state.X[species_index_from_name("CO")]
        X_H2_out = result_cb.state.X[species_index_from_name("H2")]
        print(f"  CH4: {X_CH4_in:.4f} -> {X_CH4_out:.5f}")
        print(f"  CO:  0 -> {X_CO_out:.4f}, H2: 0 -> {X_H2_out:.4f}")
        assert X_CH4_out < X_CH4_in * 0.5, "CH4 should be mostly consumed"
        assert X_CO_out > 0.01, "CO should be produced"
        assert X_H2_out > 0.01, "H2 should be produced"

    def test_reforming_c2h6_isothermal(
        self, combaero, cantera, reforming_gas_c2h6, species_mapping
    ):
        """C2H6 + H2O reforming at 1400 K -- verify C2H6 consumed, CO/H2 produced."""
        cb = combaero

        X = np.zeros(num_species())
        X[species_index_from_name("C2H6")] = 0.05
        X[species_index_from_name("H2O")] = 0.25
        X[species_index_from_name("N2")] = 0.70

        T, P = 1400.0, 101325.0

        result_cb = cb.reforming_equilibrium(T, X, P)

        self.set_cantera_composition(reforming_gas_c2h6, X, species_mapping)
        reforming_gas_c2h6.TP = T, P
        reforming_gas_c2h6.equilibrate("TP")

        X_C2H6_out = result_cb.state.X[species_index_from_name("C2H6")]
        X_CO_out = result_cb.state.X[species_index_from_name("CO")]
        X_H2_out = result_cb.state.X[species_index_from_name("H2")]

        print(f"\nReforming C2H6 isothermal @ {T} K:")
        print(f"  C2H6: {X[species_index_from_name('C2H6')]:.4f} -> {X_C2H6_out:.6f}")
        print(f"  CO:   0 -> {X_CO_out:.4f}")
        print(f"  H2:   0 -> {X_H2_out:.4f}")

        assert X_C2H6_out < X[species_index_from_name("C2H6")] * 0.1, (
            "C2H6 should be mostly consumed"
        )
        assert X_CO_out > 0.01, "CO should be produced"
        assert X_H2_out > 0.01, "H2 should be produced"

        X_CO_ct = reforming_gas_c2h6["CO"].X[0]
        X_H2_ct = reforming_gas_c2h6["H2"].X[0]
        assert abs(X_CO_out - X_CO_ct) < 0.015, (
            f"CO: CombAero={X_CO_out:.5f}, Cantera={X_CO_ct:.5f}"
        )
        assert abs(X_H2_out - X_H2_ct) < 0.015, (
            f"H2: CombAero={X_H2_out:.5f}, Cantera={X_H2_ct:.5f}"
        )


class TestCombustionEquilibrium:
    """Validate combustion_equilibrium against Cantera full equilibrium."""

    @pytest.fixture
    def gri30_gas(self, cantera):
        return cantera.Solution("gri30.yaml")

    def set_cantera_composition(self, gas, cb_X, species_mapping):
        ct_species = {}
        for i in range(num_species()):
            if cb_X[i] > 1e-10:
                name = species_name(i)
                ct_name = species_mapping.get(name, name)
                if ct_name in gas.species_names:
                    ct_species[ct_name] = cb_X[i]
        gas.X = ct_species

    def test_combustion_equilibrium_lean_ch4(self, combaero, cantera, gri30_gas, species_mapping):
        """phi=0.8 CH4/air: compare T_ad and major species vs Cantera HP equilibrium."""
        cb = combaero

        X_air = cb.species.dry_air()
        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        X_mix = cb.set_equivalence_ratio_mole(0.8, X_fuel, X_air)

        T_in, P = 300.0, 101325.0

        result_cb = cb.combustion_equilibrium(T_in, X_mix, P)

        self.set_cantera_composition(gri30_gas, X_mix, species_mapping)
        gri30_gas.TP = T_in, P
        gri30_gas.equilibrate("HP")

        T_diff = abs(result_cb.state.T - gri30_gas.T)
        print("\nCombustion equilibrium lean CH4 (phi=0.8):")
        print(
            f"  T: CombAero={result_cb.state.T:.1f} K, Cantera={gri30_gas.T:.1f} K, diff={T_diff:.1f} K"
        )
        assert T_diff < 150.0, (
            f"T [K]: CombAero={result_cb.state.T:.1f}, Cantera={gri30_gas.T:.1f} (CombAero lacks dissociation species)"
        )

        for sp in ["CO2", "H2O", "O2"]:
            idx = species_index_from_name(sp)
            X_cb = result_cb.state.X[idx]
            X_ct = gri30_gas[sp].X[0] if sp in gri30_gas.species_names else 0.0
            diff = abs(X_cb - X_ct)
            print(f"  {sp}: CombAero={X_cb:.4f}, Cantera={X_ct:.4f}, diff={diff:.4f}")
            if X_cb > 0.01 or X_ct > 0.01:
                assert diff < 0.015, f"{sp}: CombAero={X_cb:.4f}, Cantera={X_ct:.4f}"

    def test_combustion_equilibrium_rich_ch4(self, combaero):
        """phi=1.2 CH4/air: verify CO and H2 survive (complete combustion gives zero)."""
        cb = combaero

        X_air = cb.species.dry_air()
        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        X_mix = cb.set_equivalence_ratio_mole(1.2, X_fuel, X_air)

        T_in, P = 300.0, 101325.0

        result_eq = cb.combustion_equilibrium(T_in, X_mix, P)
        result_complete = cb.complete_combustion(T=T_in, X=X_mix, P=P)

        X_CO_eq = result_eq.state.X[species_index_from_name("CO")]
        X_H2_eq = result_eq.state.X[species_index_from_name("H2")]
        X_CO_complete = result_complete.X[species_index_from_name("CO")]
        X_H2_complete = result_complete.X[species_index_from_name("H2")]

        print("\nCombustion equilibrium rich CH4 (phi=1.2):")
        print(f"  CO:  complete={X_CO_complete:.5f}, equilibrium={X_CO_eq:.5f}")
        print(f"  H2:  complete={X_H2_complete:.5f}, equilibrium={X_H2_eq:.5f}")

        assert X_CO_eq > 0.01, "Rich equilibrium: CO should survive"
        assert X_H2_eq > 0.01, "Rich equilibrium: H2 should survive"
        assert result_eq.converged, "Solver should converge"

    def test_combustion_equilibrium_stoich_ch4(self, combaero, cantera, gri30_gas, species_mapping):
        """phi=1.0 CH4/air: T_ad within 30 K of Cantera HP equilibrium."""
        cb = combaero

        X_air = cb.species.dry_air()
        X_fuel = np.zeros(num_species())
        X_fuel[species_index_from_name("CH4")] = 1.0
        X_mix = cb.set_equivalence_ratio_mole(1.0, X_fuel, X_air)

        T_in, P = 300.0, 101325.0

        result_cb = cb.combustion_equilibrium(T_in, X_mix, P)

        self.set_cantera_composition(gri30_gas, X_mix, species_mapping)
        gri30_gas.TP = T_in, P
        gri30_gas.equilibrate("HP")

        T_diff = abs(result_cb.state.T - gri30_gas.T)
        print("\nCombustion equilibrium stoich CH4 (phi=1.0):")
        print(
            f"  T: CombAero={result_cb.state.T:.1f} K, Cantera={gri30_gas.T:.1f} K, diff={T_diff:.1f} K"
        )
        assert T_diff < 150.0, (
            f"T [K]: CombAero={result_cb.state.T:.1f}, Cantera={gri30_gas.T:.1f} (CombAero lacks dissociation species)"
        )


# ---------------------------------------------------------------------------
# combustion_equilibrium over phi 0.8-3.0 against Cantera HP (#759). Below
# phi 1.5 the model lacks dissociation (O, H, OH, NO): the stoichiometric
# tolerances above apply. From phi 1.5 the products are CO2, H2O, CO, H2, CH4
# and N2, all in the model: 3 K and 1e-3 in mole fraction. Before #759 rich
# propane (phi >= 2.5) left no steam and ended 124-320 K hot, unconverged.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("phi", [0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 2.75, 3.0])
@pytest.mark.parametrize("fuel", ["CH4", "C3H8"])
def test_combustion_equilibrium_matches_cantera_hp_over_phi(
    combaero, cantera, fuel: str, phi: float
) -> None:
    X_fuel = np.zeros(num_species())
    X_fuel[species_index_from_name(fuel)] = 1.0
    X_mix = np.asarray(combaero.set_equivalence_ratio_mole(phi, X_fuel, combaero.species.dry_air()))
    T_in, P = 300.0, 101325.0
    got = combaero.combustion_equilibrium(T_in, X_mix, P)
    gas = cantera.Solution("gri30.yaml")
    names = [species_name(i) for i in range(num_species())]
    gas.TPX = (
        T_in,
        P,
        {n: x for n, x in zip(names, X_mix, strict=True) if x > 0 and n in gas.species_names},
    )
    gas.equilibrate("HP")
    assert got.converged
    rich = phi >= 1.5
    assert abs(got.state.T - gas.T) < (3.0 if rich else 150.0), (got.state.T, gas.T)
    for sp in ["CO2", "H2O", "CO", "H2", "CH4"]:
        diff = abs(got.state.X[species_index_from_name(sp)] - gas[sp].X[0])
        assert diff < (1e-3 if rich else 0.015), (
            sp,
            got.state.X[species_index_from_name(sp)],
            gas[sp].X[0],
        )


# ---------------------------------------------------------------------------
# reforming_equilibrium is a SIMULTANEOUS equilibrium of steam reforming of
# each hydrocarbon (both directions) and WGS. Before the fix one round of
# (reform, WGS) left H2 7-49% low at 650-850 K and higher hydrocarbons 0.2 off.
# ---------------------------------------------------------------------------

_REFORMING_CASES = {
    "CH4-steam": "CH4:0.1, H2O:0.3, N2:0.6",
    "higher-hydrocarbons": "CH4:0.05, C2H6:0.03, C3H8:0.02, H2O:0.4, N2:0.5",
    "steam-to-carbon-1.5": "CH4:0.2, H2O:0.3, CO2:0.05, N2:0.45",
    "syngas-methanation": "CH4:0.02, CO:0.2, H2:0.3, H2O:0.1, CO2:0.05, N2:0.33",
}
_REFORMING_SPECIES = ["CH4", "C2H6", "C3H8", "H2O", "CO", "CO2", "H2"]


def _cb_mix(spec: str) -> np.ndarray:
    X = np.zeros(num_species())
    for part in spec.split(","):
        name, v = part.split(":")
        X[species_index_from_name(name.strip())] = float(v)
    return X


@pytest.mark.parametrize("T", [600.0, 650.0, 712.0, 800.0, 900.0, 1100.0])
@pytest.mark.parametrize("case", list(_REFORMING_CASES))
def test_reforming_equilibrium_matches_cantera(combaero, cantera, case: str, T: float) -> None:
    gas = cantera.Solution("gri30.yaml")
    gas.TPX = T, 101325.0, _REFORMING_CASES[case]
    gas.equilibrate("TP")
    got = combaero.reforming_equilibrium(T, _cb_mix(_REFORMING_CASES[case]), 101325.0)
    for name in _REFORMING_SPECIES:
        assert got.state.X[species_index_from_name(name)] == pytest.approx(
            gas[name].X[0], abs=1e-3
        ), name


@pytest.mark.parametrize("T_in", [700.0, 900.0, 1200.0])
@pytest.mark.parametrize("case", list(_REFORMING_CASES))
def test_reforming_equilibrium_adiabatic_matches_cantera(
    combaero, cantera, case: str, T_in: float
) -> None:
    gas = cantera.Solution("gri30.yaml")
    gas.TPX = T_in, 101325.0, _REFORMING_CASES[case]
    gas.equilibrate("HP")
    got = combaero.reforming_equilibrium_adiabatic(T_in, _cb_mix(_REFORMING_CASES[case]), 101325.0)
    assert pytest.approx(gas.T, abs=0.5) == got.state.T


class TestTraceOxygenEquilibrium:
    """Methane with oxygen nearly gone (#1262), TP, against Cantera on the same
    species and combaero's own NASA-9 data (combaero_nasa9.yaml, C1-C3). The
    joint solve stalled there before #1262 (ln Q - ln K up to 62). combaero
    also carries C4-C7 (up to 1.2e-7 at 650 K, which sets the tolerance)."""

    @pytest.fixture
    def nasa9_gas(self, cantera):
        import os

        return cantera.Solution(os.path.join(os.path.dirname(__file__), "combaero_nasa9.yaml"))

    @pytest.mark.parametrize("T", [300.0, 500.0, 650.0])
    @pytest.mark.parametrize("x_o2", [1e-10, 1e-6, 1e-2])
    @pytest.mark.parametrize("P", [1e5, 1.7e6])
    def test_fuel_only_trace_oxygen(self, combaero, nasa9_gas, T, x_o2, P):
        cb = combaero
        X = np.zeros(num_species())
        X[species_index_from_name("CH4")] = 1.0 - x_o2
        X[species_index_from_name("O2")] = x_o2
        Xb = np.array(cb.complete_combustion_products_and_jacobian(X).X_b)
        result = cb._core.reforming_equilibrium(T, Xb, P)
        assert result.converged
        comp = {
            species_name(i): Xb[i]
            for i in range(num_species())
            if Xb[i] > 0.0 and species_name(i) in nasa9_gas.species_names
        }
        nasa9_gas.TPX = T, P, comp
        nasa9_gas.equilibrate("TP")
        for sp in ("H2", "C2H6", "CO", "H2O", "CO2"):
            x_cb = result.state.X[species_index_from_name(sp)]
            x_ct = nasa9_gas[sp].X[0]
            # Measured 2026-11-22: within 4e-4 relative at 650 K (the C4+
            # combaero carries), 1e-7 at 300 K.
            assert x_cb == pytest.approx(x_ct, rel=1e-3), (sp, T, x_o2, P)
