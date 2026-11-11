"""Tests for combustor pressure-loss correlation classes."""

import pytest

from combaero.network import (
    ConstantFractionLoss,
    LinearThetaFractionLoss,
)

# ---------------------------------------------------------------------------
# Minimal duck-typed context for pure-Python tests (no C++ required)
# ---------------------------------------------------------------------------


class _Ctx:
    def __init__(self, theta: float = 0.0) -> None:
        self.theta = theta


# ---------------------------------------------------------------------------
# ConstantFractionLoss
# ---------------------------------------------------------------------------


class TestConstantFractionLoss:
    def test_returns_xi_and_zero_derivative(self) -> None:
        loss = ConstantFractionLoss(xi=0.03)
        xi, dxi_dtheta = loss(_Ctx(theta=0.0))
        assert xi == pytest.approx(0.03)
        assert dxi_dtheta == pytest.approx(0.0)

    def test_independent_of_theta(self) -> None:
        loss = ConstantFractionLoss(xi=0.05)
        for theta in [-1.0, 0.0, 0.5, 2.0, 10.0]:
            xi, dxi = loss(_Ctx(theta=theta))
            assert xi == pytest.approx(0.05)
            assert dxi == pytest.approx(0.0)

    def test_default_xi(self) -> None:
        loss = ConstantFractionLoss()
        xi, _ = loss(_Ctx())
        assert xi == pytest.approx(0.03)

    def test_zero_loss(self) -> None:
        loss = ConstantFractionLoss(xi=0.0)
        xi, dxi = loss(_Ctx(theta=1.5))
        assert xi == pytest.approx(0.0)
        assert dxi == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# LinearThetaFractionLoss
# ---------------------------------------------------------------------------


class TestLinearThetaFractionLoss:
    def test_at_zero_theta(self) -> None:
        loss = LinearThetaFractionLoss(k=0.5, xi0=0.02)
        xi, dxi = loss(_Ctx(theta=0.0))
        assert xi == pytest.approx(0.02)
        assert dxi == pytest.approx(0.5)

    def test_linear_slope(self) -> None:
        loss = LinearThetaFractionLoss(k=0.4, xi0=0.01)
        theta = 2.0
        xi, dxi = loss(_Ctx(theta=theta))
        assert xi == pytest.approx(0.4 * theta + 0.01)
        assert dxi == pytest.approx(0.4)

    def test_derivative_constant_wrt_theta(self) -> None:
        loss = LinearThetaFractionLoss(k=0.6, xi0=0.03)
        for theta in [0.0, 0.5, 1.0, 3.0]:
            _, dxi = loss(_Ctx(theta=theta))
            assert dxi == pytest.approx(0.6)

    def test_default_values(self) -> None:
        loss = LinearThetaFractionLoss()
        xi, dxi = loss(_Ctx(theta=1.0))
        assert xi == pytest.approx(0.5 * 1.0 + 0.02)
        assert dxi == pytest.approx(0.5)

    def test_fd_consistency(self) -> None:
        """Analytical derivative matches finite difference."""
        loss = LinearThetaFractionLoss(k=0.35, xi0=0.015)
        theta = 1.5
        eps = 1e-6
        xi0_val, dxi = loss(_Ctx(theta=theta))
        xi_p, _ = loss(_Ctx(theta=theta + eps))
        fd = (xi_p - xi0_val) / eps
        assert dxi == pytest.approx(fd, rel=1e-5)


# ---------------------------------------------------------------------------
# Integration test: PressureLossElement residual enforces correct ratio
# ---------------------------------------------------------------------------


class TestPressureLossElementResidual:
    """Unit-level checks on ``PressureLossElement.residuals`` for each correlation.

    Verifies the element residual ``Pt_out - Pt_in * (1 - xi) = 0``
    vanishes at the physically correct P ratio, for every supported correlation.
    """

    def _make_air_state(self, Pt: float = 150000.0, T: float = 700.0, m_dot: float = 1.0):
        """Helper: build an air-composition NetworkMixtureState."""
        import combaero as cb
        from combaero.network import NetworkMixtureState

        ns = cb.species.num_species
        Y = [0.0] * ns
        for k in range(ns):
            name = cb.species_name(k)
            if name == "N2":
                Y[k] = 0.767
            elif name == "O2":
                Y[k] = 0.233
        return NetworkMixtureState(P=Pt, Pt=Pt, Tt=T, m_dot=m_dot, Y=Y)

    def test_constant_fraction_residual_zero_at_exact_ratio(self) -> None:
        """res == 0 iff Pt_out / Pt_in == (1 - xi)."""
        from combaero.network import PressureLossElement

        xi = 0.05
        P_in = 150000.0
        state_in = self._make_air_state(Pt=P_in)
        state_out = self._make_air_state(Pt=P_in * (1.0 - xi))

        elem = PressureLossElement("loss", "a", "b", correlation=ConstantFractionLoss(xi=xi))
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)

    def test_constant_fraction_residual_nonzero_off_ratio(self) -> None:
        from combaero.network import PressureLossElement

        xi = 0.05
        P_in = 150000.0
        state_in = self._make_air_state(Pt=P_in)
        state_out = self._make_air_state(Pt=P_in)  # no drop

        elem = PressureLossElement("loss", "a", "b", correlation=ConstantFractionLoss(xi=xi))
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(P_in * xi, rel=1e-6)

    def test_linear_theta_residual_uses_cold_flow_without_source(self) -> None:
        """Without a theta source, linear-theta correlation degrades to xi0."""
        import warnings as _w

        from combaero.network import PressureLossElement

        k, xi0 = 0.5, 0.02
        P_in = 150000.0
        state_in = self._make_air_state(Pt=P_in)
        state_out = self._make_air_state(Pt=P_in * (1.0 - xi0))

        elem = PressureLossElement(
            "loss", "a", "b", correlation=LinearThetaFractionLoss(k=k, xi0=xi0)
        )
        # residuals() itself does not warn (warning fires in resolve_topology);
        # here we bypass resolve_topology and hit the no-graph path.
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)

    def test_constant_head_residual_xi_from_dynamic_head(self) -> None:
        """ConstantHeadLoss: xi = zeta * 0.5 * rho * v^2 / P_in."""
        from combaero.network import ConstantHeadLoss, PressureLossElement

        zeta = 5.0
        area = 0.1
        P_in = 150000.0
        mdot = 1.0
        state_in = self._make_air_state(Pt=P_in, m_dot=mdot)

        corr = ConstantHeadLoss(zeta=zeta, area=area)
        elem = PressureLossElement("loss", "a", "b", correlation=corr)
        # Evaluate xi manually via the same ctx builder used by the element
        ctx = elem._build_ctx(state_in, state_in, None)
        xi, _ = corr(ctx)
        assert 0.0 < xi < 1.0

        state_out = self._make_air_state(Pt=P_in * (1.0 - xi), m_dot=mdot)
        res, _ = elem.residuals(state_in, state_out)
        assert res[0] == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# PressureLossElement convective surface + heat transfer
# ---------------------------------------------------------------------------


def _air_mixture_state(P: float = 150000.0, T: float = 700.0, m_dot: float = 1.0):
    """Helper: 76.7/23.3 N2/O2 air NetworkMixtureState at requested conditions."""
    import combaero as cb
    from combaero.network import NetworkMixtureState

    ns = cb.species.num_species
    Y = [0.0] * ns
    for k in range(ns):
        name = cb.species_name(k)
        if name == "N2":
            Y[k] = 0.767
        elif name == "O2":
            Y[k] = 0.233
    return NetworkMixtureState(P=P, Pt=P, Tt=T, m_dot=m_dot, Y=Y)


class TestPressureLossElementConvectiveSurface:
    """Covers the convective-surface extension of :class:`PressureLossElement`."""

    def test_default_has_no_convective_surface(self) -> None:
        """Without a surface kwarg, element remains thermally passive."""
        from combaero.network import PressureLossElement

        elem = PressureLossElement("loss", "a", "b", correlation=ConstantFractionLoss(xi=0.03))
        assert elem.has_convective_surface is False
        assert elem.htc_and_T(_air_mixture_state()) is None

    def test_surface_with_zero_area_disabled(self) -> None:
        """A surface with ``area == 0`` must behave as disabled."""
        from combaero.heat_transfer import ConvectiveSurface, SmoothModel
        from combaero.network import PressureLossElement

        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantFractionLoss(xi=0.03),
            area=0.01,
            surface=ConvectiveSurface(area=0.0, model=SmoothModel()),
        )
        assert elem.has_convective_surface is False
        assert elem.htc_and_T(_air_mixture_state()) is None

    def test_surface_with_area_enables_htc(self) -> None:
        """With ``surface.area > 0`` and a flow area, HTC must be positive and finite."""
        from combaero.heat_transfer import ConvectiveSurface, SmoothModel
        from combaero.network import PressureLossElement

        flow_area = 0.01  # [m^2]
        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantFractionLoss(xi=0.03),
            area=flow_area,
            surface=ConvectiveSurface(area=flow_area, model=SmoothModel()),
        )
        assert elem.has_convective_surface is True

        res = elem.htc_and_T(_air_mixture_state(P=150000.0, T=700.0, m_dot=1.0))
        assert res is not None
        # Sanity ranges: forced convection in a short hot duct at 1 kg/s.
        assert res.h > 0.0 and res.h < 1e6
        assert res.Nu > 0.0
        # T_aw close to recovery temperature (~= static T for Pr~=0.7 at low Mach).
        assert 600.0 < res.T_aw < 900.0

    def test_htc_none_when_flow_area_missing(self) -> None:
        """Flow area unset => velocity undefined => htc_and_T returns None."""
        from combaero.heat_transfer import ConvectiveSurface, SmoothModel
        from combaero.network import PressureLossElement

        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantFractionLoss(xi=0.03),
            area=None,
            surface=ConvectiveSurface(area=0.02, model=SmoothModel()),
        )
        # has_convective_surface is True (surface.area > 0) but htc returns None
        # because flow area is missing. Intentional degrade-gracefully behavior.
        assert elem.has_convective_surface is True
        assert elem.htc_and_T(_air_mixture_state()) is None


# ---------------------------------------------------------------------------
# PressureLossElement diagnostics: verify enriched key set
# ---------------------------------------------------------------------------


class TestPressureLossElementDiagnostics:
    """Locks in the key set exposed via :meth:`diagnostics`.

    Protects the GUI ``Live Telemetry`` panel (which discovers fields through
    ``QUANTITY_CATALOGUE``) from silent regressions.
    """

    @staticmethod
    def _expected_common_keys() -> set[str]:
        return {
            # Loss / pressure
            "xi",
            "theta",
            "dPt",
            "dP",
            "P_in",
            "P_out",
            "T_in",
            "T_out",
            "Pt_in",
            "Pt_out",
            "Tt_in",
            "Tt_out",
            "mach_in",
            "mach_out",
            "pr_total",
            "pr_static",
            "Re",
            "velocity",
            "ref_location",
            # Thermo + transport at inlet
            "rho_ref",
            "h_ref",
            "s_ref",
            "u_ref",
            "gamma_ref",
            "a_ref",
            "cp_ref",
            "cv_ref",
            "mw_ref",
            "mu_ref",
            "k_ref",
            "Pr_ref",
            "nu_ref",
        }

    def test_constant_fraction_exposes_all_common_keys(self) -> None:
        from combaero.network import PressureLossElement

        state_in = _air_mixture_state(P=150000.0, T=700.0, m_dot=1.0)
        state_out = _air_mixture_state(P=150000.0 * 0.97, T=700.0, m_dot=1.0)
        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantFractionLoss(xi=0.03),
            area=0.01,
        )
        diag = elem.diagnostics(state_in, state_out)
        missing = self._expected_common_keys() - diag.keys()
        assert not missing, f"Missing diagnostic keys: {sorted(missing)}"
        # All values finite floats (except ref_location)
        for k, v in diag.items():
            if k in ("ref_location", "area_source"):
                assert isinstance(v, str)
            else:
                assert isinstance(v, float), f"{k} is not float: {v!r}"

    def test_head_loss_exposes_common_keys_including_htc(self) -> None:
        """With a convective surface, diagnostics must also include Nu/htc/T_aw/f."""
        from combaero.heat_transfer import ConvectiveSurface, SmoothModel
        from combaero.network import ConstantHeadLoss, PressureLossElement

        flow_area = 0.01
        state_in = _air_mixture_state(P=150000.0, T=700.0, m_dot=1.0)
        state_out = _air_mixture_state(P=140000.0, T=700.0, m_dot=1.0)
        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantHeadLoss(zeta=5.0, area=flow_area),
            area=flow_area,
            surface=ConvectiveSurface(area=flow_area, model=SmoothModel()),
        )
        diag = elem.diagnostics(state_in, state_out)
        missing = self._expected_common_keys() - diag.keys()
        assert not missing, f"Missing diagnostic keys: {sorted(missing)}"
        for key in ("Nu", "htc", "T_aw", "f"):
            assert key in diag, f"Expected HTC key {key!r} missing"
            assert diag[key] >= 0.0

    def test_diagnostics_no_htc_without_surface(self) -> None:
        """Sanity: Nu/htc/T_aw/f must be absent when no surface attached."""
        from combaero.network import PressureLossElement

        state = _air_mixture_state()
        elem = PressureLossElement(
            "loss",
            "a",
            "b",
            correlation=ConstantFractionLoss(xi=0.03),
            area=0.01,
        )
        diag = elem.diagnostics(state, state)
        for key in ("Nu", "htc", "T_aw", "f"):
            assert key not in diag


# ---------------------------------------------------------------------------
# End-to-end wall coupling with PressureLossElement on both sides
# ---------------------------------------------------------------------------


class TestPressureLossElementWallCoupling:
    """Guards that two wall-coupled PressureLossElements exchange heat.

    Was the root symptom behind the session: discrete loss could be connected
    via a thermal wall, but ``has_convective_surface`` defaulted to False so
    the wall silently no-op'd. This regression locks the fix in place.
    """

    def test_hot_and_cold_discrete_losses_exchange_heat(self) -> None:
        import numpy as np

        import combaero as cb
        from combaero.heat_transfer import ConvectiveSurface, SmoothModel
        from combaero.network import (
            ConstantFractionLoss,
            FlowNetwork,
            MassFlowBoundary,
            NetworkSolver,
            PlenumNode,
            PressureBoundary,
            PressureLossElement,
            ThermalWall,
            WallLayer,
        )

        Y_air = cb.species.dry_air_mass()

        # Hot side: inlet -> loss -> outlet
        hot_in = MassFlowBoundary("hot_in", m_dot=0.1, Tt=800.0, Y=Y_air)
        hot_out = PressureBoundary("hot_out", Pt=2e5, Tt=800.0, Y=Y_air)
        # Cold side: inlet -> loss -> outlet
        cold_in = MassFlowBoundary("cold_in", m_dot=0.05, Tt=400.0, Y=Y_air)
        cold_out = PressureBoundary("cold_out", Pt=2e5, Tt=400.0, Y=Y_air)

        # Intermediate plenums allow the solver to settle Pt freely.
        hot_mid = PlenumNode("hot_mid")
        cold_mid = PlenumNode("cold_mid")

        area = float(np.pi * 0.04**2 / 4.0)  # 0.04 m Dh
        surf_area = float(np.pi * 0.04 * 0.04)  # A_conv = pi*D*L with L=D

        hot_loss = PressureLossElement(
            "hot_loss",
            from_node="hot_in",
            to_node="hot_mid",
            correlation=ConstantFractionLoss(xi=0.02),
            area=area,
            surface=ConvectiveSurface(area=surf_area, model=SmoothModel()),
        )
        cold_loss = PressureLossElement(
            "cold_loss",
            from_node="cold_in",
            to_node="cold_mid",
            correlation=ConstantFractionLoss(xi=0.02),
            area=area,
            surface=ConvectiveSurface(area=surf_area, model=SmoothModel()),
        )

        # Outlet connections (lossless link from mid plenum to pressure boundary)
        from combaero.network import LosslessConnectionElement

        hot_exit = LosslessConnectionElement("hot_exit", "hot_mid", "hot_out")
        cold_exit = LosslessConnectionElement("cold_exit", "cold_mid", "cold_out")

        net = FlowNetwork()
        for n in (hot_in, hot_mid, hot_out, cold_in, cold_mid, cold_out):
            net.add_node(n)
        for e in (hot_loss, cold_loss, hot_exit, cold_exit):
            net.add_element(e)

        wall = ThermalWall(
            id="coupling_wall",
            element_a="hot_loss",
            element_b="cold_loss",
            layers=[WallLayer(thickness=0.002, conductivity=25.0)],
        )
        net.add_wall(wall)

        net.thermal_coupling_enabled = True
        solver = NetworkSolver(net)
        sol = solver.solve(method="hybr", options={"xtol": 1e-8})

        assert sol["__success__"], f"Did not converge: |F|={sol['__final_norm__']:.2e}"

        # Wall must carry non-trivial heat.
        Q = abs(float(sol["coupling_wall.Q"]))
        assert Q > 1.0, f"Expected non-negligible Q through wall, got {Q:.3e} W"

        # Hot-side plenum should cool, cold-side should heat.
        T_hot_mid = solver._derived_states["hot_mid"][0]
        T_cold_mid = solver._derived_states["cold_mid"][0]
        assert T_hot_mid < 800.0, f"Hot-side should cool: T_hot_mid={T_hot_mid:.1f} K"
        assert T_cold_mid > 400.0, f"Cold-side should heat: T_cold_mid={T_cold_mid:.1f} K"


def test_builtin_losses_evaluate_in_cpp_with_their_partials() -> None:
    """evaluate() is the C++ value with every partial (#718), composed by
    pressure_loss_compose (#725); calling the class gives the same (xi,
    dxi/dtheta). A context with only X (no Y) is converted, so a head loss
    reads the same density either way."""
    from types import SimpleNamespace

    import combaero as cb
    from combaero.network.pressure_loss import (
        ConstantFractionLoss,
        ConstantHeadLoss,
        LinearThetaFractionLoss,
        LinearThetaHeadLoss,
    )

    Y = list(cb.species.dry_air_mass())
    ctx = SimpleNamespace(
        theta=1.2, mdot_total=0.4, Tt_in=700.0, Pt_in=2e5, X_in=list(cb.mass_to_mole(Y))
    )
    for corr in (
        ConstantFractionLoss(xi=0.03, m_ref=0.5),
        LinearThetaFractionLoss(k=0.01, xi0=0.02, m_ref=0.5),
        ConstantHeadLoss(zeta=4.0, area=0.05),
        LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05),
    ):
        r = corr.evaluate(ctx, Y)
        assert isinstance(r, cb.PressureLossResult)
        assert r.dxi_dtheta == (r.dxi_dx[cb.PRESSURE_LOSS_INPUT_THETA] if r.dxi_dx else 0.0)
        assert corr(ctx) == pytest.approx((r.xi, r.dxi_dtheta), rel=1e-12)
        assert corr.evaluate(ctx).xi == pytest.approx(r.xi, rel=1e-12)
        # Below the design flow / for a head loss the flow partial is live.
        assert r.dxi_dmdot > 0.0


def test_the_head_loss_is_total_total() -> None:
    """Lefebvre's pressure-loss factor is total-total at low Mach (#782
    audit): the element evaluates the head loss at the inlet's (Tt, Pt). It
    used the static P with Tt, +8% in the loss at an inlet Mach of 0.3."""

    import combaero as cb
    from combaero.network.components import NetworkMixtureState, PressureLossElement
    from combaero.network.pressure_loss import ConstantHeadLoss

    Y = list(cb.species.dry_air_mass())
    corr = ConstantHeadLoss(zeta=4.0, area=0.05)
    el = PressureLossElement("loss", "a", "b", correlation=corr)
    # A moving inlet: static P well below Pt.
    a = NetworkMixtureState(P=1.6e5, Pt=2.0e5, Tt=700.0, m_dot=8.0, Y=Y)
    b = NetworkMixtureState(P=1.9e5, Pt=1.9e5, Tt=700.0, m_dot=8.0, Y=Y)
    res, _ = el.residuals(a, b)
    expected = cb.combustor_head_loss(0.0, 8.0, 700.0, 2.0e5, Y, 0.0, 4.0, 0.05).xi
    assert res[0] == pytest.approx(b.Pt - a.Pt * (1.0 - expected), rel=1e-12)


def test_a_0_7_context_is_refused_naming_the_new_fields() -> None:
    """A hand-built context with the 0.7 names T_in/P_in (P_in was the
    static P) was evaluated with a DeprecationWarning in 0.8; since 0.9.0 it
    is refused with a message naming Tt_in/Pt_in."""
    from types import SimpleNamespace

    import combaero as cb
    from combaero.network.pressure_loss import ConstantHeadLoss, PressureLossContext

    Y = list(cb.species.dry_air_mass())
    X = list(cb.mass_to_mole(Y))
    corr = ConstantHeadLoss(zeta=4.0, area=0.05)
    corr.evaluate(SimpleNamespace(theta=0.0, mdot_total=0.4, Tt_in=700.0, Pt_in=2e5, X_in=X), Y)
    old_ctx = SimpleNamespace(theta=0.0, mdot_total=0.4, T_in=700.0, P_in=2e5, X_in=X)
    with pytest.raises(AttributeError, match="Tt_in and Pt_in"):
        corr.evaluate(old_ctx, Y)
    ctx = PressureLossContext(theta=0.0, mdot_total=0.4, Tt_in=700.0, Pt_in=2e5, X_in=X)
    with pytest.raises(AttributeError):
        ctx.P_in  # noqa: B018
