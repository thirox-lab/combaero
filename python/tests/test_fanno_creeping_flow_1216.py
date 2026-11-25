"""Compressible ducts at creeping flow match Hagen-Poiseuille (#1216).

Below kFannoFlowDriveFloor (1 Pa), fanno_channel_flow and
fanno_rayleigh_channel_flow previously scaled the 1 Pa flux linearly to zero
with the odd cubic. On short ducts (L/D = 1..10), 1 Pa is inertia-dominated
(static exit) or turbulent (total exit), so creeping flow was up to 10x below
Hagen-Poiseuille independent of drive.

Viscous friction models now solve the exact low-Mach viscous/inertial closure
directly in drive and hand off C1 (both G and all five Jacobian partials) to
the compressible kernel at kFannoFlowDriveFloor.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import ChannelElement, FlowNetwork, NetworkSolver, PressureBoundary

X_AIR = list(cb.species.dry_air())
Y_AIR = list(cb.mole_to_mass(X_AIR))
PT0 = 1.0e5
TT0 = 300.0
D = 0.01


def _hp_flux(dp: float, length: float, diameter: float = D) -> float:
    rho = cb.density(TT0, PT0, X_AIR)
    mu = cb.viscosity(TT0, PT0, X_AIR)
    return dp * rho * diameter**2 / (32.0 * mu * length)


@pytest.mark.parametrize("length_over_d", [1.0, 10.0, 100.0, 1000.0])
@pytest.mark.parametrize("model", ["haaland", "serghides", "colebrook", "petukhov"])
def test_creeping_flow_matches_hagen_poiseuille_across_l_over_d(
    length_over_d: float, model: str
) -> None:
    length = length_over_d * D
    dp = 1.0e-6
    p_target = PT0 - dp
    dp_actual = PT0 - p_target
    hp = _hp_flux(dp_actual, length)
    # Physical exit head q / dp at L/D = 1, dp = 1e-6 Pa is ~1.6e-4.
    tol_s = 3.0e-4 if length_over_d == 1.0 else 3.0e-5
    f_stat = cb.fanno_channel_flow(PT0, TT0, X_AIR, p_target, False, length, D, 0.0, model)
    f_tot = cb.fanno_channel_flow(PT0, TT0, X_AIR, p_target, True, length, D, 0.0, model)
    r_stat = cb.fanno_rayleigh_channel_flow(
        PT0, TT0, TT0, X_AIR, p_target, False, length, D, 0.0, model
    )
    r_tot = cb.fanno_rayleigh_channel_flow(
        PT0, TT0, TT0, X_AIR, p_target, True, length, D, 0.0, model
    )
    assert f_stat.G / hp == pytest.approx(1.0, rel=tol_s)
    assert f_tot.G / hp == pytest.approx(1.0, rel=3.0e-5)
    assert r_stat.G / hp == pytest.approx(1.0, rel=tol_s)
    assert r_tot.G / hp == pytest.approx(1.0, rel=3.0e-5)


@pytest.mark.parametrize("length", [0.01, 0.1, 1.0])
def test_compressible_and_incompressible_channel_agree_at_creeping_flow(length: float) -> None:
    """Issue #1216 network reproducer: PressureBoundary -> ChannelElement ->
    PressureBoundary at dp = 1e-6 Pa, D = 1 cm."""
    dp = 1.0e-6
    area = 0.25 * math.pi * D**2
    hp_mdot = area * _hp_flux(PT0 - (PT0 - dp), length)
    mdots: dict[str, float] = {}
    for regime in ("incompressible", "compressible"):
        g = FlowNetwork()
        b_in = PressureBoundary("in")
        b_in.Pt, b_in.Tt, b_in.Y = PT0, TT0, Y_AIR
        b_out = PressureBoundary("out")
        b_out.Pt, b_out.Tt, b_out.Y = PT0 - dp, TT0, Y_AIR
        g.add_node(b_in)
        g.add_node(b_out)
        g.add_element(
            ChannelElement(
                "c",
                "in",
                "out",
                length=length,
                diameter=D,
                roughness=0.0,
                friction_model="haaland",
                regime=regime,
            )
        )
        res = NetworkSolver(g).solve()
        assert res["__success__"], f"{regime} L={length}: {res['__message__']}"
        mdots[regime] = float(res["c.m_dot"])
    assert mdots["compressible"] == pytest.approx(mdots["incompressible"], rel=2.0e-4)
    assert mdots["compressible"] == pytest.approx(hp_mdot, rel=3.0e-4)


@pytest.mark.parametrize("heated", [False, True])
@pytest.mark.parametrize("exit_total", [False, True])
def test_value_and_jacobians_continuous_across_drive_floor(heated: bool, exit_total: bool) -> None:
    """Continuity of G and all partials across FANNO_FLOW_DRIVE_FLOOR."""
    pt0 = 2.0e5
    tt0 = 400.0
    tte = 550.0
    length = 0.05
    floor = cb.FANNO_FLOW_DRIVE_FLOOR
    eps = 1.0e-7 * floor
    v_o2 = [float(i == cb.species_index_from_name("O2")) for i in range(cb.num_species())]

    def call(drive: float) -> cb.FannoChannelFlow:
        if heated:
            return cb.fanno_rayleigh_channel_flow(
                pt0,
                tt0,
                tte,
                X_AIR,
                pt0 - drive,
                exit_total,
                length,
                D,
                0.0,
                "haaland",
                1.0,
                -1.0,
                [v_o2],
                True,
            )
        return cb.fanno_channel_flow(
            pt0,
            tt0,
            X_AIR,
            pt0 - drive,
            exit_total,
            length,
            D,
            0.0,
            "haaland",
            1.0,
            -1.0,
            [v_o2],
            True,
        )

    lo = call(floor - eps)
    at = call(floor)
    hi = call(floor + eps)
    assert pytest.approx(hi.G, rel=1.0e-6) == lo.G
    assert lo.dG_dP_target == pytest.approx(hi.dG_dP_target, rel=1.0e-5)
    assert lo.dG_dPt0 == pytest.approx(hi.dG_dPt0, rel=1.0e-5)
    assert lo.dG_dTt0 == pytest.approx(hi.dG_dTt0, rel=1.0e-4)
    if heated:
        assert lo.dG_dTt_exit == pytest.approx(hi.dG_dTt_exit, rel=1.0e-4)
    assert lo.dG_ddir[0] == pytest.approx(at.dG_ddir[0], rel=1.0e-6)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "Known limit: adiabatic static-exit raw_flow at 1 Pa carries pre-existing "
        "d_ddir low-Mach cancellation noise (+8.5e-5 / 1.28x at 1.000 Pa on "
        "Pt0=2e5, Tt0=400 K, L=0.5 m, D=0.01 m), which the C1 sub-floor handoff "
        "inherits proportionally (analytic/FD = 0.967, 0.923, 1.93 at 0.3, 0.5, 0.9 Pa)."
    ),
)
@pytest.mark.parametrize("drive", [0.3, 0.5, 0.9])
def test_adiabatic_static_exit_subfloor_dY_known_raw_flow_limit(drive: float) -> None:
    """Pins the known adiabatic static-exit sub-floor dY limit inherited from
    raw_flow's 1 Pa d_ddir cancellation noise (will XPASS-fail once fixed)."""
    pt0 = 2.0e5
    tt0 = 400.0
    length = 0.5
    k_o2 = cb.species_index_from_name("O2")
    v_o2 = [float(i == k_o2) for i in range(cb.num_species())]
    pb = pt0 - drive
    r = cb.fanno_channel_flow(
        pt0, tt0, X_AIR, pb, False, length, D, 0.0, "haaland", 1.0, -1.0, [v_o2], True
    )
    hy = 1.0e-3
    yp = list(Y_AIR)
    ym = list(Y_AIR)
    yp[k_o2] += hy
    ym[k_o2] -= hy
    gp = cb.fanno_channel_flow(
        pt0, tt0, cb.mass_to_mole(yp), pb, False, length, D, 0.0, "haaland", 1.0, -1.0, [], False
    ).G
    gm = cb.fanno_channel_flow(
        pt0, tt0, cb.mass_to_mole(ym), pb, False, length, D, 0.0, "haaland", 1.0, -1.0, [], False
    ).G
    fd_y = (gp - gm) / (2.0 * hy)
    assert r.dG_ddir[0] / fd_y == pytest.approx(1.0, rel=1.0e-2)
