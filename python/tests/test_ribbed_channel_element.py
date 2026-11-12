"""Ribbed channels in ConvectiveSurface: wall weighting and the documented gap.

The correlation gives the RIBBED SIDE. Combining it with the smooth walls is
the element's job, because how many walls are ribbed is a design choice rather
than a property of the correlation.

Friction and heat transfer are weighted alike (#1106):

- the set's `f` is Han's FOUR-SIDED FANNING factor; the channel's drop uses
  the Darcy `4 [phi f + (1 - phi) f_s]`, Han's own area weighting (Han and
  Park 1988 Eq. 6), which is his measured two-wall `fbar` for n = 2
- heat transfer is the same area weighting of ribbed and smooth walls
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network.components import ConvectiveSurface, RibbedModel, SmoothModel

X_AIR = cb.species.dry_air()
STATE = {"T": 600.0, "P": 2.0e5, "X": X_AIR}
FLOW = {"velocity": 60.0, "diameter": 0.025, "length": 0.6, "T_hot": 900.0}


def _surface(**kw) -> ConvectiveSurface:
    return ConvectiveSurface(area=0.1, model=RibbedModel(e_D=0.06, p_e=10.0, **kw))


def _run(surface: ConvectiveSurface):
    return surface.htc_and_T(**STATE, **FLOW)


def test_ribbed_side_and_smooth_side_are_both_reported() -> None:
    """The channel average hides a 20% modelling choice; both halves are visible."""
    r = _run(_surface())
    assert r.h_ribbed > r.h_smooth > 0.0
    assert r.h_smooth < r.h < r.h_ribbed


@pytest.mark.parametrize(("n", "expect"), [(1, "least"), (2, "middle"), (4, "all ribbed")])
def test_wall_count_moves_the_average_monotonically(n: int, expect: str) -> None:
    values = {k: _run(_surface(n_ribbed_walls=k)).h for k in (1, 2, 4)}
    assert values[1] < values[2] < values[4]


def test_four_ribbed_walls_is_exactly_the_ribbed_side() -> None:
    """With no smooth wall left, the weighting must be an identity."""
    r = _run(_surface(n_ribbed_walls=4))
    assert r.h == pytest.approx(r.h_ribbed, rel=1e-12)


@pytest.mark.parametrize("n", [3, 0, 5, -2, 2.5, 3.0, float("nan"), True, "2", None])
def test_an_unrepresentable_wall_count_is_rejected_at_construction(n: object) -> None:
    """Three ribbed walls has no unambiguous geometry, and a fractional, bool
    or non-numeric count names no walls: refused when the model is built,
    not at solve time, naming the allowed values."""
    with pytest.raises(ValueError, match="1, 2 or 4"):
        RibbedModel(e_D=0.06, p_e=10.0, n_ribbed_walls=n)


@pytest.mark.parametrize(("given", "count"), [(1.0, 1), (2.0, 2), (4.0, 4), (2, 2)])
def test_an_integral_float_wall_count_is_the_count_it_names(given: float, count: int) -> None:
    """2.0 (e.g. from JSON) is two walls, stored as an int so the C++ wall
    rule accepts it, and gives the same channel as the int."""
    model = RibbedModel(e_D=0.06, p_e=10.0, n_ribbed_walls=given)
    assert model.n_ribbed_walls == count
    assert type(model.n_ribbed_walls) is int
    as_float = ConvectiveSurface(area=0.1, model=model)
    as_int = _surface(n_ribbed_walls=count)
    assert _run(as_float).dP == _run(as_int).dP


def test_the_documented_gap_and_the_knob_that_closes_it() -> None:
    """Plain smooth walls under-predict Han's own channel average by ~20%.

    Han reports the average directly, implying h_s/h_r ~ 0.70, because ribs
    enhance the adjacent smooth wall by 10-50% as well -- which no correlation
    here covers. Using the plain smooth correlation gives ~0.42.

    This is the regression test for that gap: it pins the size of the
    under-prediction, and that `smooth_wall_Nu_multiplier` is the supported way
    to close it rather than a constant invented inside the library.
    """
    plain = _run(_surface())
    assert plain.h_smooth / plain.h_ribbed == pytest.approx(0.42, abs=0.02)

    tuned = _run(_surface(smooth_wall_Nu_multiplier=1.67))
    assert tuned.h_smooth / tuned.h_ribbed == pytest.approx(0.70, abs=0.02)
    assert tuned.h / plain.h == pytest.approx(1.20, abs=0.02)


GEOM = cb.RibGeometry(e_D=0.06, p_e=10.0, W_H=1.0, alpha_deg=90.0)


def _han_fbar(re: float, n_walls: int, w_h: float = 1.0) -> float:
    """Han's passage-average FANNING friction, written out from the source.

    Independent of `rib_channel_friction`: Han and Park (1988) form the
    four-sided f from the measured two-wall fbar by area weighting (Eq. 6),
    f = fbar + (H/W)(fbar - f_s), f_s = 0.046 Re^-0.2 (Eq. 8). Ribbed walls
    of width W; perimeter fraction W/(W+H) per opposite pair.
    """
    geom = cb.RibGeometry(e_D=0.06, p_e=10.0, W_H=w_h, alpha_deg=90.0)
    f4 = cb.evaluate_rib(cb.han_1988_orthogonal(), geom, re).f
    f_s = 0.046 * re**-0.2
    phi = {1: 0.5 * w_h / (w_h + 1.0), 2: w_h / (w_h + 1.0), 4: 1.0}[n_walls]
    return phi * f4 + (1.0 - phi) * f_s


@pytest.mark.parametrize("n", [1, 2, 4])
def test_channel_friction_is_hans_area_weighting_in_darcy_form(n: int) -> None:
    """f (Darcy) = 4 x Han's wall-weighted Fanning fbar, for every wall count."""
    r = _run(_surface(n_ribbed_walls=n))
    assert r.f == pytest.approx(4.0 * _han_fbar(r.Re, n), rel=1e-9)


def test_two_ribbed_walls_inverts_hans_eq6_exactly() -> None:
    """n = 2 recovers the measured fbar Han converted from (Eq. 6 inverted)."""
    for w_h in (1.0, 2.0, 4.0):
        geom = cb.RibGeometry(e_D=0.06, p_e=10.0, W_H=w_h, alpha_deg=90.0)
        re = 3.0e4
        f4 = cb.evaluate_rib(cb.han_1988_orthogonal(), geom, re).f
        f_s = 0.046 * re**-0.2
        fbar = (f4 * w_h + f_s) / (w_h + 1.0)
        assert f4 == pytest.approx(fbar + (fbar - f_s) / w_h, rel=1e-12)
        f_d, _ = cb.rib_channel_friction(cb.han_1988_orthogonal(), geom, re, 2)
        assert f_d == pytest.approx(
            4.0 * fbar, rel=1e-9
        )  # C++ floors |Re| smoothly: 1e-10 at Re 3e4


def test_four_ribbed_walls_is_the_four_sided_factor_and_re_free() -> None:
    f_d, df = cb.rib_channel_friction(cb.han_1988_orthogonal(), GEOM, 3.0e4, 4)
    assert f_d == pytest.approx(4.0 * cb.evaluate_rib(cb.han_1988_orthogonal(), GEOM, 3.0e4).f)
    assert df == 0.0


def test_drop_matches_hans_two_wall_passage_drop() -> None:
    """The #1106 case: dP = 4 fbar (L/D) rho v^2 / 2 with Han's own fbar.

    The element gave 2884 Pa here before (four-sided Fanning f applied as
    Darcy), about 2.2x below Han's two-wall passage drop.
    """
    r = _run(_surface(n_ribbed_walls=2))
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    v = FLOW["velocity"]
    expect = 4.0 * _han_fbar(r.Re, 2) * FLOW["length"] / FLOW["diameter"] * 0.5 * rho * v * v
    assert r.dP == pytest.approx(expect, rel=1e-9)
    assert r.dP > 5000.0


@pytest.mark.parametrize("re", [-5.0e4, -10.0, 0.0, 10.0, 3.0e4])
@pytest.mark.parametrize("n", [1, 2, 4])
def test_channel_friction_derivative_matches_central_differences(re: float, n: int) -> None:
    """df_D/dRe is the smooth walls' alone; finite and odd through zero flow."""
    s = cb.han_1988_orthogonal()
    f, df = cb.rib_channel_friction(s, GEOM, re, n)
    h = max(1e-6 * abs(re), 1e-3)
    fd = (
        cb.rib_channel_friction(s, GEOM, re + h, n)[0]
        - cb.rib_channel_friction(s, GEOM, re - h, n)[0]
    ) / (2.0 * h)
    assert f > 0.0
    assert df == pytest.approx(fd, rel=1e-6, abs=1e-12)
    assert cb.rib_channel_friction(s, GEOM, -re, n)[1] == pytest.approx(-df, abs=1e-15)


def test_friction_rises_with_the_wall_count() -> None:
    f_values = {n: _run(_surface(n_ribbed_walls=n)).f for n in (1, 2, 4)}
    assert f_values[1] < f_values[2] < f_values[4]


@pytest.mark.parametrize(("n", "flagged"), [(1, True), (2, False), (4, True)])
def test_wall_counts_han_did_not_measure_are_flagged(n: int, flagged: bool) -> None:
    """Every rib set was measured on two opposite ribbed walls; one or four
    walls is the same weighting outside that configuration: reported, not
    refused. The probe point is inside han_1988_orthogonal's box otherwise."""
    s = ConvectiveSurface(area=0.1, model=RibbedModel(e_D=0.06, p_e=10.0, n_ribbed_walls=n))
    flow = dict(FLOW, velocity=15.0)
    r = s.htc_and_T(**STATE, **flow)
    assert 1.0e4 < r.Re < 6.0e4
    assert r.extrapolated is flagged


def test_the_user_knobs_still_apply_on_top() -> None:
    base = _surface()
    scaled = _surface()
    scaled.Nu_multiplier = 2.0
    scaled.f_multiplier = 3.0
    r0, r1 = _run(base), _run(scaled)
    assert r1.h == pytest.approx(2.0 * r0.h, rel=1e-12)
    assert r1.f == pytest.approx(3.0 * r0.f, rel=1e-12)


def test_a_user_supplied_set_is_used_and_G_is_inverse() -> None:
    """A user set replaces the shipped one, and G runs the other way.

    G is the heat-transfer ROUGHNESS function: a larger G means a lower Stanton
    number and so LESS heat transfer. Pinned because the sign is easy to get
    backwards -- this test was written the wrong way round first, and the
    implementation was right.
    """
    worse = cb.han_1988_orthogonal()
    worse.name, worse.source = "rig_3", "measured 2026-09"
    worse.provenance = cb.RibProvenance.User
    worse.C_G = 4.1  # above Han's 3.7

    better = cb.han_1988_orthogonal()
    better.C_G = 3.3  # below it

    default = _run(_surface()).h_ribbed
    assert _run(_surface(correlation_set=worse)).h_ribbed < default
    assert _run(_surface(correlation_set=better)).h_ribbed > default


def test_out_of_range_reports_but_still_returns() -> None:
    """Validity is advisory: the band belongs to the source's rig."""
    r = _run(_surface())
    assert isinstance(r.extrapolated, bool)
    assert r.h > 0.0


def test_smooth_model_is_unaffected() -> None:
    s = ConvectiveSurface(area=0.1, model=SmoothModel())
    assert s.htc_and_T(**STATE, **FLOW).h > 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# ---------------------------------------------------------------------------
# ChannelElement residual path
# ---------------------------------------------------------------------------

from combaero.network.components import ChannelElement, NetworkMixtureState  # noqa: E402

Y_AIR = cb.mole_to_mass(X_AIR)


def _element(model) -> ChannelElement:
    return ChannelElement(
        id="ch",
        from_node="A",
        to_node="B",
        length=0.6,
        diameter=0.025,
        roughness=1e-5,
        surface=ConvectiveSurface(area=0.1, model=model),
    )


def _element_states(m_dot: float):
    a = NetworkMixtureState(P=2e5, Pt=2.1e5, Tt=610.0, m_dot=m_dot, Y=Y_AIR)
    b = NetworkMixtureState(P=1.95e5, Pt=1.95e5, Tt=605.0, m_dot=m_dot, Y=Y_AIR)
    return a, b


def _rib_element() -> ChannelElement:
    return _element(RibbedModel(e_D=0.06, p_e=10.0, W_H=1.0, n_ribbed_walls=2))


@pytest.mark.parametrize("m_dot", [0.05, 0.35, 1.0])
def test_ribbed_jacobian_matches_central_differences(m_dot: float) -> None:
    el = _rib_element()
    _, jac = el.residuals(*_element_states(m_dot))
    h = 1e-7
    fd = (
        el.residuals(*_element_states(m_dot + h))[0][0]
        - el.residuals(*_element_states(m_dot - h))[0][0]
    ) / (2.0 * h)
    assert jac[0]["ch.m_dot"] == pytest.approx(fd, rel=1e-6)


def test_ribbed_drop_follows_the_flow_but_its_derivative_does_not() -> None:
    """The odd-quantity case: magnitude from the correlation, sign from the
    flow. dP is odd in m_dot, so d(dP)/d(m_dot) is EVEN and must not flip at
    zero. Getting that backwards puts a sign error in the Jacobian exactly at
    the crossing."""
    el = _rib_element()
    out = {}
    for m in (0.35, -0.35):
        a, b = _element_states(m)
        res, jac = el.residuals(a, b)
        out[m] = (a.Pt - b.Pt - res[0], jac[0]["ch.m_dot"])

    assert out[0.35][0] == pytest.approx(-out[-0.35][0], rel=1e-12)
    assert out[0.35][1] == pytest.approx(out[-0.35][1], rel=1e-12)


def test_ribbed_drop_exceeds_smooth_at_the_same_flow() -> None:
    a, b = _element_states(0.35)
    rib = a.Pt - b.Pt - _rib_element().residuals(a, b)[0][0]
    smooth = a.Pt - b.Pt - _element(SmoothModel()).residuals(a, b)[0][0]
    assert rib > 2.0 * smooth


def test_friction_multiplier_scales_the_ribbed_drop() -> None:
    plain = _rib_element()
    scaled = _rib_element()
    scaled.surface.f_multiplier = 2.5
    a, b = _element_states(0.35)
    dp0 = a.Pt - b.Pt - plain.residuals(a, b)[0][0]
    dp1 = a.Pt - b.Pt - scaled.residuals(a, b)[0][0]
    assert dp1 == pytest.approx(2.5 * dp0, rel=1e-12)


def test_a_ribbed_channel_solves_in_a_network() -> None:
    """The test that actually matters: does it converge in a real solve?

    Unit tests confirm the residual and its Jacobian in isolation. This
    confirms the element works where it is used -- the solver probes states the
    unit tests never construct, which is the reason the guards in the
    correlation exist at all.

    Ribs restrict flow, so at a fixed pressure drop a ribbed channel must pass
    less than a smooth one.
    """
    from combaero.network import FlowNetwork, NetworkSolver
    from combaero.network.components import PressureBoundary

    def solve(model) -> dict:
        g = FlowNetwork()
        a = PressureBoundary("A")
        a.Pt, a.Tt, a.Y = 2.4e5, 600.0, Y_AIR
        b = PressureBoundary("B")
        b.Pt, b.Tt, b.Y = 2.0e5, 600.0, Y_AIR
        g.add_node(a)
        g.add_node(b)
        g.add_element(_element(model))
        return NetworkSolver(g).solve()

    smooth = solve(SmoothModel())
    ribbed = solve(RibbedModel(e_D=0.06, p_e=10.0, W_H=1.0, n_ribbed_walls=2))

    assert smooth["__success__"], smooth.get("__message__")
    assert ribbed["__success__"], ribbed.get("__message__")
    assert 0.0 < ribbed["ch.m_dot"] < smooth["ch.m_dot"]


def test_ribbed_result_carries_the_wall_coupling_derivatives() -> None:
    """A ribbed channel joined by a ThermalWall must not die mid-solve.

    The solver's wall-coupling path reads `dh_dmdot`, `dh_dT`, `dT_aw_dmdot`
    and `dT_aw_dT` off the channel result. The first version of
    `_RibbedChannelResult` had none of them, and nothing in this file caught it
    -- not one test coupled a wall. The coupled combustor example did, which is
    the argument for examples being executed rather than only read.
    """
    r = _run(_surface())
    for field in ("dh_dmdot", "dh_dT", "dT_aw_dmdot", "dT_aw_dT"):
        assert hasattr(r, field), f"missing {field}"
        assert isinstance(getattr(r, field), float)

    # dh/dmdot must be real, not a placeholder zero: heat transfer rises with
    # flow, so the ribbed side has to contribute.
    assert r.dh_dmdot != 0.0


@pytest.mark.parametrize("area", [0.01, 0.1])
@pytest.mark.parametrize("n_walls", [1, 2, 4])
def test_ribbed_wall_coupling_derivative_matches_finite_difference(
    area: float, n_walls: int
) -> None:
    """dh/dmdot against central differences on the SOLVER's mass flow.

    The existing check above only asked for a nonzero value, which is how two
    defects survived (#456): the mass flow was taken on the convective area,
    and the St d(rho v cp)/dmdot term of h_r = St rho v cp was missing -- the
    analytic value read 2320 against a true 7652 here.
    """
    import math

    s = ConvectiveSurface(area=area, model=RibbedModel(e_D=0.06, p_e=10.0, n_ribbed_walls=n_walls))
    v, dv = FLOW["velocity"], 1e-4
    flow = {k: val for k, val in FLOW.items() if k != "velocity"}
    up = s.htc_and_T(**STATE, **flow, velocity=v + dv).h
    dn = s.htc_and_T(**STATE, **flow, velocity=v - dv).h
    rho = cb.density(STATE["T"], STATE["P"], STATE["X"])
    rho_a = rho * math.pi / 4.0 * FLOW["diameter"] ** 2
    assert _run(s).dh_dmdot == pytest.approx((up - dn) / (2 * dv) / rho_a, rel=1e-4)
