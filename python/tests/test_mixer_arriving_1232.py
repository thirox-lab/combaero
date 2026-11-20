"""A mixer mixes the streams that arrive (#1232).

``mixer_from_streams_and_jacobians`` used to form Y, h and Pt as the signed
mean sum_i m_i q_i / sum_i m_i over every stream. With streams of both signs
that mean leaves the range of its inputs and has a pole at zero net flow
(Tt -1700 K at -0.02 kg/s net, 1e10 K near zero, a negative mean Y clamped
so its partials no longer matched). Now the mix is the mean over the
arriving streams (m_i > 0), or over every stream when none arrives (the
#1224 reverse mean). Over streams of one sign the mean is a convex
combination: Y stays in the simplex and h and Pt within the inputs at any
nonzero total, so there is no pole and no near-zero band where stream 0's
state stood in. Only a mixer with no flow at all (every m_i = 0) falls back
to stream 0, with zero partials.
"""

from __future__ import annotations

import random

import pytest
from test_mixer_net_reverse_1224 import _Y, AIR, N, _analytic, _check_partials, _mix, _out

import combaero as cb

FUEL = _Y({"CH4": 1.0})


def test_the_issue_case_is_the_arriving_stream() -> None:
    # In: +0.1 kg/s air at 800 K, 2.2e5 Pa. Out: -0.12 kg/s at 400 K, 2e5 Pa.
    # The signed mean gave Tt -1700 K and Pt 1.0e5 Pa.
    s = [cb.MassStream(0.1, 800.0, 2.2e5, AIR), cb.MassStream(-0.12, 400.0, 2.0e5, AIR)]
    r = _mix(s)
    assert r.Tt_mix == pytest.approx(800.0, rel=1e-12)
    assert r.Pt_mix == 2.2e5
    # A leaving stream moves nothing; one arriving stream's own state is
    # exact, so its flow partials are exactly 0 too.
    assert all(v == 0.0 for i in (0, 1) for v in _analytic(r, i, "mdot"))
    assert all(v == 0.0 for v in _analytic(r, 1, "Tt"))
    _check_partials(s)


@pytest.mark.parametrize("net", [1e-9, -1e-9, 1e-13, -1e-13, 0.0])
def test_no_pole_at_zero_net_flow(net: float) -> None:
    # The signed mean gave 1.18e10 K at -1e-9 and -2.1e10 K at +1e-9.
    s = [cb.MassStream(0.1, 800.0, 2.2e5, AIR), cb.MassStream(net - 0.1, 400.0, 2.0e5, AIR)]
    r = _mix(s)
    assert r.Tt_mix == pytest.approx(800.0, rel=1e-12)
    assert r.Pt_mix == 2.2e5


def _random_streams(rng: random.Random, n: int) -> list:
    names = ["N2", "O2", "CO2", "H2O", "AR", "CH4", "H2", "CO"]
    out = []
    for _ in range(n):
        picks = rng.sample(names, rng.randint(1, 3))
        w = [rng.uniform(0.1, 1.0) for _ in picks]
        Y = _Y({p: v / sum(w) for p, v in zip(picks, w, strict=True)})
        m = rng.uniform(0.005, 0.1) * rng.choice([1.0, -1.0])
        out.append([m, rng.uniform(300.0, 1500.0), rng.uniform(1e5, 5e5), Y])
    return out


def test_mixed_signs_stay_within_the_arriving_streams() -> None:
    # The #1230 review's random mixers: 2-4 streams of both signs and up to
    # three gases, at a net flow down to 1e-11. The signed mean raised in
    # calc_T_from_h (17 of 20,000) and clamped a negative Y (116 of 150).
    rng = random.Random(7)
    n_mixed = 0
    for _ in range(2000):
        ss = _random_streams(rng, rng.randint(2, 4))
        net = rng.choice([1.0, -1.0]) * 10 ** rng.uniform(-11, -1)
        ss[-1][0] += net - sum(t[0] for t in ss)
        ms = [t[0] for t in ss]
        if not (min(ms) < 0.0 < max(ms)):
            continue
        n_mixed += 1
        r = _mix([cb.MassStream(*t) for t in ss])
        arr = [t for t in ss if t[0] > 0.0]
        assert min(t[1] for t in arr) - 1e-6 <= r.Tt_mix <= max(t[1] for t in arr) + 1e-6
        assert (
            min(t[2] for t in arr) * (1 - 1e-14) <= r.Pt_mix <= max(t[2] for t in arr) * (1 + 1e-14)
        )
        m_arr = sum(t[0] for t in arr)
        Y = [sum(t[0] * t[3][k] for t in arr) / m_arr for k in range(N)]
        assert r.Y_mix == pytest.approx(Y, abs=1e-15)
        assert min(r.Y_mix) >= 0.0
    assert n_mixed > 1500


@pytest.mark.parametrize(("Q", "fraction"), [(0.0, 0.0), (250.0, 0.0), (0.0, 0.15), (-80.0, 0.1)])
def test_mixed_sign_multi_species_partials_match(Q: float, fraction: float) -> None:
    # Two arriving gases, two leaving, more out than in.
    s = [
        cb.MassStream(-0.03, 450.0, 2.0e5, _Y({"N2": 0.7, "O2": 0.2, "CO2": 0.1})),
        cb.MassStream(0.01, 800.0, 2.3e5, _Y({"N2": 0.75, "O2": 0.2, "CO2": 0.05})),
        cb.MassStream(-0.02, 600.0, 2.1e5, _Y({"N2": 0.7, "O2": 0.15, "H2O": 0.15})),
        cb.MassStream(0.015, 500.0, 2.2e5, _Y({"N2": 0.6, "O2": 0.1, "H2O": 0.3})),
    ]
    ref = _mix([s[1], s[3]], Q, fraction)
    r = _mix(s, Q, fraction)
    assert _out(r) == _out(ref)
    for i in (0, 2):
        for var in ("mdot", "Tt", "Pt"):
            assert all(v == 0.0 for v in _analytic(r, i, var))
        for j in range(N):
            assert all(v == 0.0 for v in _analytic(r, i, "Y", j))
    _check_partials(s, Q=Q, fraction=fraction)


@pytest.mark.parametrize("Q", [0.0, 50.0])
def test_an_arriving_stream_crossing_zero_is_continuous(Q: float) -> None:
    # Another stream still arrives: the crossing one enters with weight 0.
    def streams(m: float) -> list:
        return [
            cb.MassStream(m, 900.0, 2.4e5, _Y({"N2": 0.75, "O2": 0.15, "H2O": 0.1})),
            cb.MassStream(0.02, 500.0, 2.0e5, AIR),
            cb.MassStream(-0.025, 650.0, 2.1e5, FUEL),
        ]

    # Within the slope times the step either side (|dq/dm| d < 1e-8 here).
    at0 = _out(_mix(streams(0.0), Q))
    for d in (1e-12, -1e-12):
        assert _out(_mix(streams(d), Q)) == pytest.approx(at0, rel=1e-10, abs=1e-10)
    _check_partials(streams(0.004), Q=Q)
    _check_partials(streams(-0.004), Q=Q)


def test_the_last_arriving_stream_crossing_zero_switches_to_the_reverse_mean() -> None:
    # With nothing arriving the mix is the leaving streams' mean (#1224). As
    # the last arriving stream crosses zero the state goes from its own state
    # to that mean: a jump, as at a node whose last inflow vanishes. No
    # network in the stress sweep, the tee sweep, the junction runner or the
    # test suite reaches a mixer with streams of both signs (#1232 audit).
    def streams(m: float) -> list:
        return [
            cb.MassStream(m, 900.0, 2.4e5, AIR),
            cb.MassStream(-0.02, 500.0, 2.0e5, AIR),
            cb.MassStream(-0.01, 800.0, 2.1e5, AIR),
        ]

    up, down = _mix(streams(1e-9)), _mix(streams(-1e-9))
    assert up.Tt_mix == pytest.approx(900.0, rel=1e-12) and up.Pt_mix == 2.4e5
    rev = _mix(streams(-1e-9)[1:])
    assert down.Pt_mix == pytest.approx((0.02 * 2e5 + 0.01 * 2.1e5 + 1e-9 * 2.4e5) / (0.03 + 1e-9))
    assert down.Tt_mix == pytest.approx(rev.Tt_mix, rel=1e-6)
    # Either side, the partials are the state's own.
    _check_partials(streams(0.003))
    _check_partials(streams(-0.003))


@pytest.mark.parametrize("scale", [1e-11, 1e-13, 1e-16, -1e-13, -1e-16])
def test_a_vanishing_flow_keeps_the_same_state(scale: float) -> None:
    # The old band |m| <= 1e-12 took stream 0's Y with the averaged h for
    # 0 < m <= 1e-12 (CH4 as stream 0: Tt 1360 K against 300 K just below
    # zero and 664 K outside the band). The mean of one-sign streams does not
    # depend on their common scale, so h, Y and Pt hold together down to zero.
    def streams(k: float) -> list:
        return [
            cb.MassStream(0.01 * k, 300.0, 2e5, FUEL),
            cb.MassStream(0.09 * k, 700.0, 2.2e5, AIR),
        ]

    ref, r = _mix(streams(1.0 if scale > 0 else -1.0)), _mix(streams(scale))
    assert _out(r) == pytest.approx(_out(ref), rel=1e-12, abs=1e-15)
    _check_partials(streams(scale), rel=1e-4)


def test_no_flow_is_stream_zero_with_zero_partials() -> None:
    s = [cb.MassStream(0.0, 500.0, 2e5, AIR), cb.MassStream(0.0, 900.0, 2.4e5, FUEL)]
    r = _mix(s)
    assert r.Tt_mix == pytest.approx(500.0, rel=1e-12)
    assert r.Pt_mix == 2e5
    assert r.Y_mix == pytest.approx(AIR, abs=1e-15)
    for i in range(2):
        for var in ("mdot", "Tt", "Pt"):
            assert all(v == 0.0 for v in _analytic(r, i, var))


@pytest.mark.parametrize(("a", "b", "w"), [("AR", "CO", 0.5), ("CH4", "CO", 0.7)])
def test_t_from_h_converges_where_tol_is_below_t_resolution(a: str, b: str, w: float) -> None:
    # A mixer still reaches h of order 1e14 J/kg through its heat: Q phi is up
    # to 1.05 Q / 1e-6 kg/s at the heat floor's peak, so a 1e8 W probe on a
    # stagnant node. The inverse's Newton step from one side could not get
    # below its 1e-6 K tol where T's own resolution is coarser (T > 1.1e9 K),
    # and it threw "no convergence in 50 iterations".
    X = [0.0] * N
    X[cb.species.species_index(a)] = w
    X[cb.species.species_index(b)] = 1.0 - w
    T = 5e11
    assert cb.calc_T_from_h_mass(cb.h_mass(T, X), X, 300.0) == pytest.approx(T, rel=1e-12)
