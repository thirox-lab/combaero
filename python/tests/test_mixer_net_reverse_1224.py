"""A mixer with net reverse flow is the same mass-weighted mean (#1224).

``mixer_from_streams_and_jacobians`` mixes signed streams into
q_mix = sum_i m_i q_i / m for Y, h and Pt. At a net reverse total
(m < -1e-12) it kept Y averaged but took h and Pt from stream 0, while its
partials were the mean's: they disagreed with a difference (dTt/dm -4085
analytic, 0 by difference), and a non-flowing stream 0 of another gas gave
Tt_mix = -4333 K (its h with the averaged Y). The state is now the mean on
either side of zero, so it and every partial come from one definition.
Reversing every stream gives the same state bit for bit and flips only the
flow partials. Streams of both signs mix only the arriving ones (#1232,
test_mixer_arriving_1232.py), and the near-zero band is gone.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero import _solver_tools as st

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
N = len(AIR)


def _Y(fractions: dict[str, float]) -> list[float]:
    Y = [0.0] * N
    for name, v in fractions.items():
        Y[cb.species.species_index(name)] = v
    return Y


def _mix(streams: list, Q: float = 0.0, fraction: float = 0.0) -> object:
    return st.mixer_from_streams_and_jacobians(streams, Q=Q, fraction=fraction)


def _out(r: object) -> list[float]:
    return [r.Tt_mix, r.Pt_mix, *r.Y_mix]


def _analytic(r: object, i: int, var: str, j: int | None = None) -> list[float]:
    """d(Tt, Pt, Y...)/d(stream i's var) as the mixer reports it."""

    def pick(sj: object) -> float:
        return sj.d_Y[j] if var == "Y" else getattr(sj, f"d_{var}")

    return [
        pick(r.dTt_mix_d_stream[i]),
        pick(r.dPt_mix_d_stream[i]),
        *(pick(row[i]) for row in r.dY_mix_d_stream),
    ]


def _perturbed(streams: list, i: int, var: str, d: float, j: int | None = None) -> list:
    s = list(streams)
    o = s[i]
    m, Tt, Pt, Y = o.m_dot, o.Tt, o.Pt, list(o.Y)
    if var == "mdot":
        m += d
    elif var == "Tt":
        Tt += d
    elif var == "Pt":
        Pt += d
    else:
        Y[j] += d
    s[i] = cb.MassStream(m, Tt, Pt, Y)
    return s


def _step(streams: list, i: int, var: str) -> float:
    o = streams[i]
    return {
        "mdot": 1e-6 * abs(o.m_dot) if o.m_dot != 0.0 else 1e-9,
        "Tt": 1e-3,
        "Pt": 1.0,
        "Y": 1e-7,
    }[var]


def _check_partials(
    streams: list,
    Q: float = 0.0,
    fraction: float = 0.0,
    rel: float = 1e-5,
    one_sided: dict[int, float] | None = None,
) -> None:
    """Every partial of stream i against a central difference, or a
    second-order one-sided one (sign ``one_sided[i]``) where a central step
    would cross into a negative fraction."""
    r = _mix(streams, Q, fraction)
    scale = [max(abs(v), 1.0) for v in _out(r)]
    scale[0] = max(abs(r.Tt_mix), 1.0)
    for i in range(len(streams)):
        for var in ("mdot", "Tt", "Pt", "Y"):
            for j in range(N) if var == "Y" else [None]:
                if var == "Y" and streams[i].m_dot == 0.0:
                    continue
                h = _step(streams, i, var)
                # A mixed fraction is clamped below 0 (normalize_fractions):
                # step a stream's zero fraction, or a flow whose weight
                # m_i / m would turn negative, only the way that keeps it >= 0.
                # A stream in the mix has a positive weight (#1232), so its
                # zero fraction steps up.
                side = None
                if var == "mdot" and one_sided and i in one_sided:
                    side = one_sided[i]
                elif var == "Y" and streams[i].Y[j] == 0.0:
                    side = 1.0
                if side is not None:
                    h *= side
                    p1 = _out(_mix(_perturbed(streams, i, var, h, j), Q, fraction))
                    p2 = _out(_mix(_perturbed(streams, i, var, 2 * h, j), Q, fraction))
                    base = _out(r)
                    fd = [
                        (-3 * b + 4 * a - c) / (2 * h) for b, a, c in zip(base, p1, p2, strict=True)
                    ]
                else:
                    p = _out(_mix(_perturbed(streams, i, var, h, j), Q, fraction))
                    m = _out(_mix(_perturbed(streams, i, var, -h, j), Q, fraction))
                    fd = [(a - b) / (2 * h) for a, b in zip(p, m, strict=True)]
                an = _analytic(r, i, var, j)
                for q, (a, f) in enumerate(zip(an, fd, strict=True)):
                    # Absolute floor: the difference's own rounding, ~ eps
                    # times the output over the step.
                    floor = 1e-12 * scale[q] / abs(h)
                    assert a == pytest.approx(f, rel=rel, abs=floor), (i, var, j, q, a, f)


def test_all_reverse_one_gas_is_the_mean_and_its_partials_match() -> None:
    # Case E of the #1223 review: dTt/dm was -4085 analytic, 0 by difference.
    s = [cb.MassStream(-0.03, 500.0, 2e5, AIR), cb.MassStream(-0.02, 700.0, 2.1e5, AIR)]
    r = _mix(s)
    assert 500.0 < r.Tt_mix < 700.0
    assert r.Pt_mix == pytest.approx((0.03 * 2e5 + 0.02 * 2.1e5) / 0.05, rel=1e-14)
    assert r.dTt_mix_d_stream[1].d_mdot != 0.0
    _check_partials(s)


def test_a_non_flowing_stream_zero_of_another_gas_no_longer_sets_h() -> None:
    # Case F: CH4 at m = 0 as stream 0, two reverse air streams of one state.
    # Tt_mix was -4333 K (CH4's h with air's Y); now it is the air's.
    fuel = _Y({"CH4": 1.0})
    s = [
        cb.MassStream(0.0, 300.0, 3e5, fuel),
        cb.MassStream(-0.02, 700.0, 2.1e5, AIR),
        cb.MassStream(-0.01, 700.0, 2.1e5, AIR),
    ]
    r = _mix(s)
    assert r.Tt_mix == pytest.approx(700.0, rel=1e-12)
    assert r.Pt_mix == 2.1e5
    assert r.Y_mix == pytest.approx(AIR, abs=1e-15)
    # One flowing gas and state: its streams' flow partials are exactly 0.
    for i in (1, 2):
        assert r.dTt_mix_d_stream[i].d_mdot == 0.0
        assert r.dPt_mix_d_stream[i].d_mdot == 0.0
        assert all(row[i].d_mdot == 0.0 for row in r.dY_mix_d_stream)
    # Stream 0 enters with a positive weight m_0 / m only for m_0 < 0.
    _check_partials(s, one_sided={0: -1.0})


@pytest.mark.parametrize(("Q", "fraction"), [(0.0, 0.0), (300.0, 0.2)])
def test_reversing_every_stream_gives_the_same_state_bit_for_bit(Q: float, fraction: float) -> None:
    fwd = [
        cb.MassStream(0.03, 450.0, 2.0e5, _Y({"N2": 0.7, "O2": 0.2, "CO2": 0.1})),
        cb.MassStream(0.02, 800.0, 2.3e5, _Y({"N2": 0.7, "O2": 0.1, "H2O": 0.2})),
        cb.MassStream(0.01, 600.0, 2.1e5, AIR),
    ]
    rev = [cb.MassStream(-s.m_dot, s.Tt, s.Pt, list(s.Y)) for s in fwd]
    a, b = _mix(fwd, Q, fraction), _mix(rev, Q, fraction)
    assert _out(a) == _out(b)
    assert a.dTt_mix_dQ == b.dTt_mix_dQ
    for i in range(len(fwd)):
        for var, sign in (("mdot", -1.0), ("Tt", 1.0), ("Pt", 1.0)):
            assert _analytic(b, i, var) == [sign * v for v in _analytic(a, i, var)]
        for j in range(N):
            assert _analytic(b, i, "Y", j) == _analytic(a, i, "Y", j)


@pytest.mark.parametrize("scale", [1.0, 1e-3, 1e-7, 1e-10])
def test_all_reverse_state_holds_as_the_flow_vanishes(scale: float) -> None:
    # Streams of one sign: the mean does not depend on their common scale, so
    # it holds all the way down to the near-zero band, as forward flow does.
    def streams(k: float) -> list:
        return [
            cb.MassStream(-0.03 * k, 500.0, 2e5, AIR),
            cb.MassStream(-0.02 * k, 900.0, 2.4e5, _Y({"N2": 0.75, "O2": 0.15, "H2O": 0.1})),
        ]

    ref, r = _mix(streams(1.0)), _mix(streams(scale))
    assert _out(r) == pytest.approx(_out(ref), rel=1e-12, abs=1e-15)
    _check_partials(streams(scale))


def test_the_near_zero_band_is_the_mean() -> None:
    # |m| <= 1e-12 took stream 0's state until #1232; now it is the mean of
    # the one-sign streams, as anywhere else, with the band's zero partials.
    s = [cb.MassStream(-4e-13, 500.0, 2e5, AIR), cb.MassStream(-3e-13, 900.0, 2.4e5, AIR)]
    r = _mix(s)
    assert r.Pt_mix == pytest.approx((4 * 2e5 + 3 * 2.4e5) / 7, rel=1e-14)
    assert 500.0 < r.Tt_mix < 900.0
    for i in range(2):
        for var in ("mdot", "Tt", "Pt"):
            assert all(v == 0.0 for v in _analytic(r, i, var))


def test_reverse_one_gas_one_state_partials_are_exactly_zero() -> None:
    s = [cb.MassStream(-0.0235, 595.3, 2e5, AIR), cb.MassStream(-0.0117, 595.3, 2e5, AIR)]
    r = _mix(s)
    for i in range(2):
        assert r.dTt_mix_d_stream[i].d_mdot == 0.0
        assert r.dPt_mix_d_stream[i].d_mdot == 0.0
        assert all(row[i].d_mdot == 0.0 for row in r.dY_mix_d_stream)


def test_reverse_one_gas_different_state_keeps_dY_exactly_zero() -> None:
    s = [cb.MassStream(-0.0235, 595.3, 2e5, AIR), cb.MassStream(-0.0117, 700.0, 2.1e5, AIR)]
    r = _mix(s)
    for i in range(2):
        assert all(row[i].d_mdot == 0.0 for row in r.dY_mix_d_stream)
        assert r.dTt_mix_d_stream[i].d_mdot != 0.0
        assert r.dPt_mix_d_stream[i].d_mdot != 0.0
    _check_partials(s)


# A legacy merging tee with a momentum chamber at its common node: while a
# leg draws flow out of the chamber the tee handed it a net reverse total. The
# chamber's closure reads that mixed state, and the stream-0 h and Pt with
# averaged partials stalled the solve far from the root (case 65: m_com -3.9
# kg/s, norm 889; case 176: m_com -6.0, norm 564). Found by the #1230 review's
# random legacy-tee sweep (5 of 200 cases went no_progress -> converged, none
# worse). #1230 converged both to the forward root.
#
# Since #1231 the tee feeds its port nodes by the sign of its legs, so the
# reversed region carries real states instead of a signed mean, and from the
# same start the root finder lands on a root with the common and straight legs
# reversed: the outlet boundary feeding the higher-pressure straight inlet.
# That root is outside the tee's model (Bassett's merging coefficients at a
# negative flow ratio) and is refused, and the solve retries from a forward
# start (#1234), which converges to the forward root #1230 found.
LEGACY_TEE_CASES = [
    # (P_straight, P_branch, P_out, Tt_straight, Tt_branch, Tt_out, F_C, psi)
    (150000.0, 148760.15882950363, 148260.333624841, 500.0, 900.0, 700.0, 2e-3, 4.0),
    (600000.0, 598294.6280038909, 597033.104368154, 500.0, 900.0, 300.0, 2e-3, 4.0),
]


@pytest.mark.parametrize("case", LEGACY_TEE_CASES, ids=["review_case_65", "review_case_176"])
def test_a_legacy_merging_tee_into_a_momentum_chamber_converges(case: tuple) -> None:
    import math
    import warnings

    from combaero.network import FlowNetwork, NetworkSolver
    from combaero.network.components import (
        LosslessConnectionElement,
        MomentumChamberNode,
        PressureBoundary,
        TeeJunctionElement,
    )

    Pis, Pib, Po, Tis, Tib, To, F_C, psi = case
    g = FlowNetwork()
    g.add_node(PressureBoundary("is", Pt=Pis, Tt=Tis, Y=AIR, coupling="total"))
    g.add_node(PressureBoundary("ib", Pt=Pib, Tt=Tib, Y=AIR, coupling="total"))
    g.add_node(PressureBoundary("o", Pt=Po, Tt=To, Y=AIR, coupling="total"))
    g.add_node(MomentumChamberNode("ps", area=F_C))
    g.add_node(MomentumChamberNode("pb", area=F_C / psi))
    g.add_node(MomentumChamberNode("pc", area=F_C))
    g.add_element(LosslessConnectionElement("ls", "is", "ps"))
    g.add_element(LosslessConnectionElement("lb", "ib", "pb"))
    g.add_element(LosslessConnectionElement("lc", "pc", "o"))
    g.add_element(
        TeeJunctionElement(
            id="tee",
            common_node="pc",
            straight_node="ps",
            branch_node="pb",
            theta=math.radians(90.0),
            F_C=F_C,
            psi=psi,
            tee_type="merging",
        )
    )
    s = NetworkSolver(g)
    s.homotopy_retry = True
    s.bc_ramp_retry = True
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve(timeout=20)
    # The forward root (#1230), reached from a forward start after the cold
    # solve's reversed root was refused (#1231, #1234).
    assert r["__outcome__"] == "converged", (r["__outcome__"], r["__message__"])
    m_com, m_b = r["tee.m_dot_com"], r["tee.m_dot_branch"]
    assert m_com > 0.0 and 0.0 < m_b < m_com
    assert r["lc.m_dot"] == pytest.approx(m_com, rel=1e-9)
    assert min(Tis, Tib) < r["pc.Tt"] < max(Tis, Tib)
    assert "forward start" in r["__message__"]
    v = r["__second_law__"]["tee"]
    assert v["verdict"] == "pass" and v["sigma_norm"] > 0.0
    L = r["__ledger__"]
    assert abs(L["mass_closure_rel"]) < 1e-12 and abs(L["energy_closure_rel"]) < 1e-12
