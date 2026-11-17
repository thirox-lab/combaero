"""A one-stream or one-gas mixer's zero flow partials are exactly 0 (#1198).

``mixer_from_streams_and_jacobians`` forms dY/dm_i, dTt/dm_i and dPt/dm_i
as differences against the mixed state, q_i - q_mix with
q_mix = sum_j m_j q_j / m. Where every flowing stream carries the same q
(one stream, or streams of one gas), the exact partial is 0, but m q / m
rounds, and the partials came out at ~1e-16 relative. They were relayed
into every downstream Jacobian row. The difference is now taken against the
common input value, which cancels exactly. Real partials (a different Tt,
heat, a non-flowing stream of another gas) are unchanged and still match a
central difference.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero import _solver_tools as st

AIR = list(cb.mole_to_mass(cb.species.dry_air()))


def _fuel() -> list[float]:
    Y = [0.0] * len(AIR)
    Y[cb.species.species_index("CH4")] = 1.0
    return Y


def _dY_dm(r: object, i: int) -> list[float]:
    return [row[i].d_mdot for row in r.dY_mix_d_stream]


# Inputs whose m q / m does not round back to q (the issue's reproducer).
SINGLE = [(0.0235, 595.3, 2e5), (0.0208, 700.66, 1.3e5), (3.7, 1234.5, 9.87e5)]


@pytest.mark.parametrize(("m", "Tt", "Pt"), SINGLE)
def test_single_stream_flow_partials_are_exactly_zero(m: float, Tt: float, Pt: float) -> None:
    r = st.mixer_from_streams_and_jacobians([cb.MassStream(m, Tt, Pt, AIR)])
    assert r.dTt_mix_d_stream[0].d_mdot == 0.0
    assert r.dPt_mix_d_stream[0].d_mdot == 0.0
    assert all(v == 0.0 for v in _dY_dm(r, 0))


def test_same_gas_same_state_partials_are_exactly_zero() -> None:
    s = [cb.MassStream(0.0235, 595.3, 2e5, AIR), cb.MassStream(0.0117, 595.3, 2e5, AIR)]
    r = st.mixer_from_streams_and_jacobians(s)
    for i in range(2):
        assert r.dTt_mix_d_stream[i].d_mdot == 0.0
        assert r.dPt_mix_d_stream[i].d_mdot == 0.0
        assert all(v == 0.0 for v in _dY_dm(r, i))


def _fd_dm(streams: list, i: int, Q: float = 0.0) -> tuple[float, float]:
    h = 1e-6 * streams[i].m_dot

    def at(dm: float) -> object:
        s = list(streams)
        o = s[i]
        s[i] = cb.MassStream(o.m_dot + dm, o.Tt, o.Pt, list(o.Y))
        return st.mixer_from_streams_and_jacobians(s, Q=Q)

    rp, rm = at(h), at(-h)
    return (rp.Tt_mix - rm.Tt_mix) / (2 * h), (rp.Pt_mix - rm.Pt_mix) / (2 * h)


def test_same_gas_different_state_keeps_the_real_partials() -> None:
    s = [cb.MassStream(0.0235, 595.3, 2e5, AIR), cb.MassStream(0.0117, 700.0, 2.1e5, AIR)]
    r = st.mixer_from_streams_and_jacobians(s)
    for i in range(2):
        assert all(v == 0.0 for v in _dY_dm(r, i))
        dT, dP = _fd_dm(s, i)
        assert r.dTt_mix_d_stream[i].d_mdot == pytest.approx(dT, rel=1e-6)
        assert r.dPt_mix_d_stream[i].d_mdot == pytest.approx(dP, rel=1e-6)


def test_single_stream_with_heat_keeps_the_heat_partial() -> None:
    s = [cb.MassStream(0.0235, 595.3, 2e5, AIR)]
    r = st.mixer_from_streams_and_jacobians(s, Q=150.0)
    assert all(v == 0.0 for v in _dY_dm(r, 0))
    dT, _ = _fd_dm(s, 0, Q=150.0)
    assert r.dTt_mix_d_stream[0].d_mdot != 0.0
    assert r.dTt_mix_d_stream[0].d_mdot == pytest.approx(dT, rel=1e-6)


def test_a_non_flowing_stream_of_another_gas_keeps_its_partials() -> None:
    # Only flowing streams set the common value: a stream at m = 0 of another
    # gas still moves the mix as soon as it flows.
    s = [cb.MassStream(0.0235, 595.3, 2e5, AIR), cb.MassStream(0.0, 300.0, 2e5, _fuel())]
    r = st.mixer_from_streams_and_jacobians(s)
    assert all(v == 0.0 for v in _dY_dm(r, 0))
    k = cb.species.species_index("CH4")
    assert _dY_dm(r, 1)[k] == pytest.approx(1.0 / 0.0235, rel=1e-12)
    assert r.dTt_mix_d_stream[1].d_mdot != 0.0
