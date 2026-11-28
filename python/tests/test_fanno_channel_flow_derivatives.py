"""fanno_channel_flow derivatives: low drive, the drive floor, composition.

At low drive the exit and inlet Mach nearly coincide, (Me - M_in)/Me ~
drive/P, and the implicit derivatives combine terms ~P/drive times larger than
their result. The slopes of the isentropic inlet flux and of the exit pressure
per unit flux were differenced, and that rounding, amplified ~1e5, put up to
8e-3 into dG/dPt0 and dG/dP at a few Pa of drive and 1e-5 noise into G below
the 1 Pa floor (through the blend exponent). Both are exact now.

dG_ddir is dG/ds along mass-fraction directions v (Y -> Y + s v). A unit
vector e_k gives dG/dY_k with the other fractions fixed; a network element
passes the one or two directions its feeding node's composition moves in.
"""

from __future__ import annotations

import numpy as np
import pytest

import combaero as cb

N = cb.num_species()
PT, TT, L, D, ROUGH = 2.0e5, 700.0, 1.0, 0.02, 1e-5


def _X(d: dict[str, float]) -> list[float]:
    X = [0.0] * N
    for k, v in d.items():
        X[cb.species_index_from_name(k)] = v
    s = sum(X)
    return [x / s for x in X]


MIXES = {
    "air": list(cb.species.dry_air()),
    "products": _X({"N2": 0.71, "H2O": 0.18, "CO2": 0.09, "O2": 0.02}),
    "co2_rich": _X({"N2": 0.5, "O2": 0.13, "CO2": 0.3, "H2O": 0.07}),
}


def _flow(
    X: list[float],
    P_target: float,
    exit_total: bool = False,
    Pt0: float = PT,
    Tt0: float = TT,
    dirs: list[list[float]] | None = None,
    derivatives: bool = True,
) -> cb.FannoChannelFlow:
    return cb.fanno_channel_flow(
        Pt0,
        Tt0,
        X,
        P_target,
        exit_total,
        L,
        D,
        ROUGH,
        "haaland",
        1.0,
        -1.0,
        dirs or [],
        derivatives,
    )


def _unit(k: int) -> list[float]:
    return [float(j == k) for j in range(N)]


@pytest.mark.parametrize("drive", [1.5, 3.0, 10.0, 100.0, 1.0e4])
@pytest.mark.parametrize("exit_total", [False, True])
def test_pressure_and_temperature_derivatives_hold_at_low_drive(
    drive: float, exit_total: bool
) -> None:
    """<= 1e-5 down to 1.5 Pa on 2 bar (was up to 8e-3)."""
    X = MIXES["products"]
    P = PT - drive
    f = _flow(X, P, exit_total)
    hp = 2e-3 * drive
    ht = 1e-4 * TT
    fd_pt = (_flow(X, P, exit_total, Pt0=PT + hp).G - _flow(X, P, exit_total, Pt0=PT - hp).G) / (
        2 * hp
    )
    fd_p = (_flow(X, P + hp, exit_total).G - _flow(X, P - hp, exit_total).G) / (2 * hp)
    fd_t = (_flow(X, P, exit_total, Tt0=TT + ht).G - _flow(X, P, exit_total, Tt0=TT - ht).G) / (
        2 * ht
    )
    assert f.dG_dPt0 == pytest.approx(fd_pt, rel=1e-5)
    assert f.dG_dP_target == pytest.approx(fd_p, rel=1e-5)
    assert f.dG_dTt0 == pytest.approx(fd_t, rel=1e-5)


@pytest.mark.parametrize("drive", [0.4, 0.8])
def test_the_flux_is_smooth_below_the_drive_floor(drive: float) -> None:
    """Second differences over 1 mK steps in Tt: rounding level (~1e-9 of G).
    The blend exponent sigma carried the differenced slopes' noise, 1e-5 of G."""
    X = MIXES["air"]
    T = TT + np.linspace(0.0, 0.05, 51)
    G = np.array([_flow(X, PT - drive, Tt0=t).G for t in T])
    assert np.max(np.abs(np.diff(G, 2))) / G.mean() < 1e-7


def _dG_dv_fd(X: list[float], P: float, exit_total: bool, v: list[float], h0: float) -> float:
    """Richardson-extrapolated central difference along v (Y -> Y + s v)."""
    Y = list(cb.mole_to_mass(X))

    def at(s_: float) -> float:
        return _flow(
            list(cb.mass_to_mole([y + s_ * vk for y, vk in zip(Y, v, strict=True)])), P, exit_total
        ).G

    def D(h: float) -> float:
        return (at(h) - at(-h)) / (2 * h)

    return (4 * D(h0 / 2) - D(h0)) / 3


def _dG_dY_fd(X: list[float], P: float, exit_total: bool, k: int) -> float:
    """Per species, others fixed."""
    Y = cb.mole_to_mass(X)
    return _dG_dv_fd(X, P, exit_total, _unit(k), min(2e-4, 0.4 * Y[k]))


@pytest.mark.parametrize("mix", sorted(MIXES))
@pytest.mark.parametrize("P_target", [1.6e5, 1.0e5, 0.5e5], ids=["unchoked", "near", "choked"])
@pytest.mark.parametrize("exit_total", [False, True])
def test_the_composition_derivative_matches_differences(
    mix: str, P_target: float, exit_total: bool
) -> None:
    X = MIXES[mix]
    Y = cb.mole_to_mass(X)
    present = [k for k in range(N) if Y[k] > 0.0]
    f = _flow(X, P_target, exit_total, dirs=[_unit(k) for k in present])
    assert len(f.dG_ddir) == len(present)
    for d, k in zip(f.dG_ddir, present, strict=True):
        fd = _dG_dY_fd(X, P_target, exit_total, k)
        assert d == pytest.approx(fd, rel=1e-5, abs=1e-6 * f.G)


@pytest.mark.parametrize("P_target", [1.6e5, 0.5e5], ids=["unchoked", "choked"])
def test_a_mixing_direction_is_the_projected_gradient(P_target: float) -> None:
    """Along v = Y_a - Y_b (two streams mixing) the directional derivative is
    the per-species gradient dotted with v -- one evaluation instead of one
    per species."""
    Ya, Yb = cb.mole_to_mass(MIXES["air"]), cb.mole_to_mass(MIXES["products"])
    v = [a - b for a, b in zip(Ya, Yb, strict=True)]
    Y = [0.6 * a + 0.4 * b for a, b in zip(Ya, Yb, strict=True)]  # the node holds both streams
    X = list(cb.mass_to_mole(Y))
    present = [k for k in range(N) if Y[k] > 0.0]
    d_v = _flow(X, P_target, dirs=[v]).dG_ddir[0]
    grad = _flow(X, P_target, dirs=[_unit(k) for k in present]).dG_ddir
    assert d_v == pytest.approx(sum(g * v[k] for g, k in zip(grad, present, strict=True)), rel=1e-5)
    assert d_v == pytest.approx(_dG_dv_fd(X, P_target, False, v, 1e-3), rel=1e-5)


def test_a_direction_into_an_absent_species_is_one_sided() -> None:
    """H2 absent from the mixture but carried by an inflow: the central step
    would make Y_H2 negative, so the difference is one-sided."""
    X = MIXES["air"]
    Y = cb.mole_to_mass(X)
    k_h2 = cb.species_index_from_name("H2")
    assert Y[k_h2] == 0.0
    v = [-y for y in Y]
    v[k_h2] += 1.0  # toward pure H2
    f = _flow(X, 1.6e5, dirs=[v])
    Yv = list(Y)

    def at(s_: float) -> float:
        return _flow(
            list(cb.mass_to_mole([y + s_ * vk for y, vk in zip(Yv, v, strict=True)])), 1.6e5
        ).G

    fd = (4 * (at(5e-5) - at(0.0)) / 5e-5 - (at(1e-4) - at(0.0)) / 1e-4) / 3
    assert f.dG_ddir[0] == pytest.approx(fd, rel=1e-3)


@pytest.mark.parametrize("drive", [3.0, 30.0])
def test_the_composition_derivative_holds_at_low_drive(drive: float) -> None:
    """The cancellation is ~P/drive there: <= 3e-3 of the column at 3 Pa on
    2 bar, <= 1e-4 at 30 Pa, measured against the column's own scale."""
    X = MIXES["products"]
    P = PT - drive
    Y = cb.mole_to_mass(X)
    present = [k for k in range(N) if Y[k] > 0.0]
    f = _flow(X, P, dirs=[_unit(k) for k in present])
    for d, k in zip(f.dG_ddir, present, strict=True):
        fd = _dG_dY_fd(X, P, False, k)
        scale = max(abs(fd), 0.05 * f.G)
        assert abs(d - fd) / scale < (3e-3 if drive < 10 else 1e-4)


def test_the_composition_derivative_is_opt_in() -> None:
    X = MIXES["products"]
    v = [_unit(cb.species_index_from_name("CO2"))]
    assert len(_flow(X, 1.6e5).dG_ddir) == 0
    zero = _flow(X, PT + 10.0, dirs=v)
    assert zero.G == 0.0 and list(zero.dG_ddir) == [0.0]
    floor = _flow(X, PT - 0.5, dirs=v)
    assert len(floor.dG_ddir) == 1 and floor.dG_ddir[0] != 0.0


@pytest.mark.parametrize("P_target", [1.6e5, PT - 0.5, 0.5e5], ids=["unchoked", "floor", "choked"])
def test_without_derivatives_the_flow_is_the_same(P_target: float) -> None:
    """A residual-only evaluation skips the derivatives, not the flow."""
    X = MIXES["products"]
    full, bare = _flow(X, P_target), _flow(X, P_target, derivatives=False)
    assert pytest.approx(full.G, rel=1e-12) == bare.G
    assert bare.choked == full.choked and bare.M_exit == pytest.approx(full.M_exit, rel=1e-12)
    if P_target != PT - 0.5:
        assert bare.dG_dPt0 == 0.0 and bare.dG_dTt0 == 0.0


@pytest.mark.parametrize("length", [0.0127, 1.0], ids=["short", "long"])
@pytest.mark.parametrize("drive", [5.0, 29.0, 400.0, 2.0e4])
def test_the_flow_does_not_depend_on_the_warm_start(length: float, drive: float) -> None:
    """The exit-Mach guess only seeds the bracket. A guess far below the root
    could not be widened past ~3.6x itself, so Brent ran without a sign change
    and returned the bracket's end: a short duct warm-started at M 0.01 gave
    G = 13.5 for every drive from 5 to 400 Pa (#585)."""
    X = MIXES["air"]

    def g(guess: float) -> float:
        return cb.fanno_channel_flow(
            1.0e5,
            300.0,
            X,
            1.0e5 - drive,
            True,
            length,
            0.028,
            0.0,
            "haaland",
            1.0,
            guess,
            [],
            False,
        ).G

    cold = g(-1.0)
    for guess in (1e-4, 0.01, 0.1, 0.5, 0.99):
        assert g(guess) == pytest.approx(cold, rel=1e-9), guess


def _air_without_co2() -> list[float]:
    """Dry air mass fractions with CO2 moved into N2: CO2 absent (Y = 0)."""
    Y = list(cb.mole_to_mass(cb.species.dry_air()))
    i_co2 = cb.species_index_from_name("CO2")
    Y[cb.species_index_from_name("N2")] += Y[i_co2]
    Y[i_co2] = 0.0
    return Y


def _trade(a: str, b: str) -> np.ndarray:
    v = np.zeros(N)
    v[cb.species_index_from_name(a)] = 1.0 / np.sqrt(2.0)
    v[cb.species_index_from_name(b)] = -1.0 / np.sqrt(2.0)
    return v


@pytest.mark.parametrize("P_target", [1.5e5, 0.5e5])
def test_a_roundoff_component_on_an_absent_species_moves_dG_ddir_by_roundoff(
    P_target: float,
) -> None:
    """A round-off component on an absent species moves dG_ddir linearly, by
    round-off. Under d_ddir its sign alone chose a one-sided stencil and moved
    the entry by 2.3e-5 relative in a network (#1199, stress random 30); the
    exact tangent (#1311) is linear in v, so nothing switches."""
    Y = _air_without_co2()
    X = list(cb.mass_to_mole(Y))
    v = _trade("O2", "N2")
    exact = _flow(X, P_target, dirs=[list(v)]).dG_ddir[0]
    assert np.isfinite(exact) and exact != 0.0
    for c in (-2.2e-20, 2.2e-20, -1e-14):
        w = v.copy()
        w[cb.species_index_from_name("CO2")] = c
        assert _flow(X, P_target, dirs=[list(w)]).dG_ddir[0] == pytest.approx(exact, rel=1e-12)


def test_ddir_one_sided_is_second_order() -> None:
    """Along a genuine direction into an absent species only s > 0 is feasible.
    The one-sided stencil is second order: within 3e-7 of central derivatives
    at s = 2h, 4h, 6h extrapolated back to 0 (the first-order one was off by
    3.3e-6)."""
    Y = np.asarray(_air_without_co2())
    v = _trade("CO2", "N2")
    h = 1e-5 * np.sqrt(2.0)

    def at(s: float) -> float:
        X = list(cb.mass_to_mole(list(Y + s * v)))
        return _flow(X, 1.5e5, dirs=[list(v)]).dG_ddir[0]

    ref = 3.0 * at(2 * h) - 3.0 * at(4 * h) + at(6 * h)
    assert at(0.0) == pytest.approx(ref, rel=3e-7)


def test_ddir_halved_step_where_doubled_step_infeasible() -> None:
    """CO2 absent blocks s < 0; H2O at 1.5e-5 (h|v_H2O| = 1e-5) blocks 2h.
    The one-sided stencil runs at half the step, still second order (#1199
    review: the first-order fallback was off by 1.2e-5, a jump in dG/dv as
    Y_H2O crossed 2e-5). Reference: full-step values at Y_H2O = 3, 4, 5e-5
    extrapolated quadratically. -v gives exactly the negated result."""
    v = _trade("CO2", "H2O")

    def d(y: float, w: np.ndarray) -> float:
        Y = np.asarray(_air_without_co2())
        Y[cb.species_index_from_name("H2O")] = y
        Y[cb.species_index_from_name("N2")] -= y
        return _flow(list(cb.mass_to_mole(list(Y))), 1.5e5, dirs=[list(w)]).dG_ddir[0]

    ref = 4.375 * d(3e-5, v) - 5.25 * d(4e-5, v) + 1.875 * d(5e-5, v)
    assert d(1.5e-5, v) == pytest.approx(ref, rel=1e-6)
    assert d(1.5e-5, -v) == -d(1.5e-5, v)
