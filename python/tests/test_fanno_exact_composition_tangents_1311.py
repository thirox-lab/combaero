"""Exact composition tangents of the Fanno kernels (#1311 adiabatic, #1310
Fanno-Rayleigh).

dG/dv along a mass-fraction direction v comes from cancelling sums at low
drive (the inlet and exit Mach numbers differ by ~drive / P). Until #1311 the
pieces were differenced (d_ddir), and the cancellation amplified their ~1e-8
error: on main e427b904 this grid's analytic/difference ratio near the 1 Pa
floor was off by up to 8 % (adiabatic) and 21 % (heated Fanno-Rayleigh, pure
O2 at 1.05 Pa). Forward mode along v carries round-off only.

The reference is G itself, differenced along v by a two-level Richardson
ladder (Y steps 8e-2 to 1e-2). Near the floor G's root carries ~2e-9
relative jitter, and dlnG/dv is only ~1.5e-3 there, so no difference of G
resolves dG/dv much below ~3e-5. The bound there is #1311's own, 1e-4 (or
three times the ladder's measured gap, where that is larger, and the gap
must stay below 5e-4 so a poor reference cannot pass a poor tangent); at
300 Pa, where the reference resolves to ~1e-8, it is 1e-6.
"""

from __future__ import annotations

import pytest

import combaero as cb

NS = cb.num_species()
I_N2, I_O2, I_AR = (cb.species_index_from_name(n) for n in ("N2", "O2", "AR"))
PT0, TT0, TT_EXIT, L, D = 2e5, 400.0, 500.0, 0.5, 0.01


def _air() -> list[float]:
    Y = [0.0] * NS
    Y[I_N2], Y[I_O2], Y[I_AR] = 0.7553, 0.2314, 0.0133
    s = sum(Y)
    return [y / s for y in Y]


Y_AIR = _air()


def _dir(**c: float) -> list[float]:
    v = [0.0] * NS
    for name, x in c.items():
        v[cb.species_index_from_name(name)] = x
    return v


O2_N2 = _dir(O2=0.5, N2=-0.5)
O2 = _dir(O2=1.0)


def _flow(Y: list[float], drive: float, total: bool, heated: bool, dirs: list[list[float]]):
    X = cb.mass_to_mole(Y)
    with_d = bool(dirs)
    if heated:
        return cb.fanno_rayleigh_channel_flow(
            PT0, TT0, TT_EXIT, X, PT0 - drive, total, L, D, 1e-5, "haaland", 1.0, -1.0, dirs, with_d
        )
    return cb.fanno_channel_flow(
        PT0, TT0, X, PT0 - drive, total, L, D, 1e-5, "haaland", 1.0, -1.0, dirs, with_d
    )


def _reference(v: list[float], drive: float, total: bool, heated: bool) -> tuple[float, float]:
    def G(h: float) -> float:
        return _flow(
            [y + h * vk for y, vk in zip(Y_AIR, v, strict=True)], drive, total, heated, []
        ).G

    def central(h: float) -> float:
        return (G(h) - G(-h)) / (2.0 * h)

    d = [central(h) for h in (8e-2, 4e-2, 2e-2, 1e-2)]
    r = [(4.0 * d[i + 1] - d[i]) / 3.0 for i in range(3)]
    q = [(16.0 * r[i + 1] - r[i]) / 15.0 for i in range(2)]
    return q[1], abs(q[1] - q[0]) / abs(q[1])


CASES = [
    # (heated, direction name, drive [Pa]); the 1 Pa floor and just above it,
    # where the cancellation is worst, then a resolved drive.
    *[(False, name, drive) for name in ("O2/N2", "O2") for drive in (1.0, 1.05, 1.5, 300.0)],
    *[(True, "O2/N2", drive) for drive in (1.0, 1.05, 1.5, 300.0)],
    # Pure O2 into the heated duct: near 1 Pa dG/dY_O2 crosses zero (4.6e-7
    # against 4.5e-4 for N2), where no difference resolves it; from 1.5 Pa on.
    *[(True, "O2", drive) for drive in (1.5, 300.0)],
]


@pytest.mark.parametrize("total", [False, True], ids=["static-exit", "total-exit"])
@pytest.mark.parametrize(("heated", "name", "drive"), CASES)
def test_dG_ddir_matches_a_resolved_difference(
    heated: bool, name: str, drive: float, total: bool
) -> None:
    v = O2_N2 if name == "O2/N2" else O2
    got = _flow(Y_AIR, drive, total, heated, [v]).dG_ddir[0]
    ref, gap = _reference(v, drive, total, heated)
    assert gap < 5e-4, f"reference unresolved: {gap:.1e}"
    bound = 1e-6 if drive >= 300.0 else max(1e-4, 3.0 * gap)
    assert got == pytest.approx(ref, rel=bound)


@pytest.mark.parametrize("heated", [False, True], ids=["adiabatic", "heated"])
@pytest.mark.parametrize("drive", [1.0, 1.5, 10.0])
def test_dG_ddir_is_linear_in_the_direction(heated: bool, drive: float) -> None:
    """dG/dv is linear in v: (e_O2 - e_N2) / 2 gives the O2/N2 trade. A
    stencil chosen per direction (d_ddir) broke this; forward mode cannot."""
    r = _flow(Y_AIR, drive, False, heated, [O2_N2, O2, _dir(N2=1.0)])
    a, b, c = r.dG_ddir
    assert a == pytest.approx(0.5 * (b - c), rel=1e-7)
