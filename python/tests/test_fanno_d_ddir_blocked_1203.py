"""dG/dv along a composition direction that no stencil can step (#1203).

`d_ddir` (src/fanno_mach.cpp) differences a duct's flux along a mass-fraction
direction v. Where v has opposite-sign components on two species that are
absent or below h|v_k|, neither side of the step keeps every fraction >= 0,
and it returned 0: the duct's J entry was 0 where the true dG/dv is 58-66.
dG/dv is linear in v, so the direction is now split into its blocking species
(each a unit vector, feasible on its positive side) and the rest, and the
parts are summed. The reference here is that linearity: sum_k v_k dG/dY_k from
per-species unit directions.
"""

import numpy as np
import pytest

import combaero as cb

_N = [cb.species_name(k) for k in range(cb.num_species())]
_I = {n: k for k, n in enumerate(_N)}


def _gas(trace: dict[str, float]) -> list[float]:
    """Dry air with its CO2 moved into N2 (as in #1203), plus mass-fraction traces."""
    Y = np.zeros(len(_N))
    for name, y in {"N2": 0.7557, "O2": 0.2314, "AR": 0.0129, **trace}.items():
        Y[_I[name]] = y
    Y /= Y.sum()
    return list(cb.mass_to_mole(list(Y)))


def _vec(d: dict[str, float]) -> list[float]:
    v = [0.0] * len(_N)
    for name, x in d.items():
        v[_I[name]] = x
    return v


def _dG_ddir(X: list[float], dirs: list[list[float]]) -> list[float]:
    r = cb._core.fanno_channel_flow(
        2e5, 600.0, X, 1.5e5, False, 0.5, 0.02, 0.0, "haaland", 1.0, -1.0, dirs
    )
    return list(r.dG_ddir)


@pytest.mark.parametrize(
    ("trace", "direction"),
    [
        ({}, {"CO2": 0.5, "H2O": -0.5, "O2": 0.5, "N2": -0.5}),
        ({"CO2": 1e-7, "H2O": 1e-7}, {"CO2": 0.5, "H2O": -0.5, "O2": 0.5, "N2": -0.5}),
        ({"CO2": 1e-7, "H2O": 1e-7}, {"CO2": 0.5, "H2O": -0.5}),
        ({"H2O": 0.9e-5}, {"CO2": 0.5, "H2O": -0.5}),
    ],
    ids=["both-absent", "both-trace", "pure-trace-trade", "co2-absent-h2o-trace"],
)
def test_a_direction_no_stencil_can_step_takes_its_linear_value(
    trace: dict[str, float], direction: dict[str, float]
) -> None:
    X = _gas(trace)
    got = _dG_ddir(X, [_vec(direction)])[0]
    parts = _dG_ddir(X, [_vec({name: 1.0}) for name in direction])
    ref = sum(direction[name] * p for name, p in zip(direction, parts, strict=True))
    assert abs(ref) > 10.0  # an order-one entry, not a vanishing one
    assert got == pytest.approx(ref, rel=1e-7)


def test_no_step_where_a_trace_crosses_the_stencil() -> None:
    """dG/dv is continuous as the H2O trace rises past the one-sided step
    (0.9e-5 split, 1.2e-5 one-sided): on main it went 0 -> 58."""
    v = _vec({"CO2": 0.5, "H2O": -0.5})
    below = _dG_ddir(_gas({"H2O": 0.9e-5}), [v])[0]
    above = _dG_ddir(_gas({"H2O": 1.2e-5}), [v])[0]
    assert below == pytest.approx(above, rel=1e-5)
