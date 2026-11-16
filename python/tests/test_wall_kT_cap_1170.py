"""#1170: a ThermalWall's k(T) fixed-point loop settles well inside its cap,
and reports the cap when it does reach it.

The slow walls are a metal on the hot side over alumina, gas at 2000 K and
more, coolant at 300 K and less: the plain loop contracts by ~0.7 a step
there and stopped at its 30 (radiative: 60) steps, silently, with Q
depending on the conductivity the previous call left (up to 5e-8 relative,
2e-5 with a 100 K coolant). After a few plain steps the loop now takes
Anderson steps, and "{wall}.kT_converged" is 0 when the cap still ends it.
"""

from __future__ import annotations

import math
import warnings

import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    GasRadiation,
    HeatLoss,
    MassFlowBoundary,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
    components,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))
STARTS = (300.0, 1125.0, 2000.0, 120.0)


def _stack(metal: str, t: float) -> list[WallLayer]:
    return [WallLayer(t, 20.0, material=metal), WallLayer(t, 20.0, material="alumina")]


def _set_k(w: ThermalWall, T: float) -> None:
    for L in w.layers:
        L.conductivity = cb.get_material_conductivity(L.material, T)


# The issue's reproducer, and the 100 K coolant: both reached the cap from
# every start before.
SLOW = [
    pytest.param(
        "inconel625", 3e-3, (5000.0, 2000.0, 0.01, 1e5, 250.0, 0.01), None, id="reproducer"
    ),
    pytest.param("ss304", 3e-3, (5000.0, 2000.0, 0.01, 1e5, 100.0, 0.01), None, id="coolant_100K"),
    pytest.param("haynes230", 0.1, (200.0, 2000.0, 0.01, 1e4, 100.0, 0.01), "gas", id="radiative"),
]


def _wall(metal: str, t: float, rad: str | None) -> ThermalWall:
    kw = {"radiation_a": GasRadiation(beam_length=0.05)} if rad == "gas" else {}
    return ThermalWall("w", "a", "b", layers=_stack(metal, t), contact_area=0.01, **kw)


def _call(w: ThermalWall, args: tuple, rad: str | None, C_ab=(math.inf, math.inf)):  # noqa: ANN202
    gas = {"gas_a": (args[1], 2e6, AIR)} if rad == "gas" else {}
    return w.compute_coupling(*args, *C_ab, False, **gas)


def _reference_Q(metal: str, t: float, args: tuple, rad: str | None) -> float:
    """The plain loop run to its own settling, far past any cap."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(components, "_KT_PLAIN_STEPS", 10**9)
        mp.setattr(components, "_KT_MAX_STEPS", 5000)
        mp.setattr(components, "_KT_MAX_STEPS_RADIATIVE", 5000)
        w = _wall(metal, t, rad)
        _set_k(w, 0.5 * (args[1] + args[4]))
        Q = float(_call(w, args, rad).Q)
        assert w._kT_converged
    return Q


def _count_steps(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    steps: list[int] = []
    orig = components._KtFixedPoint.step

    def step(
        self: components._KtFixedPoint, x: list[float], g: list[float]
    ) -> tuple[list[float], bool]:
        out = orig(self, x, g)
        if out[1]:
            steps.append(self.steps)
        return out

    monkeypatch.setattr(components._KtFixedPoint, "step", step)
    return steps


@pytest.mark.parametrize(("metal", "t", "args", "rad"), SLOW)
def test_a_wall_that_reached_the_cap_now_settles(
    metal: str, t: float, args: tuple, rad: str | None, monkeypatch: pytest.MonkeyPatch
) -> None:
    """From four starts: settled in at most half the cap, Q at the reference
    to round-off, and the starts agree to round-off."""
    Q_ref = _reference_Q(metal, t, args, rad)
    steps = _count_steps(monkeypatch)
    Qs = []
    for T0 in STARTS:
        w = _wall(metal, t, rad)
        _set_k(w, T0)
        Qs.append(float(_call(w, args, rad).Q))
        assert w._kT_converged
    assert len(steps) == len(STARTS)
    cap = components._KT_MAX_STEPS_RADIATIVE if rad else components._KT_MAX_STEPS
    assert max(steps) <= cap // 2, steps
    assert max(abs(Q - Q_ref) for Q in Qs) <= 1e-11 * abs(Q_ref)
    assert (max(Qs) - min(Qs)) <= 1e-11 * abs(Q_ref)


@pytest.mark.parametrize(("metal", "t", "args", "rad"), SLOW)
def test_the_plain_loop_reached_the_cap_on_them(
    metal: str, t: float, args: tuple, rad: str | None, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The cases are the slow ones: without the acceleration (main's loop)
    every start ends at the cap, now reported."""
    monkeypatch.setattr(components, "_KT_PLAIN_STEPS", 10**9)
    for T0 in (120.0, 2000.0):
        w = _wall(metal, t, rad)
        _set_k(w, T0)
        _call(w, args, rad)
        assert not w._kT_converged


def test_a_call_that_settles_within_the_plain_steps_is_unchanged(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The acceleration starts after _KT_PLAIN_STEPS: a warm or easy call is
    the plain loop's bit for bit, Jacobian included."""
    steps = _count_steps(monkeypatch)
    args = (800.0, 1300.0, 0.05, 150.0, 400.0, 0.05)
    out = []
    for plain in (components._KT_PLAIN_STEPS, 10**9):
        monkeypatch.setattr(components, "_KT_PLAIN_STEPS", plain)
        w = ThermalWall("w", "a", "b", layers=_stack("ss310", 2e-3), contact_area=0.05)
        _set_k(w, 300.0)
        r = w.compute_coupling(*args, 50.0, 20.0, False)
        out.append((r.Q, r.T_hot, r.dQ_dh_a, r.dQ_dT_aw_a, r.dQ_dh_b, r.dQ_dT_aw_b))
    assert 1 < steps[0] <= components._KT_PLAIN_STEPS + 1
    assert out[0] == out[1]


def test_the_coupling_derivatives_hold_at_the_accelerated_exit() -> None:
    """The reproducer's wall between two finite streams: every input
    derivative against central differences of the converged loop."""
    w = ThermalWall("w", "a", "b", layers=_stack("inconel625", 3e-3), contact_area=0.01)
    args = [5000.0, 2000.0, 0.01, 1e5, 250.0, 0.01, 50.0, 80.0]
    r = w.compute_coupling(*args)
    assert w._kT_converged
    for i, name in ((0, "dQ_dh_a"), (1, "dQ_dT_aw_a"), (3, "dQ_dh_b"), (4, "dQ_dT_aw_b")):
        # The h steps well above the loop's 1e-12 settling noise (Q moves
        # little with h_b = 1e5), the T steps small for the k(T) curvature.
        e = (1e-4 if name.startswith("dQ_dh") else 1e-5) * abs(args[i])
        hi, lo = list(args), list(args)
        hi[i] += e
        lo[i] -= e
        fd = (w.compute_coupling(*hi).Q - w.compute_coupling(*lo).Q) / (2 * e)
        assert getattr(r, name) == pytest.approx(fd, rel=1e-7), name


@pytest.mark.parametrize("rad", [None, "gas"])
@pytest.mark.parametrize("cap", [1, 9])
def test_a_capped_call_is_reported(
    rad: str | None, cap: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A cap of one step (plain) and one that ends among the Anderson steps:
    the loop ends unsettled and says so."""
    metal, t, args = "inconel625", 3e-3, (5000.0, 2000.0, 0.01, 1e5, 250.0, 0.01)
    monkeypatch.setattr(components, "_KT_MAX_STEPS", cap)
    monkeypatch.setattr(components, "_KT_MAX_STEPS_RADIATIVE", cap)
    w = _wall(metal, t, rad)
    _set_k(w, 2000.0)
    _call(w, args, rad)
    assert not w._kT_converged
    # The same wall with the default cap settles, and says so.
    monkeypatch.undo()
    _call(w, args, rad)
    assert w._kT_converged


@pytest.mark.parametrize("rad", [None, "gas"])
def test_a_capped_exit_uses_the_layers_own_k(
    rad: str | None, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ending among the Anderson steps, the exchange and the profile are the
    layers' last k (G's image), not the extrapolated iterate: the same wall
    with constant-k layers at those k gives the same Q and profile."""
    monkeypatch.setattr(components, "_KT_MAX_STEPS", 9)
    monkeypatch.setattr(components, "_KT_MAX_STEPS_RADIATIVE", 9)
    metal, t, args = "inconel625", 3e-3, (5000.0, 2000.0, 0.01, 1e5, 250.0, 0.01)
    w = _wall(metal, t, rad)
    _set_k(w, 2000.0)
    r = _call(w, args, rad)
    assert not w._kT_converged
    frozen = _wall(metal, t, rad)
    frozen.layers = [WallLayer(L.thickness, L.conductivity) for L in w.layers]
    r_frozen = _call(frozen, args, rad)
    assert r.Q == r_frozen.Q
    assert w._last_profile == frozen._last_profile


def _pair(layers: list[WallLayer], T_hot: float, T_cold: float) -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=T_hot, Y=AIR))
    n.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=T_hot, Y=AIR))
    n.add_node(MassFlowBoundary("c_in", m_dot=0.02, Tt=T_cold, Y=AIR))
    n.add_node(PressureBoundary("c_out", Pt=1.5e5, Tt=T_cold, Y=AIR))
    for eid, a, b, d in (("hot", "h_in", "h_out", 0.04), ("cold", "c_in", "c_out", 0.01)):
        n.add_element(
            ChannelElement(
                eid, a, b, length=0.3, diameter=d, surface=ConvectiveSurface(model=SmoothModel())
            )
        )
    n.add_wall(ThermalWall("w", "hot", "cold", layers=layers))
    return n


def _solve(net: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(net).solve()
    assert r["__success__"], r["__message__"]
    return r


def test_the_solution_reports_kT_converged(monkeypatch: pytest.MonkeyPatch) -> None:
    """1 for a settled k(T) wall and for a constant-k one; 0 when the report's
    own call ends at the cap (forced: a tolerance no step meets)."""
    r = _solve(_pair(_stack("inconel625", 3e-3), 2000.0, 250.0))
    assert r["w.kT_converged"] == 1.0
    r = _solve(_pair([WallLayer(0.003, 20.0)], 2000.0, 250.0))
    assert r["w.kT_converged"] == 1.0
    monkeypatch.setattr(components, "_KT_TOL", 0.0)
    r = _solve(_pair(_stack("inconel625", 3e-3), 2000.0, 250.0))
    assert r["w.kT_converged"] == 0.0
    r = _solve(_pair([WallLayer(0.003, 20.0)], 2000.0, 250.0))
    assert r["w.kT_converged"] == 1.0


def _one_sided(side_b: object) -> FlowNetwork:
    n = FlowNetwork()
    n.add_node(PressureBoundary("h_in", Pt=1.2e5, Tt=2000.0, Y=AIR))
    n.add_node(PressureBoundary("h_out", Pt=1.1e5, Tt=2000.0, Y=AIR))
    n.add_element(
        ChannelElement(
            "hot",
            "h_in",
            "h_out",
            length=0.3,
            diameter=0.04,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    n.add_wall(
        ThermalWall(
            "w",
            "hot",
            None,
            side_b=side_b,
            layers=_stack("inconel625", 3e-3),
            contact_area=0.01,
        )
    )
    return n


def test_a_heat_loss_wall_reports_kT_converged(monkeypatch: pytest.MonkeyPatch) -> None:
    """The one-sided report: 1 when settled, 0 under a forced cap."""
    r = _solve(_one_sided(HeatLoss(T_ext=250.0, h=1e4)))
    assert r["w.kT_converged"] == 1.0
    monkeypatch.setattr(components, "_KT_TOL", 0.0)
    r = _solve(_one_sided(HeatLoss(T_ext=250.0, h=1e4)))
    assert r["w.kT_converged"] == 0.0


def test_an_adiabatic_wall_reports_kT_converged(monkeypatch: pytest.MonkeyPatch) -> None:
    """An adiabatic wall has the whole stack at one T: no loop, always 1."""
    monkeypatch.setattr(components, "_KT_TOL", 0.0)
    r = _solve(_one_sided("adiabatic"))
    assert r["w.kT_converged"] == 1.0


# Review R2 (#1170): two thick alumina layers, gas radiation on side A, a
# 21 K coolant, from k(120 K): an Anderson extrapolant goes non-positive, which
# the profile would refuse ("t/k values must be non-negative").
R2_LAYERS = [
    WallLayer(0.010649317909138194, 10.0, "alumina"),
    WallLayer(0.017786748253306506, 10.0, "alumina"),
]
R2_ARGS = (619731.6841837598, 1719.8458929850024, 0.01, 965613.0711133536, 21.45645833596489, 0.01)


def test_a_non_positive_extrapolant_falls_back_to_the_plain_step(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The guard fires (Anderson step refused, g taken) and the call settles."""
    rejected: list[int] = []
    orig = components._KtFixedPoint.step

    def step(self: components._KtFixedPoint, x: list[float], g: list[float]) -> tuple:
        out = orig(self, x, g)
        if out[1] or self.steps >= components._KT_MAX_STEPS_RADIATIVE:
            rejected.append(self.rejected)
        return out

    monkeypatch.setattr(components._KtFixedPoint, "step", step)
    w = ThermalWall(
        "w",
        "a",
        "b",
        layers=[WallLayer(L.thickness, L.conductivity, L.material) for L in R2_LAYERS],
        contact_area=0.01,
        radiation_a=GasRadiation(beam_length=0.05),
    )
    _set_k(w, 120.0)
    r = w.compute_coupling(*R2_ARGS, math.inf, math.inf, False, gas_a=(R2_ARGS[1], 2e6, AIR))
    assert w._kT_converged
    assert rejected and rejected[-1] >= 1
    assert math.isfinite(r.Q) and r.Q > 0.0
