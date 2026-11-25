"""A cooled duct can be bistable: a subsonic root and choked flow (#1307).

Cooling raises a duct's Pt (the Rayleigh gain), so its needed drive
Pt0 - Pt_exit(G) is not monotone in G: it rises from zero flow to a peak and
falls back. At a drive below that peak the duct has a subsonic root on the
rising branch, and choked flow holds too wherever the exit pressure at the
choked flux is at or above the target. As Tt_exit or L moves the peak below
the drive, the subsonic roots merge and vanish (a saddle-node), and the flux
jumps to the choke: no continuous root choice removes that.

The kernel keeps the start-up root (the first stable one, which a slow ramp
from rest reaches). fanno_rayleigh_choked_alternative says whether choked
flow also holds; a duct's diagnostics and the solve's __bistable_ducts__
report it. Case from the issue: Pt0 1e5 Pa, Tt0 400 K, L = D = 1 cm, exit
total, air.
"""

from __future__ import annotations

import math
import warnings
from collections.abc import Callable

import numpy as np
import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    CrossflowSegmentElement,
    FlowNetwork,
    NetworkSolver,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

X = list(cb.species.dry_air())
Y = list(cb.mole_to_mass(cb.species.dry_air()))
PT0, TT0, D = 1e5, 400.0, 0.01


def _flow(Te: float, drive: float, L: float = 0.01) -> cb.FannoChannelFlow:
    return cb.fanno_rayleigh_channel_flow(
        PT0, TT0, Te, X, PT0 - drive, True, L, D, 0.0, "haaland", 1.0, -1.0, [], False
    )


def _alt(Te: float, drive: float, L: float = 0.01) -> cb.FannoChokedAlternative:
    return cb.fanno_rayleigh_choked_alternative(
        PT0, TT0, Te, X, PT0 - drive, True, L, D, 0.0, "haaland"
    )


@pytest.mark.parametrize("drive", [1.5, 3.0, 10.0])
def test_the_jump_to_choke_is_flagged_from_its_subsonic_side(drive: float) -> None:
    """Bisect Tt_exit to the fold, to 1e-11 K: below it the duct chokes, above
    it a subsonic root an order of magnitude smaller in G. That root is
    flagged, and its alternative is the choked flux the other side takes, so
    a solve sitting on either side names the other state."""
    lo, hi = 380.0, 400.0
    assert _flow(lo, drive).choked and not _flow(hi, drive).choked
    for _ in range(41):
        mid = 0.5 * (lo + hi)
        if _flow(mid, drive).choked:
            lo = mid
        else:
            hi = mid
    choked, sub = _flow(lo, drive), _flow(hi, drive)
    alt = _alt(hi, drive)
    # Measured: 7.9x at 1.5 Pa, 5.2x at 3 Pa, 2.5x at 10 Pa.
    assert choked.G / sub.G > 2.0
    assert alt.admissible and alt.exit_margin > 0.0
    assert pytest.approx(choked.G, rel=1e-12) == alt.G


@pytest.mark.parametrize(
    ("label", "values", "call"),
    [
        ("Tt_exit", np.linspace(388.0, 394.0, 25), lambda v: (v, 1.5, 0.01)),
        ("L", np.linspace(0.112, 0.128, 17), lambda v: (250.0, 1.5, float(v))),
    ],
)
def test_every_branch_flip_in_a_scan_is_flagged(
    label: str, values: np.ndarray, call: Callable[[float], tuple[float, float, float]]
) -> None:
    """Across the issue's Tt_exit and L folds (above the 1 Pa drive floor),
    each step between a choked and an unchoked result has the unchoked one
    flagged bistable. The scan must cross a fold, or it tests nothing."""
    flows = [_flow(*call(v)) for v in values]
    flips = [i for i in range(len(flows) - 1) if flows[i].choked != flows[i + 1].choked]
    assert flips, f"the {label} scan crosses no fold"
    for i in flips:
        j = i if not flows[i].choked else i + 1
        assert _alt(*call(values[j])).admissible, (label, values[j])


@pytest.mark.parametrize("Te", [400.0, 420.0, 396.0])
def test_choked_flow_is_not_admissible_where_the_drive_is_monotone(Te: float) -> None:
    """Adiabatic, heated, and cooled too little for the choke's exit Pt to
    reach the target: the exit pressure at choke is below it, no flag."""
    alt = _alt(Te, 1.0)
    assert not alt.admissible
    assert alt.exit_margin < 0.0


def _pair(T_cold: float, reversed_duct: bool = False, segment: bool = False) -> FlowNetwork:
    """A duct from 400 K air at 1 bar to an exit-total boundary 2 Pa lower,
    cooled through a wall by a cold duct beside it. The cold side's T sets
    the hot duct's exit Tt: 100 K puts it inside the bistable band (391 K),
    300 K above it (398 K). reversed_duct declares the hot duct from its exit
    to its inlet (the flow then runs against it); segment makes it a
    CrossflowSegmentElement."""
    g = FlowNetwork()
    g.add_node(PressureBoundary("in", Pt=PT0, Tt=TT0, Y=Y))
    g.add_node(PressureBoundary("out", Pt=PT0 - 2.0, Tt=TT0, Y=Y, coupling="total"))
    g.add_node(PressureBoundary("k_in", Pt=1.2e5, Tt=T_cold, Y=Y))
    g.add_node(PressureBoundary("k_out", Pt=PT0, Tt=T_cold, Y=Y))
    ends = ("out", "in") if reversed_duct else ("in", "out")
    surface_area = 2.0 * math.pi * D * 0.01
    if segment:
        c = CrossflowSegmentElement(
            "c", *ends, length=0.01, area=math.pi * D * D / 4.0, Dh=D, regime="compressible"
        )
        c.surface = ConvectiveSurface(area=surface_area)
        g.add_element(c)
    else:
        g.add_element(
            ChannelElement(
                "c",
                *ends,
                length=0.01,
                diameter=D,
                roughness=0.0,
                regime="compressible",
                surface=ConvectiveSurface(area=surface_area),
            )
        )
    g.add_element(
        ChannelElement(
            "k",
            "k_in",
            "k_out",
            length=0.01,
            diameter=D,
            roughness=0.0,
            regime="compressible",
            surface=ConvectiveSurface(area=surface_area),
        )
    )
    g.add_wall(ThermalWall(id="w", element_a="c", element_b="k", layers=[WallLayer(5e-4, 200.0)]))
    return g


def _solve(g: FlowNetwork) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = NetworkSolver(g).solve(timeout=60)
    assert r["__success__"], r["__message__"]
    return r


def test_a_solved_network_names_a_bistable_duct() -> None:
    """The solve stays the start-up root, and is a success; the duct's keys
    and __bistable_ducts__ name the choked flow it could also carry, at the
    duct's solved end states (the kernel's own verdict on them)."""
    r = _solve(_pair(100.0))
    d = r["__element_diag__"]["c"]
    assert d["choked"] == 0.0 and d["rayleigh_in_momentum"] == 1.0
    assert d["bistable"] == 1.0
    alt = cb.fanno_rayleigh_choked_alternative(
        r["in.Pt"], r["in.Tt"], d["Tt_out"], X, r["out.Pt"], True, 0.01, D, 0.0, "haaland"
    )
    assert alt.admissible
    A = math.pi * D * D / 4.0
    assert d["m_dot_choked_alternative"] == pytest.approx(A * alt.G, rel=1e-12)
    assert d["m_dot_choked_alternative"] > 5.0 * d["m_dot"]
    assert r["__bistable_ducts__"] == {
        "c": {"m_dot": d["m_dot"], "m_dot_choked_alternative": d["m_dot_choked_alternative"]}
    }


def test_a_duct_outside_the_band_reports_the_keys_unset() -> None:
    r = _solve(_pair(300.0))
    for eid in ("c", "k"):
        d = r["__element_diag__"][eid]
        assert d["bistable"] == 0.0
        assert math.isnan(d["m_dot_choked_alternative"])
    assert r["__bistable_ducts__"] == {}


def test_a_static_exit_duct_can_be_bistable_too() -> None:
    """A static exit loses the exit head, which narrows the band but does not
    close it: strongly cooled, choked flow holds beside a subsonic root an
    order of magnitude smaller. Its alternative is the flux the same duct
    chokes at under a large drive."""
    args = (2.64e5, 899.6, 302.0, X)
    geom = (False, 0.738, 0.0243, 0.0, "haaland")
    sub = cb.fanno_rayleigh_channel_flow(*args, 2.64e5 - 19.6, *geom, 1.0, -1.0, [], False)
    alt = cb.fanno_rayleigh_choked_alternative(*args, 2.64e5 - 19.6, *geom)
    hard = cb.fanno_rayleigh_channel_flow(*args, 1e5, *geom, 1.0, -1.0, [], False)
    assert not sub.choked and hard.choked
    assert alt.admissible and alt.G > 10.0 * sub.G
    assert pytest.approx(hard.G, rel=1e-10) == alt.G


@pytest.mark.parametrize(("reversed_duct", "segment"), [(True, False), (False, True)])
def test_the_flag_follows_the_flow_and_the_element(reversed_duct: bool, segment: bool) -> None:
    """Declared against the flow, the duct reports its alternative with the
    flow's sign; as a crossflow segment (its own duct ends, #802), the same
    flag as the channel between the same boundaries. rel 1e-3, not exact: a
    wall's arrangement follows the declared directions (#627), so declaring
    the duct backwards moves its heat by ~3e-5 and its m_dot by ~2e-5."""
    ref = _solve(_pair(100.0))["__element_diag__"]["c"]
    r = _solve(_pair(100.0, reversed_duct=reversed_duct, segment=segment))
    d = r["__element_diag__"]["c"]
    sign = -1.0 if reversed_duct else 1.0
    assert d["bistable"] == 1.0
    assert d["m_dot"] == pytest.approx(sign * ref["m_dot"], rel=1e-3)
    assert d["m_dot_choked_alternative"] == pytest.approx(
        sign * ref["m_dot_choked_alternative"], rel=1e-3
    )
    assert set(r["__bistable_ducts__"]) == {"c"}


def test_a_failing_choke_solve_leaves_the_flag_unknown(monkeypatch: pytest.MonkeyPatch) -> None:
    """The choke solve can fail where a subsonic solve converged (long laminar
    ducts, #1337). A report never fails a solve: both keys are NaN and the
    duct is not listed."""

    def fail(*_args: object) -> None:
        raise RuntimeError("fanno (Mach) brent: no sign change")

    monkeypatch.setattr(cb, "fanno_rayleigh_choked_alternative", fail)
    r = _solve(_pair(100.0))
    d = r["__element_diag__"]["c"]
    assert math.isnan(d["bistable"]) and math.isnan(d["m_dot_choked_alternative"])
    assert r["__bistable_ducts__"] == {}
