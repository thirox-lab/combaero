"""An orifice in a pipe at zero flow (#997).

With an approach pipe (0 < beta < 1) the Jobson kernel drives from the pipe's
static pressure p1 = p(Pt, G = |m| / A_pipe). Two roundings met at zero flow:

- duct_static_state solved in ln P, and exp(log(Pt)) came back one ulp ABOVE
  Pt at G = 0: a static above its own total;
- the Jobson momentum balance wrote the jet's 1 - (p / p0)^((n-1)/n) as a
  difference near 1, which rounded to 0 a few ulp below p0, so u = 0 against a
  drive p1 - p > 0 and g = 2 c / (B + root) was inf.

The inf reached the approach static of the next evaluation as G = inf and
chamber_static_state's Newton returned NaN: every Jobson variant in a pipe
failed at equal boundary pressures (HTTP 500 under the strict suite).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

import combaero as cb
from combaero import _solver_tools

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "scripts"))

import gui_result_keys as keys  # noqa: E402
import test_gui_http_routes as http  # noqa: E402

Y_AIR = list(cb.mole_to_mass(cb.species.dry_air()))


@pytest.mark.parametrize("Pt", [1.0e5, 1.05e5, 3.7e5, 12345.678, 2.0e6])
@pytest.mark.parametrize("Tt", [300.0, 900.0])
def test_the_approach_static_at_zero_flux_is_the_total_to_the_bit(Pt: float, Tt: float) -> None:
    a = _solver_tools.approach_static_pressure(Pt, Tt, 0.0, Y_AIR, False)
    assert Pt == a.P
    for G in (1e-12, 1e-6, 1e-3, 1.0, 50.0):
        assert Pt >= _solver_tools.approach_static_pressure(Pt, Tt, G, Y_AIR, False).P


def test_the_jobson_kernel_is_finite_at_a_one_ulp_drive() -> None:
    """Finite down to one ulp of the drive, both directions. Below the 1e-3 Pa
    drive floor (#1256) the flow is the floor's times the odd cubic h(t),
    t = dp / floor, h = (3 - sigma)/2 t + (sigma - 1)/2 t^3, sigma = 1/2
    (kJetFloorExponent): linear in dp at one ulp."""
    Pt, Tt, beta = 1.0e5, 300.0, 0.5
    area = math.pi * 0.05**2 / 4.0
    floor = 1.0e-3  # kJetDriveFloor [Pa], flow_drive_floor.h

    def kernel(p_down: float) -> object:
        return _solver_tools.orifice_jobson_residuals_and_jacobian(
            0.0, Pt, Pt, Tt, Y_AIR, p_down, 0.6176, area, beta, 1.0
        )

    at = kernel(Pt - floor)
    sigma = 0.5
    k_ref = at.m_dot_calc * (3.0 - sigma) / (2.0 * floor)  # m / dp as dp -> 0
    p = Pt
    for _ in range(40):
        p = math.nextafter(p, 0.0)
        r = kernel(p)
        values = (
            r.m_dot_calc,
            r.d_mdot_dP_static_up,
            r.d_mdot_dP_static_down,
            r.d_mdot_dP_total_up,
            r.d_mdot_dTt_up,
            r.d_mdot_dCd,
        )
        assert all(math.isfinite(v) for v in values), (Pt - p, values)
        assert r.m_dot_calc > 0.0
        assert r.m_dot_calc / (Pt - p) == pytest.approx(k_ref, rel=1e-6)
    # The reverse band: the back pressure one ulp above the drive.
    r = kernel(math.nextafter(Pt, math.inf))
    assert math.isfinite(r.m_dot_calc) and r.m_dot_calc < 0.0
    assert math.isfinite(r.d_mdot_dP_static_down)


@pytest.mark.parametrize(
    "variant",
    [
        {},
        {"compressible_model": "jobson1955"},
        {"correlation": "Lichtarowicz", "plate_thickness": 0.1},
        {"correlation": "IdelchikThick", "plate_thickness": 0.1},
        {"correlation": "McGreehanSchotsch", "plate_thickness": 0.1},
    ],
    ids=["default", "jobson1955", "Lichtarowicz", "IdelchikThick", "McGreehanSchotsch"],
)
def test_an_orifice_in_a_pipe_solves_at_zero_flow(variant: dict) -> None:
    graph = keys._orifice_in_pipe(variant)
    for node in graph["nodes"]:
        if node["id"] == "in":
            node["data"]["Pt"] = 1.0e5  # the outlet's: no drive
    r = http._post("/solve", graph)
    assert r.status_code == 200, r.text[:300]
    body = r.json()
    assert body["success"], body["message"]
    m = body["element_results"]["x"]["m_dot"]
    assert abs(m) < 1e-12, m
