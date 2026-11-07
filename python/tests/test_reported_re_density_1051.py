"""The reported Reynolds number carries the density (#1051).

The elements' reference block and the chamber/plenum diagnostics reported
v Dh / mu, low by a factor of rho, and the channel's reported f was evaluated
at it. The rows were right: they form Re from the mass flux. At 5 bar of air
(rho ~ 5.8 kg/m^3) the old value is off by that factor, so the checks below
cannot pass by accident."""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    FlowNetwork,
    MomentumChamberNode,
    NetworkSolver,
    OrificeElement,
    PressureBoundary,
)

D = 0.1  # [m] pipe
D_BORE = 0.05  # [m] plate
ROUGHNESS = 1e-5  # [m], the channel default
FRICTION_MODEL = "haaland"  # the channel default


@pytest.fixture(scope="module")
def run() -> dict:
    net = FlowNetwork()
    net.add_node(PressureBoundary("supply", Pt=5.0e5, Tt=300.0))
    net.add_node(MomentumChamberNode("up"))
    net.add_node(MomentumChamberNode("down"))
    net.add_node(PressureBoundary("sink", Pt=4.95e5, Tt=300.0))
    net.add_element(
        ChannelElement("pipe", "supply", "up", length=2.0, diameter=D, roughness=ROUGHNESS)
    )
    net.add_element(OrificeElement("plate", "up", "down", diameter=D_BORE))
    net.add_element(ChannelElement("pipe_down", "down", "sink", length=1.0, diameter=D))
    r = NetworkSolver(net).solve()
    assert r["__success__"], r["__message__"]
    assert r["pipe.rho_ref"] > 5.0
    return r


def test_channel_re_is_the_mass_flux_re(run: dict) -> None:
    expected = 4.0 * abs(run["pipe.m_dot"]) / (math.pi * D * run["pipe.mu_ref"])
    assert run["pipe.Re"] == pytest.approx(expected, rel=1e-9)
    # Its surface correlation forms Re from the mass flux on its own.
    assert run["pipe.Re"] == pytest.approx(run["pipe.Re_surface"], rel=1e-2)


def test_channel_reports_f_at_that_re(run: dict) -> None:
    e_D = ROUGHNESS / D
    f = cb._core.friction_and_jacobian(FRICTION_MODEL, run["pipe.Re"], e_D).result[0]
    assert run["pipe.f"] == pytest.approx(f, rel=1e-9)


def test_orifice_re_is_the_bore_re(run: dict) -> None:
    expected = 4.0 * abs(run["plate.m_dot"]) / (math.pi * D_BORE * run["plate.mu_ref"])
    assert run["plate.Re"] == pytest.approx(expected, rel=1e-9)


def test_chamber_re_is_the_mass_flux_re(run: dict) -> None:
    m = abs(run["pipe.m_dot"])
    expected = m * run["up.Dh"] / (run["up.area"] * run["up.mu"])
    assert run["up.Re"] == pytest.approx(expected, rel=1e-9)
