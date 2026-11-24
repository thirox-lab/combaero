"""An unknown friction_model refuses instead of making a frictionless duct (#1213).

The C++ dispatcher ``solver::friction_and_jacobian`` returned ``{0, 0}`` with an
INVALID flag for an unknown tag. Every consumer read ``.result`` and dropped the
flag, so ``ChannelElement(friction_model="Haaland")`` solved with f = 0 and
twice the flow of the registered "haaland". Only the zero-flow branch threw,
through ``friction_f_re_at_zero``.
"""

from __future__ import annotations

import pytest

import combaero as cb
from combaero.network import FlowNetwork, NetworkSolver
from combaero.network import components as cbn

BAD = ("Haaland", "nope", "")
KNOWN = sorted(cb._core.Registry.get().available_friction_models())


def test_registry_lists_the_four_models() -> None:
    assert KNOWN == ["colebrook", "haaland", "petukhov", "serghides"]


@pytest.mark.parametrize("tag", BAD)
def test_dispatcher_refuses_an_unknown_tag(tag: str) -> None:
    with pytest.raises(ValueError, match="Unknown friction model"):
        cb._core.friction_and_jacobian(tag, 1.0e4, 0.0)


@pytest.mark.parametrize("tag", KNOWN)
def test_dispatcher_keeps_the_known_tags(tag: str) -> None:
    f, df = cb._core.friction_and_jacobian(tag, 1.0e4, 1.0e-4).result
    assert 0.02 < f < 0.05
    assert df < 0.0


def _channel_rows(tag: str, m_dot: float) -> object:
    # A 1 m, 1 cm air duct; the arguments of the incompressible channel row.
    return cb._core.channel_residuals_and_jacobian(
        m_dot,
        1.01e5,
        1.01e5,
        300.0,
        list(cb.mole_to_mass(cb.species.dry_air())),
        1.0e5,
        1.0,
        0.01,
        0.0,
        tag,
    )


@pytest.mark.parametrize("m_dot", [0.0, 1.0e-6, 1.0e-3])
def test_channel_kernel_refuses_at_any_flow(m_dot: float) -> None:
    # Before #1213 only m = 0 raised (friction_f_re_at_zero); m != 0 gave dP 0.
    with pytest.raises(ValueError, match="Unknown friction model"):
        _channel_rows("nope", m_dot)


@pytest.mark.parametrize(
    "make",
    [
        lambda tag: cbn.ChannelElement(
            "d", "a", "b", length=1.0, diameter=0.01, friction_model=tag
        ),
        lambda tag: cbn.CrossflowSegmentElement(
            "d", "a", "b", length=1.0, area=1.0e-4, friction_model=tag
        ),
    ],
    ids=["ChannelElement", "CrossflowSegmentElement"],
)
@pytest.mark.parametrize("tag", BAD)
def test_elements_refuse_at_construction(make: object, tag: str) -> None:
    with pytest.raises(ValueError, match=r"'d': unknown friction_model"):
        make(tag)  # type: ignore[operator]


@pytest.mark.parametrize("regime", ["incompressible", "compressible"])
def test_the_typo_no_longer_solves(regime: str) -> None:
    # The issue's lone duct: "haaland" solves, "Haaland" now refuses.
    def build(tag: str) -> NetworkSolver:
        g = FlowNetwork()
        g.add_node(cbn.PressureBoundary("a", Pt=1.01e5, Tt=300.0))
        g.add_node(cbn.PressureBoundary("b", Pt=1.0e5, Tt=300.0))
        g.add_element(
            cbn.ChannelElement(
                "d", "a", "b", length=1.0, diameter=0.01, friction_model=tag, regime=regime
            )
        )
        return NetworkSolver(g)

    res = build("haaland").solve()
    assert res["__success__"]
    with pytest.raises(ValueError, match="unknown friction_model 'Haaland'"):
        build("Haaland")
