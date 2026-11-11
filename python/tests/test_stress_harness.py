"""The stress harness (scripts/stress) keeps building valid networks.

A dev tool, not shipped: it builds every element family cold and over its
range, which is how #614 and #619-#624 were found. This only guards that its
generators still construct and validate when an element's API changes; the
sweeps themselves are run by hand (scripts/stress/run.py).
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "stress"))

from networks import (  # noqa: E402
    EXPECT_REFUSED_AT_BUILD,
    FAMILIES,
    FUEL_CASES,
    FUEL_ROUTES,
)


def _builds(family: str, seed: int) -> list[str]:
    """Build a draw and resolve and validate it as the solver does (theta
    sources resolve here); a draw declared refused before the solve must be,
    for its reason."""
    net, tags = FAMILIES[family](seed)
    assert net.elements and tags
    refused = [
        t[len(EXPECT_REFUSED_AT_BUILD) :] for t in tags if t.startswith(EXPECT_REFUSED_AT_BUILD)
    ]
    if refused:
        with pytest.raises(ValueError, match=refused[0]):
            net.resolve_all_topology()
    else:
        net.resolve_all_topology()
        net.validate()
    return tags


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_every_family_builds_valid_networks(family: str) -> None:
    for seed in range(5):
        _builds(family, seed)


def test_every_fuel_controller_case_and_route_builds() -> None:
    """#1049: each case on each route, in the first len(cases) x len(routes)
    seeds, builds and declares exactly one expected outcome."""
    seen = set()
    for seed in range(len(FUEL_CASES) * len(FUEL_ROUTES)):
        tags = _builds("fuel_controller", seed)
        assert sum(t.startswith("expect:") for t in tags) == 1, tags
        case = next(t for t in tags if t.startswith("case:"))
        route = next(t for t in tags if t.startswith("route:"))
        seen.add((case, route))
    assert len(seen) == len(FUEL_CASES) * len(FUEL_ROUTES)


@pytest.mark.parametrize(
    ("tag", "result", "extra", "verdict"),
    [
        ("expect:converged:lean", {"__success__": True}, {"fuel_phi": {"fuel": 0.6}}, None),
        (
            "expect:converged:lean",
            {"__success__": True},
            {"fuel_phi": {"fuel": 1.4}},
            "wrong branch",
        ),
        (
            "expect:converged:rich",
            {"__success__": True},
            {"fuel_phi": {"fuel": 0.9}},
            "wrong branch",
        ),
        ("expect:converged:rich", {"__success__": False}, {}, "refused, expected"),
        ("expect:converged:lean", None, {"exc": "ValueError: x"}, "raised"),
        ("expect:refused:above the peak", {"__success__": True}, {}, "solved, expected refused"),
        (
            "expect:refused:above the peak",
            {"__success__": False},
            {"message": "... is above the peak, 2545 K ..."},
            None,
        ),
        (
            "expect:refused:above the peak",
            {"__success__": False},
            {"message": "no progress"},
            "another reason",
        ),
        (
            "expect:refused_at_build:without passing",
            None,
            {"exc": "ValueError: ... without passing 'b'"},
            None,
        ),
        ("expect:refused_at_build:without passing", {"__success__": False}, {}, "before the solve"),
        ("expect:converged", {"__success__": True}, {}, "unknown expectation"),
    ],
)
def test_a_declared_outcome_is_checked(
    tag: str, result: dict | None, extra: dict, verdict: str | None
) -> None:
    """run.py's check_expected: a miss of a draw's declared outcome is named."""
    from run import check_expected

    got = check_expected({"tags": ["fuel_controller", tag], **extra}, result)
    if verdict is None:
        assert got is None
    else:
        assert got is not None and verdict in got
    assert check_expected({"tags": ["random"]}, {"__success__": False}) is None


def test_a_seed_is_reproducible() -> None:
    """A recorded finding (family, seed) must rebuild the same network."""
    for family, build in FAMILIES.items():
        a, b = build(3)[0], build(3)[0]
        assert sorted(a.elements) == sorted(b.elements), family
        assert sorted(a.nodes) == sorted(b.nodes), family


def test_ejector_draws_stay_inside_their_operating_range() -> None:
    """Incompressible ejector draws put the outlet below the jet pump's
    shut-off head (#619), dPt_o = 2 (Pt_p - Pt_s) b (1 - b/2), b = A_e/A_3,
    and every port chamber is in the draw's regime (it used to default to
    compressible)."""
    seen = 0
    for seed in range(40):
        net, tags = FAMILIES["ejector"](seed)
        regime = tags[1]
        assert {net.nodes[n].regime for n in ("mp", "ms", "mo")} == {regime}
        if regime != "incompressible":
            continue
        seen += 1
        ej = net.elements["ej"]
        b = ej.nozzle_exit_area / ej.mixing_area
        Pp, Ps, Po = (net.nodes[n].Pt for n in ("bp", "bs", "bo"))
        assert Ps < Po < Ps + 2.0 * (Pp - Ps) * b * (1.0 - 0.5 * b)
    assert seen > 5


@pytest.mark.parametrize("seed", [3329152591, 3329152610])
def test_a_loss_between_two_combustors_names_its_theta_source(seed: int) -> None:
    """Fresh seeds drew a PressureLossElement between two combustors, which
    the network refuses without a theta_source (the 0:40 regression set never
    drew one); the generator names the upstream node."""
    net, _ = FAMILIES["random"](seed)
    net.resolve_all_topology()
    net.validate()


def test_fresh_seeds_are_n_consecutive_from_a_printed_offset(
    capsys: pytest.CaptureFixture[str],
) -> None:
    from run import _seeds

    assert _seeds("0:40") == range(0, 40)
    assert _seeds("7") == range(7)
    fresh = _seeds("random:5")
    assert len(fresh) == 5
    assert f"--seeds {fresh.start}:{fresh.stop}" in capsys.readouterr().err
