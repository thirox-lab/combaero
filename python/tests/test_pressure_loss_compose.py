"""PressureLossElement's loss models compose in one C++ call (#725).

``cb.pressure_loss_compose`` replaced the element's two Python-side kernel
calls (the value at the state's mole fractions, the partials at its raw mass
fractions). The old path is kept here as the oracle, bit for bit, on a grid
of states; and the element's analytic row (the composed partials mapped onto
the unknowns) is checked against the same model run through the differenced
user-callable path, on networks with and without a theta source.
"""

from __future__ import annotations

import copy
import itertools

import numpy as np
import pytest
from test_jacobian_completeness_481 import _jac_err, _premixed_and_dilution_network

import combaero as cb
from combaero.network import (
    ConstantFractionLoss,
    ConstantHeadLoss,
    FlowNetwork,
    LinearThetaFractionLoss,
    LinearThetaHeadLoss,
    NetworkSolver,
    OrificeElement,
    PlenumNode,
    PressureBoundary,
    PressureLossElement,
)

pytestmark = pytest.mark.filterwarnings("ignore::UserWarning")

BUILTINS = [
    ConstantFractionLoss(xi=0.03),
    ConstantFractionLoss(xi=0.03, m_ref=0.5),
    LinearThetaFractionLoss(k=0.01, xi0=0.02),
    LinearThetaFractionLoss(k=0.01, xi0=0.02, m_ref=0.5),
    ConstantHeadLoss(zeta=4.0, area=0.05),
    LinearThetaHeadLoss(k=1.0, zeta0=3.0, area=0.05),
]
IDS = ["cfrac", "cfrac_mref", "ltfrac", "ltfrac_mref", "chead", "lthead"]


def _old(corr: object, theta: float, m: float, Tt: float, Pt: float, Y: list[float]) -> tuple:
    """The replaced path: (xi, dxi_dtheta) from the value call at
    mole_to_mass(X), the other partials from evaluate at the raw Y."""
    # The state's X: Y clamped into [0, 1] and normalised (NetworkMixtureState.X).
    Ys = [max(0.0, min(1.0, y)) for y in Y]
    X = list(cb.mass_to_mole([y / sum(Ys) for y in Ys]))
    Yv = list(cb.mole_to_mass(X))
    if isinstance(corr, (ConstantFractionLoss, LinearThetaFractionLoss)):
        k = corr.k if isinstance(corr, LinearThetaFractionLoss) else 0.0
        base = corr.xi0 if isinstance(corr, LinearThetaFractionLoss) else corr.xi
        th = theta if isinstance(corr, LinearThetaFractionLoss) else 0.0
        mm = m if corr.m_ref is not None else 0.0
        m_ref = corr.m_ref or 0.0
        v = cb.combustor_fraction_loss(th, mm, k, base, m_ref)
        p = v
    else:
        k = corr.k if isinstance(corr, LinearThetaHeadLoss) else 0.0
        base = corr.zeta0 if isinstance(corr, LinearThetaHeadLoss) else corr.zeta
        v = cb.combustor_head_loss(theta, m, Tt, Pt, Yv, k, base, corr.area)
        p = cb.combustor_head_loss(theta, m, Tt, Pt, Y, k, base, corr.area)
    return X, (v.xi, v.dxi_dtheta, p.dxi_dmdot, p.dxi_dTt_in, p.dxi_dPt_in, list(p.dxi_dY))


@pytest.mark.parametrize("corr", BUILTINS, ids=IDS)
def test_the_composition_is_the_old_two_call_path_bit_for_bit(corr: object) -> None:
    air = list(cb.species.dry_air_mass())
    i_ch4 = cb.species_index_from_name("CH4")
    i_co2 = cb.species_index_from_name("CO2")
    lean = [0.97 * y for y in air]
    lean[i_ch4] += 0.03
    # A perturbed Newton iterate: unnormalised, a slightly negative trace.
    off = [1.03 * y for y in lean]
    off[i_co2] = -1e-9
    model = corr.pressure_loss_model()
    head = model.kind == cb.PressureLossKind.Head
    n = 0
    for theta, m, Tt, Pt, Y in itertools.product(
        (0.0, 1.7), (0.2, 0.8, 3.0), (300.0, 650.0), (1.0e5, 1.6e6), (air, lean, off)
    ):
        X, want = _old(corr, theta, m, Tt, Pt, Y)
        inp = cb.PressureLossInputs(
            mdot=m, Tt_in=Tt, Pt_in=Pt, Y_in=Y if head else [], X_in=X if head else [], x=[theta]
        )
        r = cb.pressure_loss_compose(model, inp)
        got = (r.xi, r.dxi_dtheta, r.dxi_dmdot, r.dxi_dTt_in, r.dxi_dPt_in, list(r.dxi_dY))
        assert got == want, (theta, m, Tt, Pt)
        n += 1
    assert n == 72


def _plain_network(corr: object) -> FlowNetwork:
    """A loss between two plenums, no theta source, a CH4-bearing feed."""
    g = FlowNetwork()
    Y = [0.97 * y for y in cb.species.dry_air_mass()]
    Y[cb.species_index_from_name("CH4")] += 0.03
    for name, Pt, Tt in (("a", 2.0e5, 650.0), ("b", 1.0e5, 300.0)):
        bnd = PressureBoundary(name)
        bnd.Pt, bnd.Tt, bnd.Y = Pt, Tt, list(Y)
        g.add_node(bnd)
    g.add_node(PlenumNode("p1"))
    g.add_node(PlenumNode("p2"))
    g.add_element(OrificeElement("o1", "a", "p1", Cd=0.8, diameter=0.03, correlation="fixed"))
    g.add_element(PressureLossElement("loss", "p1", "p2", correlation=corr))
    g.add_element(OrificeElement("o2", "p2", "b", Cd=0.8, diameter=0.03, correlation="fixed"))
    return g


class _Callable:
    """The same built-in model as an opaque user callable: the element takes
    its differenced path."""

    def __init__(self, corr: object) -> None:
        object.__setattr__(self, "_corr", corr)

    def __setattr__(self, name: str, value: object) -> None:
        # The element's inferred area lands on the model it wraps.
        setattr(self._corr, name, value)

    def __call__(self, ctx: object) -> tuple[float, float]:
        return self._corr(ctx)

    def __getattr__(self, name: str) -> object:
        # area, xi, m_ref ... for the seed and area inference; never the model.
        if name == "pressure_loss_model":
            raise AttributeError(name)
        return getattr(self._corr, name)


@pytest.mark.parametrize("theta_source", [False, True], ids=["plain", "combustor"])
@pytest.mark.parametrize("corr", BUILTINS, ids=IDS)
def test_the_analytic_row_matches_the_differenced_path(corr: object, theta_source: bool) -> None:
    build = _premixed_and_dilution_network if theta_source else _plain_network
    rows = []
    for c in (copy.deepcopy(corr), _Callable(copy.deepcopy(corr))):
        g = build(c)
        g.resolve_all_topology()
        s = NetworkSolver(g)
        x0 = np.array(s._build_x0(), dtype=float)
        rng = np.random.default_rng(3)
        out = []
        for k in range(3):
            x = x0 if k == 0 else x0 * (1.0 + 0.02 * rng.standard_normal(x0.size))
            F, J = s._residuals_and_jacobian(x)
            i = [j for j, (o, _) in enumerate(s._row_kinds(x)) if o == "PressureLossElement 'loss'"]
            assert len(i) == 1
            out.append((np.asarray(F)[i[0]], J.toarray()[i[0]], np.abs(x)))
        rows.append(out)
    for (F_a, J_a, ax), (F_d, J_d, _) in zip(*rows, strict=True):
        # The value is the same C++ call either way.
        assert F_a == F_d
        # Entry by entry, scaled as the solver sees it (J_ij |x_j|): the
        # differences carry ~1e-8 relative noise, a floor at 1e-10 of the row.
        Ja, Jd = J_a * ax, J_d * ax
        floor = 1e-10 * max(np.max(np.abs(Jd)), 1e-300)
        assert np.all(np.abs(Ja - Jd) <= 1e-7 * np.abs(Jd) + floor)
        assert np.count_nonzero(Ja) == np.count_nonzero(Jd)


@pytest.mark.parametrize("theta_source", [False, True], ids=["plain", "combustor"])
@pytest.mark.parametrize("corr", BUILTINS, ids=IDS)
def test_the_loss_row_matches_differences_of_the_network(corr: object, theta_source: bool) -> None:
    """The whole mapping (theta chain onto the burned and unburned states,
    m_dot sign, composition columns) against differences of the residuals."""
    build = _premixed_and_dilution_network if theta_source else _plain_network
    g = build(copy.deepcopy(corr))
    g.resolve_all_topology()
    s = NetworkSolver(g)
    x0 = np.array(s._build_x0(), dtype=float)
    assert _jac_err(s, x0, owner="PressureLossElement 'loss'") < 1e-6
