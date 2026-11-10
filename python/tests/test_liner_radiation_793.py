"""#793: liner radiation, Lefebvre & Ballal (2010) Ch. 8 -- gas-to-wall R1
(Eq. 8.7) with the luminous Reeves emissivity (Eq. 8.11) and wall-to-casing
R2 (Eq. 8.17) on a ThermalWall's faces, opt-in.

FIDELITY on printed quantities: the book's three chained worked examples
(uncooled primary-zone liner pp. 325-328, film-cooled at x/s = 18 pp.
338-340, 0.5 mm TBC pp. 350-351), solved from the printed inputs with the
printed convective coefficients. It tests our implementation of the book's
chain, not the physics of eps_g (Reeves' fit is not checked against data
here).

The network case (combaero's own convection on the uncooled liner) is a
method comparison, not fidelity: its hot-side h is Gnielinski on combaero's
products, not the book's 0.017 Re^0.8 primary-zone constant.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pytest

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, SmoothModel
from combaero.network import (
    CasingRadiation,
    ChannelElement,
    FlowNetwork,
    GasRadiation,
    MassFlowBoundary,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    ThermalWall,
    WallLayer,
)

AIR = list(cb.mole_to_mass(cb.species.dry_air()))

# The book's example inputs (Lefebvre & Ballal 3rd ed. p. 325-327).
TG = 2280.0
P3 = 3040e3
Q_PZ = 0.0588  # nominal: eps_g = 0.61 needs it (Tg used 0.85 x 0.0588)
L_B = 0.0792  # 0.6 D_L
LUM = 1.7
EPS_W = 0.7
EPS_C = 0.4
AREA_RATIO = 0.1344 / 0.192
T3 = 880.0
H1 = 565.0  # Eq. 8.28, C1 = 565 (Tg - Tw1)
H2 = 921.0  # Eq. 8.29 (printed 291: erratum)
R_METAL = 0.0012 / 26.0
R_TBC = 0.0005 / 2.6
H1_FILM, T_AW_FILM = 1926.0, 1176.0  # p. 340, C1 = 1926 (1176 - Tw1)


def _book_wall(layers: list[WallLayer]) -> ThermalWall:
    return ThermalWall(
        id="w",
        element_a="gas",
        element_b="annulus",
        layers=layers,
        contact_area=1.0,
        radiation_a=GasRadiation(
            beam_length=L_B,
            emissivity_wall=EPS_W,
            luminosity=LUM,
            T_gas=TG,
            fuel_air_ratio=Q_PZ,
            pressure=P3,
        ),
        radiation_b=CasingRadiation(
            area_ratio=AREA_RATIO, emissivity_wall=EPS_W, emissivity_casing=EPS_C, T_casing=T3
        ),
    )


# --- fidelity: printed intermediates ---------------------------------------


def test_book_intermediates() -> None:
    """eps_g = 0.61 and R1 = 794 460 - 0.0032 Tw^2.5 (Eq. 8.26); R2 = 2.29
    (Tw/100)^4 - 13 715 (Eq. 8.27). The book rounds eps_g to two digits
    before forming R1's constant, so that constant agrees to the rounding
    (0.61 / 0.6092 = 1.0014)."""
    e = cb.gas_emissivity_luminous(TG, P3, Q_PZ, L_B, LUM)
    assert e.eps == pytest.approx(0.61, abs=0.005)
    r_hot = cb.liner_gas_radiation(TG, 0.0, e.eps, EPS_W)
    assert pytest.approx(794460.0, rel=2e-3) == r_hot.R1
    coef = r_hot.R1 / TG**2.5  # the Tw^2.5 coefficient
    assert coef == pytest.approx(0.0032, rel=0.01)
    F = cb.liner_casing_exchange_factor(EPS_W, EPS_C, AREA_RATIO)
    assert pytest.approx(2.29, abs=0.005) == F * 1e8
    R2_T3, _, _ = cb.liner_casing_radiation(T3, 0.0, EPS_W, EPS_C, AREA_RATIO)
    assert pytest.approx(13715.0, rel=1e-3) == R2_T3
    # The book's sigma is 5.67e-8; CODATA's is 0.007 % higher.
    assert pytest.approx(5.67e-8, rel=1e-4) == cb.STEFAN_BOLTZMANN


# --- fidelity: the chained worked examples, through ThermalWall -----------


@pytest.mark.parametrize(
    ("case", "h1", "T_aw", "layers", "printed", "tol"),
    [
        # From the inputs the uncooled liner is 1.6 / 1.4 K above print: the
        # book solves with eps_g rounded to 0.61 and C1 = 1 280 500 - 562 Tw1.
        ("uncooled", H1, TG, [WallLayer(0.0012, 26.0)], (1640.0, 1603.0), 2.0),
        ("film x/s=18", H1_FILM, T_AW_FILM, [WallLayer(0.0012, 26.0)], (1283.0, 1265.0), 1.0),
        (
            "film + 0.5 mm TBC",
            H1_FILM,
            T_AW_FILM,
            [WallLayer(0.0005, 2.6), WallLayer(0.0012, 26.0)],
            (1304.0, 1236.0, 1220.0),
            1.0,
        ),
    ],
)
def test_lefebvre_worked_examples(
    case: str, h1: float, T_aw: float, layers: list, printed: tuple, tol: float
) -> None:
    """Fidelity on printed quantities: between two reservoirs the wall's
    Q/A is C1 + R1 = C2 + R2 = K, the book's Eq. 8.25 (8.46 with the TBC)."""
    w = _book_wall(layers)
    r = w.compute_coupling(h1, T_aw, 1.0, H2, T3, 1.0)
    prof = w._last_profile
    assert len(prof) == len(printed)
    for got, want in zip(prof, printed, strict=True):
        assert got == pytest.approx(want, abs=tol), (case, prof)
    # Each face's balance, with the kernels on their own.
    e = cb.gas_emissivity_luminous(TG, P3, Q_PZ, L_B, LUM).eps
    R1 = cb.liner_gas_radiation(TG, prof[0], e, EPS_W).R1
    R2 = cb.liner_casing_radiation(prof[-1], T3, EPS_W, EPS_C, AREA_RATIO)[0]
    assert pytest.approx(h1 * (T_aw - prof[0]) + R1, rel=1e-10) == r.Q
    assert pytest.approx(H2 * (prof[-1] - T3) + R2, rel=1e-10) == r.Q
    assert r.radiation["q_rad_a"] == pytest.approx(R1, rel=1e-10)
    assert r.radiation["q_rad_b"] == pytest.approx(R2, rel=1e-10)


# --- the coupling's derivatives ---------------------------------------------


def _gas_radiation(model: str) -> GasRadiation:
    """The hot face's spec: Eq. 8.11 with a luminosity, or the Bordbar 2021
    WSGG (#1071), which reads the composition through x_H2O and x_CO2."""
    if model == "bordbar2021":
        return GasRadiation(beam_length=0.08, emissivity_model="bordbar2021")
    return GasRadiation(beam_length=0.08, luminosity=1.7)


def _products(q: float) -> list[float]:
    """n-heptane burnt completely at fuel/air ratio q (mass fractions)."""
    Y = [v / (1.0 + q) for v in AIR]
    Y[cb.species_index_from_name("NC7H16")] = q / (1.0 + q)
    X = cb.complete_combustion_to_CO2_H2O(list(cb.mass_to_mole(Y)))
    return list(cb.mole_to_mass(X))


@pytest.mark.parametrize(
    "wall",
    [
        {"contact_area": 0.05},
        {"area_a": 0.03, "area_b": 0.04},
    ],
)
@pytest.mark.parametrize(
    "layers",
    [
        [WallLayer(0.0012, 26.0)],
        [WallLayer(0.0003, 1.0, material="ysz"), WallLayer(0.0012, 20.0, material="hastelloy_x")],
    ],
    ids=["constant_k", "k_of_T"],
)
@pytest.mark.parametrize("C", [(math.inf, math.inf), (40.0, 25.0)], ids=["reservoirs", "streams"])
@pytest.mark.parametrize("model", ["lefebvre", "bordbar2021"])
def test_coupling_derivatives(wall: dict, layers: list, C: tuple, model: str) -> None:
    """Every input derivative of a radiating wall's Q, through the surface
    temperatures' fixed point (and the k(T) one), against central
    differences: h, T_aw per side and the radiation's stream Tt, P and Y."""
    Y = _products(0.04)
    w = ThermalWall(
        "w",
        "a",
        "b",
        layers=layers,
        radiation_a=_gas_radiation(model),
        radiation_b=CasingRadiation(area_ratio=0.75),
        **wall,
    )
    args = [900.0, 1500.0, 0.05, 700.0, 880.0, 0.05, C[0], C[1], True]
    gas_a = (1600.0, 2.0e6, Y)
    gas_b = (870.0, 2.0e6, AIR)

    def Q(a: list, ga: tuple = gas_a, gb: tuple = gas_b) -> float:
        return w.compute_coupling(*a, gas_a=ga, gas_b=gb).Q

    r = w.compute_coupling(*args, gas_a=gas_a, gas_b=gas_b)
    for i, name in ((0, "dQ_dh_a"), (1, "dQ_dT_aw_a"), (3, "dQ_dh_b"), (4, "dQ_dT_aw_b")):
        e = 1e-6 * abs(args[i])
        hi, lo = list(args), list(args)
        hi[i] += e
        lo[i] -= e
        fd = (Q(hi) - Q(lo)) / (2 * e)
        assert getattr(r, name) == pytest.approx(fd, rel=1e-6, abs=1e-9 * abs(r.Q)), name
    for side, g in (("a", gas_a), ("b", gas_b)):
        for j, name in ((0, "Tt"), (1, "P")):
            e = 1e-6 * g[j]
            hi = tuple(v + e if k == j else v for k, v in enumerate(g))
            lo = tuple(v - e if k == j else v for k, v in enumerate(g))
            if side == "a":
                fd = (Q(args, ga=hi) - Q(args, ga=lo)) / (2 * e)
            else:
                fd = (Q(args, gb=hi) - Q(args, gb=lo)) / (2 * e)
            got = getattr(r, f"rad_dQ_d{name}_{side}")
            assert got == pytest.approx(fd, rel=1e-6, abs=1e-9 * abs(r.Q)), (side, name)
    # The composition enters side A through q (Lefebvre) or x_H2O, x_CO2
    # (Bordbar); side B's casing reads none.
    assert r.rad_dQ_dY_b is None
    dY = np.asarray(r.rad_dQ_dY_a)
    for k in range(len(Y)):
        if Y[k] == 0.0:
            continue
        e = 1e-6
        hi, lo = list(Y), list(Y)
        hi[k] += e
        lo[k] -= e
        fd = (Q(args, ga=(gas_a[0], gas_a[1], hi)) - Q(args, ga=(gas_a[0], gas_a[1], lo))) / (2 * e)
        assert dY[k] == pytest.approx(fd, rel=1e-5, abs=1e-9 * abs(r.Q)), k


def test_radiation_adds_heat_and_none_changes_nothing() -> None:
    """A gray gas only adds heat to a cooler wall; without a spec the
    coupling is the convective kernel's own result object."""
    base = ThermalWall("w", "a", "b", layers=[WallLayer(0.0012, 26.0)], contact_area=1.0)
    r0 = base.compute_coupling(565.0, TG, 1.0, H2, T3, 1.0)
    assert isinstance(r0, cb.WallCouplingResult)
    r1 = _book_wall([WallLayer(0.0012, 26.0)]).compute_coupling(565.0, TG, 1.0, H2, T3, 1.0)
    assert r1.Q > r0.Q
    # A zero beam length is refused, not silently non-radiating.
    with pytest.raises(ValueError):
        GasRadiation(beam_length=0.0)
    with pytest.raises(ValueError):
        CasingRadiation(area_ratio=0.7, emissivity_casing=1.5)
    with pytest.raises(ValueError):
        ThermalWall("w", "a", "b", radiation_a="gas")  # type: ignore[arg-type]


# --- in a network ------------------------------------------------------------


def _mixing_net(
    rad: bool = True, layers: list | None = None, model: str = "lefebvre", **wall: object
) -> FlowNetwork:
    """Hot products mixed with air in a plenum feed the radiating duct, so
    the wall's gas Tt, P and composition all move with the unknowns; both
    ducts deliver into interior plenums, whose T carries the wall heat."""
    n = FlowNetwork()
    n.add_node(PressureBoundary("air_in", Pt=3.10e6, Tt=880.0, Y=AIR))
    n.add_node(PressureBoundary("prod_in", Pt=3.10e6, Tt=2400.0, Y=_products(0.08)))
    n.add_node(PlenumNode("mix"))
    n.add_node(PlenumNode("h_mid"))
    n.add_node(PressureBoundary("h_out", Pt=3.0e6, Tt=2000.0, Y=_products(0.08)))
    n.add_node(MassFlowBoundary("c_in", m_dot=0.3, Tt=880.0, Y=AIR))
    n.add_node(PlenumNode("c_mid"))
    n.add_node(PressureBoundary("c_out", Pt=3.0e6, Tt=880.0, Y=AIR))
    n.add_element(ChannelElement("a1", "air_in", "mix", length=0.1, diameter=0.03))
    n.add_element(ChannelElement("p1", "prod_in", "mix", length=0.1, diameter=0.05))
    n.add_element(
        ChannelElement(
            "hot",
            "mix",
            "h_mid",
            length=0.1,
            diameter=0.13,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    n.add_element(ChannelElement("h_exit", "h_mid", "h_out", length=0.1, diameter=0.06))
    n.add_element(
        ChannelElement(
            "cold",
            "c_in",
            "c_mid",
            length=0.1,
            diameter=0.03,
            surface=ConvectiveSurface(model=SmoothModel()),
        )
    )
    n.add_element(ChannelElement("c_exit", "c_mid", "c_out", length=0.1, diameter=0.02))
    kw: dict = {}
    if rad:
        kw = {
            "radiation_a": _gas_radiation(model),
            "radiation_b": CasingRadiation(area_ratio=0.7),
        }
    n.add_wall(
        ThermalWall("w", "hot", "cold", layers=layers or [WallLayer(0.0012, 26.0)], **kw, **wall)
    )
    return n


def _solve(n: FlowNetwork) -> tuple[NetworkSolver, dict]:
    s = NetworkSolver(n)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = s.solve()
    assert r["__success__"], r["__message__"]
    return s, r


@pytest.mark.parametrize(
    ("layers", "wall"),
    [
        (None, {}),
        (
            [
                WallLayer(0.0003, 1.0, material="ysz"),
                WallLayer(0.0012, 20.0, material="hastelloy_x"),
            ],
            {},
        ),
        (None, {"area_a": 0.04, "area_b": 0.05}),
    ],
    ids=["constant_k", "k_of_T", "per_side"],
)
@pytest.mark.parametrize("model", ["lefebvre", "bordbar2021"])
def test_network_relay_is_exact(layers: list | None, wall: dict, model: str) -> None:
    """The two wall-heated nodes' temperature sensitivities (the relay the
    Jacobian is built from) against central differences of the propagated
    temperatures, in every unknown. Dropping any one radiation term (gas Tt,
    P, composition) moves them by 0.03-3 K; this holds them to 1e-5 K."""
    s, r = _solve(_mixing_net(True, layers, model, **wall))
    assert r["w.Q_rad_a"] > 0.0 and r["w.Q_rad_b"] > 0.0
    x = np.array(r["__x_solution__"])
    relay = s._propagate_states(x, with_relay=True)

    def T(v: np.ndarray, nid: str) -> float:
        s._propagate_states(v, with_relay=False)
        return float(s._derived_states[nid][0])

    for nid in ("h_mid", "c_mid"):
        for j in range(len(x)):
            h = 1e-6 * max(abs(x[j]), 1e-3)
            e = np.zeros_like(x)
            e[j] = h
            fd = (T(x + e, nid) - T(x - e, nid)) / (2 * h)
            an = relay[nid].get(j, {}).get("T", 0.0)
            assert abs(an - fd) * abs(x[j]) < 1e-5, (nid, s.unknown_names[j], an, fd)


def test_network_reports_radiation_and_closes_energy() -> None:
    """The radiative fluxes are reported per face; the heat is booked to the
    two streams, so the network's energy balance still closes."""
    s, r = _solve(_mixing_net(True))
    for key in ("q_rad_a", "Q_rad_a", "T_surface_a", "eps_g_a", "q_rad_b", "Q_rad_b"):
        assert f"w.{key}" in r
    assert r["w.T_surface_a"] == pytest.approx(r["w.T_interface"][0])
    assert 0.0 < r["w.eps_g_a"] < 1.0
    _, r0 = _solve(_mixing_net(False))
    assert "w.q_rad_a" not in r0
    assert r["w.Q"] > r0["w.Q"]


def test_book_liner_as_a_network() -> None:
    """The uncooled example as a combaero network: products of n-heptane at
    q = 0.05 and 2280 K through a 0.132 m liner, 7.074 kg/s of 880 K air in
    the annulus (Dh 0.0576 m, 0.01476 m^2), a short segment so neither stream
    heats up. Radiation from the book's inputs; convection is combaero's own
    (Gnielinski): h_a 523 against the book's 565 (0.017 Re^0.8 with k_g 14 %
    above combaero's), h_b 911 against 921. NOT fidelity: a method
    comparison, Tw1/Tw2 1632.9/1596.6 K against the printed 1640/1603.

    What it does test: the network's wall is the local balance of its own
    h's and the book's radiation, to round-off."""
    L = 0.01
    A = math.pi * 0.132 * L
    n = FlowNetwork()
    n.add_node(MassFlowBoundary("g_in", m_dot=2.62, Tt=TG, Y=_products(0.05)))
    n.add_node(PressureBoundary("g_out", Pt=P3, Tt=TG, Y=_products(0.05)))
    n.add_node(MassFlowBoundary("a_in", m_dot=7.074, Tt=T3, Y=AIR))
    n.add_node(PressureBoundary("a_out", Pt=P3, Tt=T3, Y=AIR))
    n.add_element(
        ChannelElement(
            "liner",
            "g_in",
            "g_out",
            length=L,
            diameter=0.132,
            surface=ConvectiveSurface(area=A, model=SmoothModel()),
        )
    )
    n.add_element(
        ChannelElement(
            "annulus",
            "a_in",
            "a_out",
            length=L,
            diameter=math.sqrt(4 * 0.01476 / math.pi),
            Dh=0.0576,
            surface=ConvectiveSurface(area=math.pi * 0.1344 * L, model=SmoothModel()),
        )
    )
    n.add_wall(
        ThermalWall(
            "w",
            "liner",
            "annulus",
            layers=[WallLayer(0.0012, 26.0)],
            contact_area=A,
            radiation_a=GasRadiation(
                beam_length=L_B, emissivity_wall=EPS_W, luminosity=LUM, fuel_air_ratio=Q_PZ
            ),
            radiation_b=CasingRadiation(
                area_ratio=AREA_RATIO, emissivity_wall=EPS_W, emissivity_casing=EPS_C
            ),
        )
    )
    _, r = _solve(n)
    Tw1, Tw2 = r["w.T_interface"]
    assert Tw1 == pytest.approx(1632.9, abs=1.0)
    assert Tw2 == pytest.approx(1596.6, abs=1.0)
    # The reported surfaces are the local balance at the network's own h,
    # T_aw and the stream's Tt (the fed MassFlowBoundary's 2280 K) and P.
    R1 = cb.liner_gas_radiation(TG, Tw1, r["w.eps_g_a"], EPS_W).R1
    R2 = cb.liner_casing_radiation(Tw2, T3, EPS_W, EPS_C, AREA_RATIO)[0]
    q_wall = (Tw1 - Tw2) / R_METAL
    assert r["w.h_a"] * (r["w.T_aw_a"] - Tw1) + R1 == pytest.approx(q_wall, rel=1e-9)
    assert r["w.h_b"] * (Tw2 - r["w.T_aw_b"]) + R2 == pytest.approx(q_wall, rel=1e-6)
    # Q itself is the streams' epsilon-NTU exchange on those conductances:
    # 0.05 % below the inlet-state balance at this segment's NTU.
    assert r["w.Q"] / A == pytest.approx(q_wall, rel=1e-3)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(pytest.main([str(Path(__file__)), "-q"]))
