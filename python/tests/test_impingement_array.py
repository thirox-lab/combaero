"""ImpingementArray: the channel-exit configuration as one object (#465).

It adds rows the solver already knows (validated in
test_impingement_plate_validation.py), so the claims here are about the
assembly: the rows are wired as the validated chain, the wall reaches every
row on its own footprint, and the summary is the rows' own numbers.
"""

from __future__ import annotations

import math

import pytest

import combaero as cb
from combaero.network import (
    ChannelElement,
    ConvectiveSurface,
    FlowNetwork,
    ImpingementArray,
    ImpingementCrossflowElement,
    ImpingementPlateElement,
    NetworkSolver,
    PlenumNode,
    PressureBoundary,
    WallLayer,
)

D = 0.00254
Y_AIR = cb.mole_to_mass(cb.species.dry_air())


def _pb(name: str, Pt: float, Tt: float = 300.0) -> PressureBoundary:
    b = PressureBoundary(name)
    b.Pt, b.Tt, b.Y = Pt, Tt, Y_AIR
    return b


def _array(n_rows: int = 10) -> ImpingementArray:
    return ImpingementArray(
        "ia", n_rows=n_rows, d_jet=D, xn_d=5.0, yn_d=4.0, z_d=2.0, span=0.122, plate_thickness=D
    )


def _net(arr: ImpingementArray, hot: bool = False) -> FlowNetwork:
    net = FlowNetwork()
    net.add_node(_pb("supply", 1.03e5))
    net.add_node(_pb("exit", 1.0e5))
    arr.add_to(net, "supply", "exit")
    if hot:
        net.add_node(_pb("hot_in", 1.02e5, 1200.0))
        net.add_node(_pb("hot_out", 1.0e5, 1200.0))
        L = arr.n_rows * arr.xn_d * D
        net.add_element(
            ChannelElement(
                "hot",
                "hot_in",
                "hot_out",
                length=L,
                diameter=0.05,
                surface=ConvectiveSurface(area=L * arr.span),
            )
        )
        arr.add_wall(net, "w", "hot", [WallLayer(thickness=0.001, conductivity=20.0)])
    return net


def test_the_rows_are_the_validated_chain() -> None:
    arr = _array(3)
    net = _net(arr)
    for i in (1, 2, 3):
        p = net.elements[f"ia__p{i}"]
        x = net.elements[f"ia__x{i}"]
        assert isinstance(p, ImpingementPlateElement) and p.row == i
        assert (p.from_node, p.to_node) == ("supply", f"ia__c{i}")
        assert isinstance(x, ImpingementCrossflowElement)
        assert x.from_node == f"ia__c{i}"
        assert x.to_node == (f"ia__c{i + 1}" if i < 3 else "exit")
        assert x.area == pytest.approx(2.0 * D * 0.122)
        assert x.length == pytest.approx(5.0 * D)
        assert isinstance(net.nodes[f"ia__c{i}"], PlenumNode)


def test_it_reproduces_florschuetz_distribution() -> None:
    """Same check as the hand-built chain, through the configuration."""
    arr = _array()
    res = NetworkSolver(_net(arr)).solve()
    assert res["__success__"], res.get("__message__")
    s = arr.summarize(res["__element_diag__"])
    beta = math.sqrt(2) * cb.FLORSCHUETZ_1981_DEFAULT_CD * (math.pi / 4) / (4.0 * 2.0)
    assert s["jet_flow_nonuniformity"] == pytest.approx(
        math.cosh(9.5 * beta) / math.cosh(0.5 * beta), rel=0.03
    )
    eq8 = cb.crossflow_to_jet_ratio_at_row(4.0, 2.0, cb.FLORSCHUETZ_1981_DEFAULT_CD, 10)
    assert s["Gc_Gj_max"] == pytest.approx(eq8, rel=0.04)


def test_the_summary_is_the_rows_own_numbers() -> None:
    arr = _array(4)
    res = NetworkSolver(_net(arr)).solve()
    diag = res["__element_diag__"]
    s = arr.summarize(diag)
    rows = [diag[f"ia__p{i}"] for i in range(1, 5)]
    assert s["m_dot"] == pytest.approx(sum(r["m_dot"] for r in rows), rel=1e-12)
    assert s["rows_Nu"] == [r["Nu"] for r in rows]
    assert s["n_holes"] == 4 * rows[0]["n_holes"]
    assert s["dP"] == pytest.approx(1.03e5 - 1.0e5, rel=1e-6)
    assert s["rows_Gc_Gj"][0] == 0.0


def test_one_wall_reaches_every_row_on_its_own_footprint() -> None:
    arr = _array(5)
    net = _net(arr, hot=True)
    res = NetworkSolver(net).solve()
    assert res["__success__"], res.get("__message__")
    footprint = net.elements["ia__p1"].surface.area
    for i in range(1, 6):
        w = net.walls[f"w__r{i}"]
        assert (w.element_a, w.element_b) == ("hot", f"ia__p{i}")
        assert w.contact_area == pytest.approx(footprint)
        assert res[f"w__r{i}.Q"] > 0.0  # hot side A -> coolant B
    # Heating the crossflow lowers its density, so the jets' momentum costs
    # more pressure and the supply grows more non-uniform than isothermal.
    iso = NetworkSolver(_net(_array(5))).solve()
    s_hot = arr.summarize(res["__element_diag__"])
    s_iso = _array(5).summarize(iso["__element_diag__"])
    assert s_hot["jet_flow_nonuniformity"] > s_iso["jet_flow_nonuniformity"]


def test_n_rows_must_be_a_positive_integer() -> None:
    with pytest.raises(ValueError, match="n_rows"):
        ImpingementArray(
            "ia", n_rows=0, d_jet=D, xn_d=5, yn_d=4, z_d=2, span=0.1, plate_thickness=D
        )
    with pytest.raises(ValueError, match="n_rows"):
        ImpingementArray(
            "ia", n_rows=2.5, d_jet=D, xn_d=5, yn_d=4, z_d=2, span=0.1, plate_thickness=D
        )


def test_radiation_reaches_every_rows_wall_on_the_named_face() -> None:
    """#1071: the specs go on each row's ThermalWall, hot face and target
    face mapped to sides A/B by `array_is_side_a`; more radiation from the
    gas heats the coolant more."""
    from combaero.network import CasingRadiation, GasRadiation

    gas = GasRadiation(beam_length=0.05, fuel_air_ratio=0.03)
    casing = CasingRadiation(area_ratio=1.0)
    arr = _array(3)
    base = NetworkSolver(_net(arr, hot=True)).solve()
    for side_a in (False, True):
        arr = _array(3)
        net = _net(arr, hot=False)
        net.add_node(_pb("hot_in", 1.02e5, 1200.0))
        net.add_node(_pb("hot_out", 1.0e5, 1200.0))
        L = arr.n_rows * arr.xn_d * D
        net.add_element(
            ChannelElement(
                "hot",
                "hot_in",
                "hot_out",
                length=L,
                diameter=0.05,
                surface=ConvectiveSurface(area=L * arr.span),
            )
        )
        arr.add_wall(
            net,
            "w",
            "hot",
            [WallLayer(thickness=0.001, conductivity=20.0)],
            array_is_side_a=side_a,
            radiation_hot=gas,
            radiation_target=casing,
        )
        res = NetworkSolver(net).solve()
        assert res["__success__"], res.get("__message__")
        for i in range(1, 4):
            w = net.walls[f"w__r{i}"]
            hot, target = ("b", "a") if side_a else ("a", "b")
            assert getattr(w, f"radiation_{hot}") is gas
            assert getattr(w, f"radiation_{target}") is casing
            q = res[f"w__r{i}.Q"]
            assert abs(q) > abs(base[f"w__r{i}.Q"]), i
