"""Fidelity of cb's Bordbar et al. (2021) WSGG against the paper's own slab
results (#1071).

Bordbar, Coelho, Fraga, Franca & Hostikka (2021), Int. J. Heat Mass Transfer
173, 121207, Section 4.2, Figs. 3-5: a 1 m slab between black walls at 300 K,
T(x) = 300 K + 1500 K sin^2(pi x/L) (Eq. 18), three composition profiles
(Eqs. 19-21), at 1, 30 and 0.5 atm. The figures plot the radiative source Sr
and flux qr from their WSGG model (lines) and from line-by-line integration
(squares). The lines are the WSGG model's own output, so reproducing them is
a FIDELITY check of the coefficients and of Eqs. 7-11 as implemented in
``cb.wsgg_bordbar2021_gray_gases``. The squares are scored only to show the
model-vs-LBL error the paper reports (Table 2); they say nothing about our
implementation.

The RTE is solved as the paper describes (Section 3): discrete ordinates on
200 cells, each gray gas plus the transparent one solved separately and
summed (Eqs. 12-17). Three choices are ours, all discretisation: an exact
exponential cell march in place of the paper's finite differencing,
Gauss-Legendre ordinates in place of Lathrop & Carlson's set, and 16 of them
per hemisphere instead of the paper's 6. The paper reports its 12 ordinates
as converged; with Gauss-Legendre they are not quite (max / mean deviation
2.44 / 0.93 % at 4 per hemisphere, 1.42 / 0.44 % at 6, 1.30 / 0.29 % at 8,
1.25 / 0.28 % at 16), so the converged plateau is scored -- the model's own
answer, not a quadrature residual.

The curves were extracted from the article PDF's vector paths (see
``extractions/bordbar_2021_slab.md``).

Run: ``uv run python -m validation.radiation.bordbar2021_slab``
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass

import numpy as np

import combaero as cb
from validation.radiation import DATA_ROOT

DATA = DATA_ROOT / "bordbar2021_slab"

# Section 4.2: slab length [m], wall temperature [K], cells. Ordinates per
# hemisphere: converged (see the module docstring), the paper's is 6.
SLAB_LENGTH = 1.0
T_WALL = 300.0
N_CELLS = 200
N_MU_HEMISPHERE = 16


def profiles(case: int, x_hat: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Eqs. 18-21: T [K], x_H2O, x_CO2 and Mr = x_H2O / x_CO2 at x_hat = x / L
    (Case 1 is Mr = 2 everywhere, the walls included, where both vanish)."""
    s2 = np.sin(np.pi * x_hat) ** 2
    T = 300.0 + 1500.0 * s2
    if case == 1:
        xh = 0.2 * s2
        xc = xh / 2.0
    elif case == 2:
        xh = 0.1 * (1.0 + s2)
        xc = 1.0 - xh
    elif case == 3:
        xh = (1.0 + s2) / 3.0
        xc = 1.0 - xh
    else:
        raise ValueError(f"case must be 1, 2 or 3, not {case}")
    Mr = np.full_like(x_hat, 2.0) if case == 1 else xh / xc
    return T, xh, xc, Mr


def _gray_gases(p_atm: float, T: np.ndarray, Mr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a = np.zeros((len(T), 4))
    kappa = np.zeros((len(T), 4))
    for n in range(len(T)):
        g = cb.wsgg_bordbar2021_gray_gases(p_atm, float(T[n]), float(Mr[n]))
        a[n] = g.a
        kappa[n] = g.kappa
    return a, kappa


def solve_slab(p_atm: float, case: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Discrete-ordinates WSGG solve: (x_hat cells, Sr [kW/m^3], x_hat faces,
    qr [kW/m^2])."""
    xf = np.linspace(0.0, 1.0, N_CELLS + 1)
    xc_hat = 0.5 * (xf[1:] + xf[:-1])
    dx = SLAB_LENGTH / N_CELLS
    T, xh, xc, Mr = profiles(case, xc_hat)
    a, kp = _gray_gases(p_atm, T, Mr)
    # kappa_i = kappa_p,i (p_H2O + p_CO2), p in atm (Eq. 7)
    k = kp * ((xh + xc) * p_atm)[:, None]
    sigma = cb.STEFAN_BOLTZMANN
    Ib = sigma * T**4 / math.pi
    # Black walls: band i of the wall emits a_i(T_wall) Ib(T_wall), at the
    # wall's own composition.
    Mr_w = profiles(case, np.array([0.0, 1.0]))[3]
    aw, _ = _gray_gases(p_atm, np.array([T_WALL, T_WALL]), Mr_w)
    Ibw = sigma * T_WALL**4 / math.pi
    mu, w = np.polynomial.legendre.leggauss(2 * N_MU_HEMISPHERE)
    sel = mu > 0.0
    mu = mu[sel]
    w = w[sel] / w[sel].sum()
    qf = np.zeros(N_CELLS + 1)
    Sr = np.zeros(N_CELLS)
    for g in range(5):
        if g < 4:
            ag, kg, awg = a[:, g], k[:, g], aw[:, g]
        else:
            ag, kg, awg = 1.0 - a.sum(1), np.zeros(N_CELLS), 1.0 - aw.sum(1)
        Sg = ag * Ib
        G = np.zeros(N_CELLS)
        for m, wm in zip(mu, w, strict=True):
            tau = kg * dx / m
            e = np.exp(-tau)
            avg = np.where(tau > 1e-12, (1.0 - e) / np.maximum(tau, 1e-300), 1.0)
            Ip = np.zeros(N_CELLS + 1)
            Ip[0] = awg[0] * Ibw
            for c in range(N_CELLS):
                Ip[c + 1] = Ip[c] * e[c] + Sg[c] * (1.0 - e[c])
            Im = np.zeros(N_CELLS + 1)
            Im[N_CELLS] = awg[1] * Ibw
            for c in range(N_CELLS - 1, -1, -1):
                Im[c] = Im[c + 1] * e[c] + Sg[c] * (1.0 - e[c])
            Ip_avg = Sg + (Ip[:-1] - Sg) * avg
            Im_avg = Sg + (Im[1:] - Sg) * avg
            qf += 2.0 * math.pi * wm * m * (Ip - Im)
            G += 2.0 * math.pi * wm * (Ip_avg + Im_avg)
        Sr += kg * (G - 4.0 * math.pi * Sg)
    return xc_hat, Sr / 1e3, xf, qf / 1e3


def load(name: str) -> dict[tuple[float, int, str], tuple[np.ndarray, np.ndarray]]:
    rows: dict[tuple[float, int, str], list[tuple[float, float]]] = {}
    with open(DATA / name, newline="") as f:
        for r in csv.DictReader(f):
            key = (float(r["p_atm"]), int(r["case"]), r["quantity"])
            rows.setdefault(key, []).append((float(r["x_hat"]), float(r["value"])))
    return {k: (np.array([x for x, _ in v]), np.array([y for _, y in v])) for k, v in rows.items()}


@dataclass
class CurveScore:
    p_atm: float
    case: int
    quantity: str
    n: int
    max_norm: float  # max |ours - paper| / max |paper|
    mean_norm: float
    lbl_n: int
    lbl_max_norm: float  # max |ours - LBL| / max |LBL| (the model's own error)
    lbl_mean_norm: float


def score() -> list[CurveScore]:
    curves = load("wsgg_curves.csv")
    lbl = load("lbl_points.csv")
    cache: dict[tuple[float, int], tuple[np.ndarray, ...]] = {}
    out = []
    for (p, case, q), (x, v) in curves.items():
        if (p, case) not in cache:
            cache[(p, case)] = solve_slab(p, case)
        xc, Sr, xf, qr = cache[(p, case)]
        xs, ys = (xc, Sr) if q == "Sr" else (xf, qr)
        dev = np.abs(np.interp(x, xs, ys) - v) / np.abs(v).max()
        xl, vl = lbl[(p, case, q)]
        dl = np.abs(np.interp(xl, xs, ys) - vl) / np.abs(vl).max()
        out.append(
            CurveScore(p, case, q, len(x), float(dev.max()), float(dev.mean()), len(xl), float(dl.max()), float(dl.mean()))
        )
    return out


def main() -> None:
    res = score()
    print("Bordbar 2021 Figs. 3-5, FIDELITY (cb vs the paper's own WSGG curves)")
    print(" p_atm case  q    n   max%   mean%  | vs LBL (model error): n  max%  mean%")
    for s in res:
        print(
            f"{s.p_atm:6g} {s.case:4d} {s.quantity:>3s} {s.n:4d} {100 * s.max_norm:6.2f} {100 * s.mean_norm:6.2f}"
            f"  | {s.lbl_n:3d} {100 * s.lbl_max_norm:6.2f} {100 * s.lbl_mean_norm:6.2f}"
        )
    n = sum(s.n for s in res)
    mx = max(s.max_norm for s in res)
    mean = sum(s.mean_norm * s.n for s in res) / n
    print(f"ALL n={n}: max {100 * mx:.2f}%, mean {100 * mean:.2f}%")


if __name__ == "__main__":
    main()
