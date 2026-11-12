#!/usr/bin/env python3
"""Combustor liner convective cooling design with bypass sweep.

Flow topology (counterflow convective cooling):

  Compressor exit (T2, P2)
       |
       +--[f_cool]--> liner cooling channels --> T_cool_exit  (coolant heats up)
       |
       +--[1-f_cool]--> bypass (at T2)
       |
       [mix: coolant_exit + bypass] --> oxidizer stream at T_mix
       |
       [set_fuel_stream_for_Tad(COT_target)] --> fuel flow adjusted to hit COT
       |
       complete_combustion --> COT = const (adiabatic)

Energy balance note:
  Wall heat loss lowers effective COT; set_fuel_stream_for_Tad compensates by
  increasing fuel flow to maintain COT = const. More cooling -> more fuel needed.

Counterflow pressure budget (both in series):
  dP_total/P = dP_burner/P + dP_cool/P

Bypass sweep: vary the cooling fraction f_cool for two channel geometries
and print, per point, T_metal, the total dP/P (burner + coolant) and the
thermal performance factor eta_thp, plus the smallest f_cool that keeps
T_metal < 800 C. The example reports these numbers; it does not judge
which geometry is the better choice:
  - Smooth
  - Rib-roughened  (Han et al. 1988, e/D=0.07, P/e=8, alpha=60 deg)
    Note: channel_ribbed uses the geometry-only rib_enhancement_factor internally.
    For standalone high-Re Nu/Nu0 and f/f0 use rib_enhancement_factor_high_re
    and rib_friction_multiplier_high_re (Singh & Ekkad 2017, Re=30k-400k).

Wall configurations compared:
  - Bare Inconel 718
  - Inconel 718 + YSZ APS TBC (fresh and aged 1000 h)
  - Inconel 718 + YSZ EB-PVD TBC

Key combaero tools used:
  - humid_air_composition              (realistic coolant composition)
  - mix + set_fuel_stream_for_Tad      (closed energy loop at fixed COT)
  - channel_smooth / channel_ribbed / channel_dimpled
  - thermal_performance_factor         (Webb & Eckert 1972 eta = Nu_r/f_r^(1/3))
  - wall_temperature_profile           (multi-layer T profile)
  - k_inconel718 / k_tbc_ysz          (temperature-dependent material k)
  - list_materials                     (material database)
  - nusselt_dittus_boelter / htc_from_nusselt  (hot-side HTC)
"""

import matplotlib.pyplot as plt
import numpy as np
from plot_utils import show_or_save

import combaero as cb
from combaero.heat_transfer import ConvectiveSurface, RibbedModel

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def hot_side_htc(T_hot: float, P: float, X: list, D_h: float, u: float) -> float:
    """Hot-gas-side HTC via Dittus-Boelter on combustor annulus."""
    rho = cb.density(T_hot, P, X)
    mu = cb.viscosity(T_hot, P, X)
    k = cb.thermal_conductivity(T_hot, P, X)
    Re = rho * u * D_h / mu
    Nu = cb.nusselt_dittus_boelter(Re, cb.prandtl(T_hot, P, X), heating=False)
    return cb.htc_from_nusselt(Nu, k, D_h)


def channel_result(
    name: str, T: float, P: float, X: list, u: float, D_h: float, L: float
) -> cb.ChannelResult:
    if name == "Smooth":
        return cb.channel_smooth(T, P, X, u, D_h, L)
    if name == "Ribbed":
        # Ribs are evaluated through a ConvectiveSurface rather than a bare
        # correlation call: the correlation gives the RIBBED SIDE, and
        # combining it with the smooth walls depends on how many walls carry
        # ribs, which is the element's job. See issue #334.
        surface = ConvectiveSurface(
            area=1.0,
            model=RibbedModel(e_D=0.07, p_e=8.0, alpha_deg=60.0, W_H=1.0, n_ribbed_walls=2),
        )
        return surface.htc_and_T(T, P, X, u, D_h, L)
    raise ValueError(f"Unknown channel type: {name}")


def solve_operating_point(
    f_cool: float,
    mdot_total: float,
    T2: float,
    P2: float,
    X_cool: list,
    X_air: list,
    X_ch4: list,
    COT_target: float,
    dP_burner_frac: float,
    D_ch: float,
    L_ch: float,
    N_ch: int,
    D_ann: float,
    u_hot_ann: float,
    ch_name: str,
    wall_layers: list[tuple[float, float]],
) -> dict:
    """Solve the closed energy loop for one (f_cool, channel, wall) point.

    Iteration: T_hot estimate -> h_hot -> q_wall -> T_cool_exit -> T_mix
               -> set_fuel_stream_for_Tad -> COT -> update T_hot.
    Converges in 2-3 iterations.
    """
    A_ch_total = N_ch * np.pi * (D_ch / 2) ** 2
    A_liner = np.pi * D_ch * L_ch * N_ch  # total wetted inner area

    mdot_cool = f_cool * mdot_total
    mdot_bypass = (1.0 - f_cool) * mdot_total
    rho_cool = cb.density(T2, P2, X_cool)
    u_cool = mdot_cool / (rho_cool * A_ch_total)
    cp_cool = cb.cp_mass(T2, X_cool)

    P_hot = P2 * (1.0 - dP_burner_frac)
    T_hot_est = COT_target

    r_ch = channel_result(ch_name, T2, P2, X_cool, u_cool, D_ch, L_ch)
    h_cool = r_ch.h

    # Smooth-channel baseline for thermal performance factor
    r_smooth = cb.channel_smooth(T2, P2, X_cool, u_cool, D_ch, L_ch)
    eta_thp = cb.thermal_performance_factor(
        Nu_ratio=r_ch.Nu / r_smooth.Nu,
        f_ratio=max(r_ch.f / r_smooth.f, 1.0),
    )

    for _ in range(4):
        h_hot = hot_side_htc(T_hot_est, P_hot, X_air, D_ann, u_hot_ann)

        t_over_k = np.array([t / k for t, k in wall_layers])
        temps, q_wall = cb.wall_temperature_profile(T_hot_est, T2, h_hot, h_cool, t_over_k)

        # Coolant exit temperature from energy balance
        T_cool_exit = T2 + q_wall * A_liner / (mdot_cool * cp_cool)

        # Mix coolant_exit + bypass -> oxidizer
        cool_s = cb.Stream()
        cool_s.T, cool_s.P, cool_s.X, cool_s.mdot = T_cool_exit, P2, X_cool, mdot_cool
        byp_s = cb.Stream()
        byp_s.T, byp_s.P, byp_s.X, byp_s.mdot = T2, P2, X_air, mdot_bypass
        oxidizer = cb.mix([cool_s, byp_s])

        # Find fuel flow to hit COT_target (closed energy loop)
        fuel_s = cb.Stream()
        fuel_s.T, fuel_s.P, fuel_s.X = 300.0, P2, X_ch4
        fuel_s = cb.set_fuel_stream_for_Tad(COT_target, fuel_s, oxidizer)

        # Verify COT
        mixed_in = cb.mix([cool_s, byp_s, fuel_s])
        burned = cb.complete_combustion(mixed_in.T, mixed_in.X, mixed_in.P)
        T_hot_est = burned.T

    T_hot_surf = temps[0]
    T_metal_hot = temps[1] if len(temps) > 2 else temps[0]

    return {
        "f_cool": f_cool,
        "u_cool": u_cool,
        "h_cool": h_cool,
        "h_hot": h_hot,
        "Re": r_ch.Re,
        "Nu": r_ch.Nu,
        "dP_cool": r_ch.dP,
        "dP_total_frac": (dP_burner_frac * P2 + r_ch.dP) / P2 * 100.0,
        "T_cool_exit": T_cool_exit,
        "T_hot": T_hot_est,
        "FAR": fuel_s.mdot / mdot_total,
        "q_wall": q_wall,
        "T_hot_surf": T_hot_surf,
        "T_metal_hot": T_metal_hot,
        "eta_thp": eta_thp,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:

    print("=" * 70)
    print("Combustor Liner Convective Cooling -- Bypass Sweep (COT = const)")
    print("=" * 70)

    # -----------------------------------------------------------------------
    # Fixed operating conditions
    # -----------------------------------------------------------------------
    PR = 20.0
    P2 = PR * 101325.0
    T2 = 667.0  # K  compressor exit (isentropic, PR=20)
    RH_in = 0.60
    COT_target = 1758.0  # K  fixed combustor outlet temperature
    dP_burner_frac = 0.02  # 2% burner dP/P
    D_ann = 0.30  # m  combustor annulus hydraulic diameter
    u_hot_ann = 25.0  # m/s

    X_cool = cb.species.humid_air(288.15, 101325.0, RH_in)
    X_air = cb.species.dry_air()
    X_ch4 = cb.species.empty()
    X_ch4[cb.species.indices["CH4"]] = 1.0

    mdot_total = 1.0  # kg/s reference air flow

    # -----------------------------------------------------------------------
    # Channel geometry
    # -----------------------------------------------------------------------
    D_ch = 0.008  # m
    L_ch = 0.35  # m
    N_ch = 8  # parallel channels per liner panel

    # -----------------------------------------------------------------------
    # Materials
    # -----------------------------------------------------------------------
    T_metal_limit = 1073.15  # K = 800 C
    k_inconel = cb.k_inconel718(900.0)
    k_ysz_fresh = cb.k_tbc_ysz(1200.0, hours=0.0, is_ebpvd=False)
    k_ysz_aged = cb.k_tbc_ysz(1200.0, hours=1000.0, is_ebpvd=False)
    k_ysz_ebpvd = cb.k_tbc_ysz(1200.0, hours=0.0, is_ebpvd=True)
    t_wall = 0.003
    t_tbc = 0.0003

    wall_configs = {
        "Inconel only": [(t_wall, k_inconel)],
        "Inconel + YSZ APS (fresh)": [(t_tbc, k_ysz_fresh), (t_wall, k_inconel)],
        "Inconel + YSZ APS (1000 h)": [(t_tbc, k_ysz_aged), (t_wall, k_inconel)],
        "Inconel + YSZ EB-PVD": [(t_tbc, k_ysz_ebpvd), (t_wall, k_inconel)],
    }

    print(f"\nMaterials: {cb.list_materials()}")
    print(f"Inconel 718 k = {k_inconel:.2f} W/(m*K) @ 900 K")
    print(
        f"YSZ APS fresh k = {k_ysz_fresh:.3f}, aged k = {k_ysz_aged:.3f}, EB-PVD k = {k_ysz_ebpvd:.3f} W/(m*K) @ 1200 K"
    )
    print(f"Wall: {t_wall * 1e3:.1f} mm Inconel + {t_tbc * 1e6:.0f} um TBC")
    print(
        f"\nCOT target = {COT_target:.0f} K  |  dP_burner/P = {dP_burner_frac * 100:.0f}%  |  Metal limit = {T_metal_limit - 273.15:.0f} deg C"
    )
    print(f"Channels: {N_ch} parallel, D_h={D_ch * 1e3:.0f} mm, L={L_ch * 1e3:.0f} mm")

    # -----------------------------------------------------------------------
    # Bypass sweep
    # -----------------------------------------------------------------------
    f_cool_vals = np.linspace(0.05, 0.35, 25)
    channel_names = ["Smooth", "Ribbed"]

    # Primary wall config for sweep plots
    primary_cfg = "Inconel only"
    primary_layers = wall_configs[primary_cfg]

    sweep: dict[str, list[dict]] = {n: [] for n in channel_names}

    for f_cool in f_cool_vals:
        for ch_name in channel_names:
            pt = solve_operating_point(
                f_cool,
                mdot_total,
                T2,
                P2,
                X_cool,
                X_air,
                X_ch4,
                COT_target,
                dP_burner_frac,
                D_ch,
                L_ch,
                N_ch,
                D_ann,
                u_hot_ann,
                ch_name,
                primary_layers,
            )
            sweep[ch_name].append(pt)

    # -----------------------------------------------------------------------
    # Print summary
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print(f"Bypass Sweep  ({primary_cfg})")
    print(f"{'=' * 70}")
    hdr = f"  {'f_cool':>6}  {'u[m/s]':>7}  {'h_cool':>7}  {'T_exit[K]':>9}  {'T_metal[C]':>11}  {'FAR':>7}  {'dP_tot/P%':>10}  {'eta_thp':>7}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    for ch_name in channel_names:
        print(f"\n  Channel: {ch_name}")
        for pt in sweep[ch_name][::4]:  # every 4th point
            ok = "" if pt["T_metal_hot"] < T_metal_limit else " !"
            print(
                f"  {pt['f_cool']:6.2f}  {pt['u_cool']:7.1f}  {pt['h_cool']:7.0f}"
                f"  {pt['T_cool_exit']:9.1f}  {pt['T_metal_hot'] - 273.15:11.1f}{ok}"
                f"  {pt['FAR']:7.5f}  {pt['dP_total_frac']:10.2f}  {pt['eta_thp']:7.3f}"
            )

    # -----------------------------------------------------------------------
    # Minimum f_cool per channel x wall config
    # -----------------------------------------------------------------------
    print(f"\n{'=' * 70}")
    print("Minimum bypass fraction to achieve T_metal < 800 C")
    print(f"{'=' * 70}")
    print(f"  {'Channel':10}  {'Wall config':35}  {'f_min':>6}  {'dP_tot/P%':>10}  {'FAR':>8}")
    print(f"  {'-' * 10}  {'-' * 35}  {'-' * 6}  {'-' * 10}  {'-' * 8}")

    design_points: dict[str, dict[str, dict]] = {n: {} for n in channel_names}

    for ch_name in channel_names:
        for cfg_name, layers in wall_configs.items():
            f_min_pt = None
            for f_cool in f_cool_vals:
                pt = solve_operating_point(
                    f_cool,
                    mdot_total,
                    T2,
                    P2,
                    X_cool,
                    X_air,
                    X_ch4,
                    COT_target,
                    dP_burner_frac,
                    D_ch,
                    L_ch,
                    N_ch,
                    D_ann,
                    u_hot_ann,
                    ch_name,
                    layers,
                )
                if pt["T_metal_hot"] < T_metal_limit:
                    f_min_pt = pt
                    break
            design_points[ch_name][cfg_name] = f_min_pt
            if f_min_pt:
                print(
                    f"  {ch_name:10}  {cfg_name:35}  {f_min_pt['f_cool']:6.2f}"
                    f"  {f_min_pt['dP_total_frac']:10.2f}  {f_min_pt['FAR']:8.5f}"
                )
            else:
                print(f"  {ch_name:10}  {cfg_name:35}  {'NOT MET':>6}")

    # -----------------------------------------------------------------------
    # Plots
    # -----------------------------------------------------------------------
    colors = {"Smooth": "steelblue", "Ribbed": "darkorange"}

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        f"Combustor Liner Cooling -- Bypass Sweep  "
        f"(COT={COT_target:.0f} K, PR={PR:.0f}, dP_burner={dP_burner_frac * 100:.0f}%)\n"
        f"{N_ch} channels D_h={D_ch * 1e3:.0f} mm, L={L_ch * 1e3:.0f} mm  |  {primary_cfg}",
        fontsize=11,
    )

    # Plot 1: T_metal vs f_cool
    ax = axes[0, 0]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        T_vals = [pt["T_metal_hot"] - 273.15 for pt in sweep[ch_name]]
        ax.plot(f_vals, T_vals, color=colors[ch_name], lw=2, label=ch_name)
    ax.axhline(
        T_metal_limit - 273.15,
        color="red",
        ls="--",
        lw=1.5,
        label=f"Limit {T_metal_limit - 273.15:.0f} deg C",
    )
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Metal hot-face temperature [deg C]")
    ax.set_title("Metal Temperature vs Bypass Fraction")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 2: Total dP/P vs f_cool
    ax = axes[0, 1]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        dp_vals = [pt["dP_total_frac"] for pt in sweep[ch_name]]
        ax.plot(f_vals, dp_vals, color=colors[ch_name], lw=2, label=ch_name)
    ax.axhline(
        dP_burner_frac * 100,
        color="gray",
        ls=":",
        lw=1.5,
        label=f"Burner only {dP_burner_frac * 100:.0f}%",
    )
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Total dP/P [%]  (burner + coolant)")
    ax.set_title("Counterflow Pressure Budget")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: FAR vs f_cool (fuel penalty for wall cooling)
    ax = axes[1, 0]
    for ch_name in channel_names:
        f_vals = [pt["f_cool"] for pt in sweep[ch_name]]
        far_vals = [pt["FAR"] * 100 for pt in sweep[ch_name]]
        ax.plot(f_vals, far_vals, color=colors[ch_name], lw=2, label=ch_name)
    ax.set_xlabel("Bypass fraction f_cool [-]")
    ax.set_ylabel("Fuel/air ratio FAR [%]")
    ax.set_title("Fuel Flow vs Bypass Fraction\n(more cooling -> more fuel to maintain COT)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: T_metal at f_min per channel x wall config
    ax = axes[1, 1]
    cfg_names = list(wall_configs.keys())
    x = np.arange(len(cfg_names))
    width = 0.25
    for i, ch_name in enumerate(channel_names):
        T_vals = []
        for cfg_name in cfg_names:
            pt = design_points[ch_name].get(cfg_name)
            T_vals.append(pt["T_metal_hot"] - 273.15 if pt else float("nan"))
        ax.bar(x + i * width, T_vals, width, color=colors[ch_name], label=ch_name, alpha=0.85)
    ax.axhline(
        T_metal_limit - 273.15,
        color="red",
        ls="--",
        lw=1.5,
        label=f"Limit {T_metal_limit - 273.15:.0f} deg C",
    )
    ax.set_xticks(x + width)
    ax.set_xticklabels(
        [c.replace(" + ", "\n+ ").replace(" (", "\n(") for c in cfg_names],
        fontsize=7,
    )
    ax.set_ylabel("Metal temperature at f_min [deg C]")
    ax.set_title("Metal Temp at Minimum Bypass\n(per channel x wall config)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    plt.tight_layout()
    show_or_save(fig, "combustor_liner_cooling.png")


if __name__ == "__main__":
    main()
