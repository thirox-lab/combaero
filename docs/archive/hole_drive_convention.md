# Hole drive convention: static vs total supply reference (2026-11-12)

This is the design record behind `HoleCd(supply_drive=...)` (#1087, #1122) and the rule that every hole runs **from the upstream total to the downstream static**.

It is our own analysis of the held hole-Cd data. Every number is cross-source, nothing is fitted, and no figure is reproduced. The scripts and point tables were local (`tmp/cd_audit/drive_convention/`) and are not part of the repository.

## The rule (owner, 2026-11-12)

> "orifices / losses are Pt to static. That way they own the energy balance."

A discharge hole drives from the upstream TOTAL pressure to the downstream STATIC one. How much of a supply's dynamic head the hole recovers is its Cd model's business, not the drive's.

In a network this matters only where a hole's supply node does not carry the supply's total. That is a duct-fed `EffusionPlateElement`, such as an `EffusionLiner` panel, whose station is a `PlenumNode` holding the duct static. There the drive is rebuilt from the closure: `hole_supply_drive`, the station static plus the momentum-chamber head of the mean segment flow.

## The conversion, exact under each paper's own definitions

A hole has one measured flow m and two ideal flows, so the two Cds are related exactly:

    Cd_T = m / m_id(Pt1, Tt1 -> p2)        Cd_S = m / m_id(p1, T1 -> p2)
    Cd_T / Cd_S = m_id(p1, T1 -> p2) / m_id(Pt1, Tt1 -> p2)

Here p1, T1 is the supply station's static state, Pt1, Tt1 its stagnation state, and p2 the exit static.

- **Incompressible:** the ratio is `sqrt(1 - X)`, with `X = (Pt1 - p1)/(Pt1 - p2) = 1/VHR = u^2/(1 + u^2)` and `u = U1/Vi`.
- **Compressible:** the ratio of the isentropic fluxes, which depends on the supply Mach as well. On Gritsch's points (Ma_c to 0.6) it differs from `sqrt(1 - X)` by 0 to -4.7 %.
- **The static temperature matters.** Taking the static with the TOTAL temperature (the pre-0.12 station state) instead of the static temperature moves the ideal flow by 0 to +3.5 % on Gritsch's points.
- **Where p1 <= p2, the static Cd is undefined, although the hole flows.** This applies to 2 cylindrical and 14 shaped Gritsch points.

combaero implements the conversion as `hole_static_to_total_cd_ratio`. It takes both ideal flows from the network's own compressible nozzle kernel, so a static-referenced Cd on the total drive gives that model's own flow to round-off. An element in the incompressible regime (the solver's incompressible warm-start proxy included) takes the incompressible form `sqrt(1 - X)` instead (`incompressible=True`), so the proxy's forward kernel and its conversion are one model (review of #1122, B2).

## What each held source references

| source (held) | Cd from -> to | crossflow parameter |
|---|---|---|
| McGreehan & Schotsch 1988 | supply STATIC (pressure and temperature) -> exit static; they convert Rohde onto it | U1/Vi (Eq. 17) |
| Rohde, Richards & Metger 1969 | duct TOTAL -> jet static | VHR |
| Kaddah 1964 | approach TOTAL -> jet static, incompressible | K = VHR |
| Hay, Lampard & Benmansour 1983 | coolant stagnation -> mainstream static | Mc |
| Gritsch, Schulz & Wittig 1998 | coolant total -> mainstream static | Ma_c |
| Thole et al. 1997 | stagnation at hole entrance -> exit static | Ma_c |
| Adkins & Gueroui 1986 | annulus total -> core static | Phi = 1 - 1/VHR |

Six of the seven held sources are total-referenced; McGreehan & Schotsch alone reference the static. Rohde (p. 9) states the trade-off himself: the static suits a normal hole, because there is no head recovery, and the total suits an inclined one, which recovers some of it.

## Which convention collapses the data

- **Between labs, the convention cannot make sets agree better.** At matched X every set shifts by the same `0.5 ln(1 - X)`. Over 606 cross-lab point pairs (|dX| < 0.02), the RMS of ln(R_a/R_b) x 100 is 23.97 on the total basis and 23.31 on the static basis.
- **Within a set it does differ.** The table gives the RMS of ln(Cd/Cd0) x 100; smaller means that convention alone absorbs the supply crossflow.

| holes | static | total |
|---|---|---|
| normal: Kaddah, Rohde, MS Fig. 6, Hay (4 labs, 11 sets) | **15.8** | 32.5 |
| 30 deg leaning across the supply flow (Hay, with mainstream) | **5.9** | 16.2 |
| 30 deg leaning with the supply flow (Gritsch, 64 points) | 57.9 | **14.5** |

## The finding: a recovered-head fraction f

Define f through `(Cd_T/Cd0)^2 = (1 - X) + f X`: f = 0 means the hole sees the static, and f = 1 means it sees the total. f does not depend on the convention the data are stored in. For X >= 0.1:
- normal holes: median -0.18, i.e. about 0;
- 30 deg holes leaning across the flow: +0.16;
- 30 deg co-aligned holes (Gritsch): +1.25.

**Choosing a drive is therefore choosing f.** Under the rule, the drive is total and f belongs to the Cd model:
- McGreehan-Schotsch's Eq. 17 Cd is static-referenced. It is converted exactly, giving f = 0, which the normal-hole and cross-flow-leaning data support.
- A base without an inlet term (`inlet_crossflow="none"`) gave f = 1 as first shipped in #1122. That is the co-aligned behaviour, and on liners it is a large change.
- **Decision (owner, 2026-11-13): f is one field, `HoleCd.supply_recovery`, default 0 at every angle and for every base** (C1). Every base Cd is referenced to P1 + f q1 and converted exactly onto the total drive by `hole_supply_recovery_cd_ratio`, Tf being the isentropic temperature at P1 + f q1 on the isentrope through the static and stagnation states (both ends exact). The basis is the normal-hole row above (15.8 vs 32.5). `supply_recovery=1.0` is the co-aligned behaviour and replaces the unreleased `supply_drive="total"`; `supply_drive="static"` remains only as the exact v0.11.1 reproduction. f applies only to holes fed by a duct crossflow; every other hole refuses it.
- A lean- and direction-dependent f is #1124. It needs a second co-aligned lab (Rohde 1969 Fig. 11, a 45 deg hole leaning with the duct flow, held) and orientation data: Gritsch, Schulz & Wittig 2001, J. Turbomach. 123(4), 781-787, doi:10.1115/1.1397306.

## Consequences recorded in #1122

- **The old normal-hole comparison measured a mismatch, not the drive.** Applying the static Eq. 17 Cd unconverted to a total drive gave a normal-hole RMS of 17.5 -> 29.6 %. Converted, the result equals the static one.
- **Reversed duct flow is now demoted.** Booking a duct-fed hole's stream at the duct total (the rule) exposes a station model that manufactures total pressure under reversed duct flow (#1125). Such a liner is demoted by the second-law check until #1125 is fixed.
