"""The combustion equilibrium is partial by design: reforming and water-gas
shift over a fixed species table, N2 inert (no NOx), no dissociation.

A species added to the thermo table is not equilibrated by any reaction: it
would pass through at its fed amount with a zero Jacobian, and nothing else
would fail. This test pins the table so that an addition (NO, NO2, OH, ...)
stops here until its equilibrium treatment is decided. Scope: the Chemical
Equilibrium section of docs/API_CPP.md.
"""

import combaero as cb

# The species the reforming + water-gas-shift equilibrium was written for.
EQUILIBRIUM_SPECIES = (
    "N2",
    "O2",
    "AR",
    "CO2",
    "H2O",
    "CH4",
    "C2H6",
    "C3H8",
    "IC4H10",
    "NC5H12",
    "NC6H14",
    "NC7H16",
    "CO",
    "H2",
    "NH3",
)


def _species() -> tuple[str, ...]:
    return tuple(cb.species_name(i) for i in range(cb.num_species()))


def test_species_table_matches_equilibrium_scope() -> None:
    added = sorted(set(_species()) - set(EQUILIBRIUM_SPECIES))
    removed = sorted(set(EQUILIBRIUM_SPECIES) - set(_species()))
    assert not added and not removed, (
        f"species table changed (added {added}, removed {removed}). The "
        "combustion equilibrium only equilibrates reforming and the water-gas "
        "shift, with N2 inert: a new species is not equilibrated and passes "
        "through at its fed amount. Decide its reactions (or that it is inert), "
        "update docs/API_CPP.md (Chemical Equilibrium) and then this list."
    )
