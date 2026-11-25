# Thermo Data Generator

Tools to extract thermodynamic and transport data from external databases and
generate `thermo_transport_data.h` for the C++ library.

## Adding a Species

This is the primary use case. The complete workflow for adding one or more species:

> **The equilibrium does not pick up a new species.** `combustion_equilibrium`
> equilibrates reforming and the water-gas shift over the current table only; N2
> is inert. A new species (NO, NO2, OH, ...) needs a decision on its reactions
> first. `python/tests/test_equilibrium_scope.py` fails until that is made; see
> the Chemical Equilibrium section of `docs/API_CPP.md`.

### Step 1 — Obtain source data

Source files are **not included** in this repository (licensing). You need at
least one Cantera YAML mechanism covering the new species. If the species
appears in multiple mechanisms (common for transport data), you can supply all
of them and the extractor will merge them.

**Cantera YAML** (NASA-7 thermo + Lennard-Jones transport):
- Download from [Cantera](https://cantera.org/examples/input-files.html)
- Or convert Chemkin-format files with `ck2yaml`:

```bash
ck2yaml \
    --thermo NUIGMech1.1.THERM \
    --transport NUIGMech1.1.TRAN \
    --output NUIGMech1.1.yaml \
    --permissive
```

**NASA CEA** (NASA-9 thermo, up to 20 000 K — no transport):
- Use the [NASA CEA web interface](https://cearun.grc.nasa.gov/)
- Select species, run, save the text output

### Step 2 — Extract species data

The current species list is recorded in `.generator_state.json`. Read it to
build the `-s` argument: take all existing species and append the new one(s).

```bash
cd thermo_data_generator

# Single YAML source (transport data for all species in one mechanism)
uv run extract_species_data.py \
    --yaml JetSurf2.yaml \
    -s N2 O2 AR CO2 H2O CH4 C2H6 C3H8 IC4H10 NC5H12 NC6H14 NC7H16 CO H2 NEW_SPECIES \
    -o merged_species.json

# Multiple YAML sources — when transport data is split across mechanisms.
# First file wins on conflict; later files fill gaps.
uv run extract_species_data.py \
    --yaml JetSurf2.yaml secondary.yaml \
    -s N2 O2 AR CO2 H2O CH4 C2H6 C3H8 IC4H10 NC5H12 NC6H14 NC7H16 CO H2 NEW_SPECIES \
    -o merged_species.json

# Add NASA-9 thermo via CEA (higher accuracy, wider T range).
# YAML still supplies transport; CEA overlays thermo.
uv run extract_species_data.py \
    --yaml JetSurf2.yaml secondary.yaml \
    --cea cea_output.txt \
    -s N2 O2 AR CO2 H2O CH4 C2H6 C3H8 IC4H10 NC5H12 NC6H14 NC7H16 CO H2 NEW_SPECIES \
    -o merged_species.json
```

The output `merged_species.json` is gitignored (large derived file, regenerated
on demand).

### Step 3 — Add the common name (manual)

`include/common_names.h` maps formula → human-readable name. This is the
**only manual step** — common names are editorial choices that cannot be
derived automatically.

Edit **both** maps:

```cpp
// formula_to_name
{"NEW_SPECIES", "Human Readable Name"},

// name_to_formula
{"Human Readable Name", "NEW_SPECIES"},
```

Both maps must stay in sync. The pre-commit hook `check-common-names.sh` will
catch any missing entries.

### Step 4 — Regenerate the C++ header

Use the species list from `.generator_state.json` plus your new species:

```bash
uv run generate_thermo_data.py \
    --json merged_species.json \
    --species "N2,O2,AR,CO2,H2O,CH4,C2H6,C3H8,IC4H10,NC5H12,NC6H14,NC7H16,CO,H2,NEW_SPECIES" \
    --prefer-nasa9 \
    --check-names --strict \
    --output ../include/thermo_transport_data.h
```

`--check-names --strict` exits non-zero if any species is missing a common name
entry in `common_names.h` — catches Step 3 omissions before the build.

Then regenerate the C2 bridges of the NASA-9 joins (#549) from the new header:

```bash
uv run python thermo_data_generator/generate_c2_refit.py   # from the repo root
```

### Step 5 — Rebuild the C++ extension

```bash
cd ..
uv pip install -e .
```

### Step 6 — Verify

```bash
# Smoke-test Python bindings
uv run python -c "import combaero as cb; print(cb.species.names)"

# Full test suite
uv run pytest python/tests/ -v
```

The test suite and GUI adapt automatically — there are no hardcoded species
counts or lists in tests, Python bindings, or the GUI.

### What adapts automatically

| Layer | How it adapts |
|-------|--------------|
| C++ library | Loops use `species_names.size()` — no hardcoded counts |
| Python API (`cb.species`) | Reads count from C++ core at import time |
| Python unit tests | All use `cb.species.num_species` / `cb.num_species()` dynamically |
| Cantera validation tests | All use `num_species()` + `species_name(i)` from `_core` |
| GUI backend (`/metadata/species`) | Returns `cb.species.names` from C++ |
| GUI frontend | Fetches species list from API at startup |

### What requires a manual edit

| File | What to add | Why manual |
|------|-------------|-----------|
| `include/common_names.h` | Entry in both `formula_to_name` and `name_to_formula` | Common names are editorial — cannot be derived from formulas |

---

## Scripts

| Script | Purpose |
|--------|---------|
| `extract_species_data.py` | Extracts thermo + transport from YAML and/or CEA, merges multiple sources |
| `generate_thermo_data.py` | Generates `thermo_transport_data.h` from merged JSON |
| `generate_c2_refit.py` | Generates `include/thermo_c2_refit_data.h`, the C2 bridges of the NASA-9 joins (#549), from the coefficients already in `thermo_transport_data.h`. Rerun after regenerating that header. |

### Legacy (deprecated)

| Script | Purpose |
|--------|---------|
| `create_cpp_header.py` | Old CSV-to-header converter; superseded by `generate_thermo_data.py` |

---

## Data Sources

### Cantera YAML Mechanisms

Provide NASA-7 polynomial thermo and Lennard-Jones transport parameters.
Recommended sources (covering different species families):

| Mechanism | Coverage |
|-----------|---------|
| GRI-Mech 3.0 | Natural gas, C1–C2 hydrocarbons |
| JetSurF 2.0 | Jet fuel surrogates, C1–C16 |
| San Diego mechanism | Hydrocarbon combustion |
| NUIGMech 1.1 | Comprehensive: C0–C7 + nitrogen chemistry (NH3, NOx) |

When a species' transport data appears in a secondary mechanism but not your
primary one, pass both to `extract_species_data.py --yaml primary.yaml secondary.yaml`.
The first file takes priority on conflicts; later files fill gaps.

**Flagged transport yields (#485).** A mechanism can mark a species' transport
as an unreferenced estimate. NUIGMech1.1 does so for every species here:
`\AUTHOR: WARNING !\REF: WARNING !\COMMENT: theoret trans`. Such a set is
replaced by the first later source with unflagged parameters, for TRANSPORT
ONLY; thermo keeps first-wins. Each entry records `transport_source` and
`transport_flagged` in `merged_species.json`. Without this, NUIGMech's O2
(eps 676 K against 107.4 K in GRI30, JetSurf2, Aramco2 and San Diego) made
air's viscosity 11.8% low. The current header comes from
`--yaml NUIGMech1.1.yaml JetSurf2.yaml --cea NASA9_subset.txt`: thermo from
NASA-9, transport from JetSurf2, and NUIGMech for NH3 alone.

Check new transport data against literature, not only against another
mechanism (`python/tests/test_transport_literature.py`).

### Collision integrals (`collision_integrals.py`)

Generates `include/collision_integral_data.h`: Omega*(2,2) on a 41 T* x 8
delta* grid and its slope d ln Omega*/d ln T* at every node, for the C1 cubic
Hermite interpolation in `src/transport.cpp`. Run from the repository root
with the main environment (it needs numpy/scipy; about 3 minutes on 8 cores):

```bash
uv run python thermo_data_generator/collision_integrals.py
```

The delta* = 0 column is computed, not transcribed: classical deflection
angle (outermost turning point), cross-section split at the orbiting impact
parameter, thermal average and its T* derivative. It refuses to write if it
disagrees with Kim & Monroe (2014) by more than 1e-4 anywhere in 0.3-400.
The polar columns are Monchick & Mason (1961) scaled by the delta* = 0
column's exact/tabulated ratio, with spline-estimated slopes.

**Note — PyYAML boolean parsing**: species names like `NO` parse as `False`
in YAML 1.1. The extractor handles this automatically.

### NASA CEA Database

Provides NASA-9 polynomials with an extended temperature range (200–20 000 K).
No transport data — transport must come from a YAML mechanism.

Use the [NASA CEA web interface](https://cearun.grc.nasa.gov/), select species,
run, and save the text output. Pass it via `--cea` to overlay NASA-9 thermo on
top of the YAML-sourced transport data.

### Chemkin-format files

Convert to Cantera YAML first with `ck2yaml` (installed with Cantera):

```bash
ck2yaml --thermo file.THERM --transport file.TRAN --output out.yaml --permissive
```

---

## NASA Polynomial Formats

### NASA-7 (from YAML)

Two temperature ranges, 7 coefficients each:

```
Cp/R = a1 + a2*T + a3*T² + a4*T³ + a5*T⁴
H/RT = a1 + a2*T/2 + a3*T²/3 + a4*T³/4 + a5*T⁴/5 + a6/T
S/R  = a1*ln(T) + a2*T + a3*T²/2 + a4*T³/3 + a5*T⁴/4 + a7
```

### NASA-9 (from CEA)

Up to three temperature intervals, 10 coefficients each (wider T range):

```
Cp/R = a1/T² + a2/T + a3 + a4*T + a5*T² + a6*T³ + a7*T⁴
```

`generate_thermo_data.py --prefer-nasa9` (default) uses NASA-9 where available
and falls back to NASA-7 otherwise.

---

## Output

### JSON (`merged_species.json`)

Intermediate file, gitignored. Contains merged thermo and transport data from
all sources. Passed to `generate_thermo_data.py`.

### Generator state (`.generator_state.json`)

Tracked in git. Records the species list, JSON source path, and NASA format
used in the last successful generation — the canonical reference for what is
currently in `thermo_transport_data.h`.

### C++ Header (`include/thermo_transport_data.h`)

Auto-generated. Do not edit manually. Contains:

- `species_names` — `std::vector<std::string>`
- `species_index` — `std::unordered_map<std::string, int>` for O(1) lookup
- `molar_masses` — `std::vector<double>` [g/mol]
- `nasa_coeffs` — NASA-7 or NASA-9 structs depending on `--prefer-nasa9`
- `transport_props` — Lennard-Jones / Stockmayer transport parameters per species:
  geometry, well_depth [K], diameter [Ang], polarizability [Ang³],
  dipole_moment [Debye] (0.0 for non-polar), z_rot [-]
- `molecular_structures` — atomic composition (C, H, O, N counts)

---

## Testing

```bash
uv run pytest
```

---

## Licensing Note

Source data files (YAML mechanisms, CEA output) are not included in this
repository. Many thermodynamic databases have specific licensing terms:

- Some mechanisms incorporate Burcat's database (restricted use)
- NASA CEA data is public domain
- Cantera mechanism files have varying licenses

Always verify the license of your source data before use.
