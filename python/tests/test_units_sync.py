import inspect

import combaero as cb

# These are pure Python utilities, constants, or modules that
# do not have physical units.
IGNORE_LIST = {
    # Modules
    "compressible",
    "incompressible",
    "heat_transfer",
    "network",
    "species",
    "geometry",
    "vortex",
    "_solver_tools",
    # Utilities
    "get_warning_handler",
    "set_warning_handler",
    "is_well_behaved",
    "suppress_warnings",
    # Unit query API itself
    "has_units",
    "get_units",
    "output_units",
    "input_units",
    "all_units",
    "list_functions_with_units",
    # String manipulation / helpers
    "common_name",
    "formula_to_name",
    "name_to_formula",
    "formula",
    "list_materials",
    "num_species",
    "species_name",
    "species_index_from_name",
    "species_molar_mass",
    "species_molar_mass_from_name",
    # Mix is a special case function combining streams
    "mix",
    # Enums
    "CombustionMethod",
    "BoundaryCondition",
    "CorrelationValidity",
    "CorrelationResult",
    "ChamberResult",
    "MixtureState",
    "Stream",
    "AreaChangeResult",
    "AreaChangeLoss",
    "AreaChangeFlow",
    "AreaChangeElementFlow",
    "CompleteCombustionJacobian",
    "EquilibriumJacobian",
    "BordaCarnotStep",
    "TeeJunctionResult",
    "BranchInput",
    "CompressibleTeeResult",
    "compressible_branching_tee_rj",
    "compressible_merging_tee_rj",
    "tee_K5",
    "tee_K6",
    "tee_K11",
    "tee_K12",
    "tee_blend_weight",
    "tee_check_inputs",
    "merging_tee_K_straight",
    "merging_tee_K_branch",
    "branching_tee_K_straight",
    "branching_tee_K_branch",
    "State::DP",
    "State::HP",
    "State::PV",
    "State::SH",
    "State::SP",
    "State::SV",
    "State::UP",
    "State::UV",
    "State::VH",
    "State::set_P",
    "State::set_T",
    # Result aggregates and the unit registry itself. These carry no single
    # unit; their MEMBERS are what units_data.h describes, and pybind11
    # result structs expose those as plain attributes on an instance rather
    # than on the class, so dir(cls) does not reach them.
    "ChannelSolverResult",
    "IncompressibleFlowSolution",
    "MassStream",
    "EntropyGenerationResult",
    "MomentumChamberResult",
    "OrificeResult",
    "WallCouplingResult",
    "Registry",
}

# Special magic methods created by pybind11 that aren't public calculation points
IGNORE_METHODS = {
    "__init__",
    "__repr__",
    "__str__",
    "__copy__",
    "__deepcopy__",
    "__doc__",
    "__module__",
    "__class__",
}


def public_surface() -> list[str]:
    """Every name `import combaero` actually exposes.

    Deliberately NOT `cb.__all__`. That list is hand-maintained, so a symbol
    added without touching it escaped this check entirely -- 54 of them had,
    14 with no units entry at all. The surface a user can reach is what has to
    be covered, and it is the union: `dir()` for what is reachable, `__all__`
    for anything declared but lazily bound.
    """
    reachable = {n for n in dir(cb) if not n.startswith("_")}
    return sorted(reachable | set(cb.__all__))


def test_all_covers_the_public_surface() -> None:
    """`__all__` must not drift below what the module exposes.

    Keeps the two definitions of "public" from diverging again. A name that is
    reachable but undeclared is either API that belongs in `__all__` or an
    accidental import leak that should be underscored at its import site --
    `contextmanager`, `import_module`, `Generator` and `PackageNotFoundError`
    were all the latter.
    """
    reachable = {n for n in dir(cb) if not n.startswith("_")}
    undeclared = sorted(reachable - set(cb.__all__))
    assert not undeclared, (
        f"{len(undeclared)} names are reachable on `combaero` but absent from "
        "`__all__`. Add them if they are API; underscore the import if they "
        f"leaked: {undeclared}"
    )


def test_no_phantom_overloads() -> None:
    """A function must not advertise two identical signatures.

    pybind11 tries overloads in registration order, so a second `m.def` with
    the SAME signature can never be selected -- but it is still documented.
    `help()` prints "Overloaded function" and lists it as "2.", so a reader
    is invited to work out a difference that does not exist.

    30 functions carried one, `channel_mdot` two. Most arose from a
    correlation being bound in a general section and again in a
    domain-specific one; nothing failed, which is why they accumulated.

    Distinct signatures are left alone -- `orifice_flow_thermo(.., Cd)` and
    `orifice_flow_thermo(.., cd_fn)` are a real, useful overload pair.
    """
    import re

    offenders = {}
    for name in public_surface():
        doc = getattr(getattr(cb, name, None), "__doc__", None) or ""
        if "Overloaded function" not in doc:
            continue
        sigs = [
            m.strip()
            for m in re.findall(r"^\s*\d+\.\s+(.*)$", doc, re.M)
            if m.strip().startswith(name + "(")
        ]
        distinct = list(dict.fromkeys(sigs))
        if len(sigs) > len(distinct):
            offenders[name] = f"{len(sigs)} listed, {len(distinct)} distinct"

    assert not offenders, (
        "these functions register the same signature more than once in "
        "_core.cpp; every registration after the first is unreachable but "
        f"still shows up in help(): {offenders}"
    )


def test_api_unit_sync():
    """
    Dynamically discover all exposed functions and classes from combaero,
    and assert that they have corresponding unit metadata defined in units_data.h.
    This prevents the documentation and API from drifting out of sync.
    """
    missing = []

    for name in public_surface():
        if name in IGNORE_LIST:
            continue

        obj = getattr(cb, name, None)
        if obj is None:
            continue

        # Skip strings like __version__
        if isinstance(obj, str):
            continue

        if inspect.isclass(obj):
            # For classes, check their public members
            for attr in dir(obj):
                if attr.startswith("_"):
                    continue
                if attr in IGNORE_METHODS:
                    continue

                method_name = f"{name}::{attr}"
                if method_name in IGNORE_LIST:
                    continue
                if not cb.has_units(method_name):
                    missing.append(method_name)
        else:
            # Standalone function or property exported directly
            if not cb.has_units(name):
                missing.append(name)

    # Use standard assert for pytest
    assert len(missing) == 0, (
        f"The following {len(missing)} API members are missing from units_data.h:\n"
        + "\n".join(missing)
    )


def test_docstring_energy_basis_matches_units_data() -> None:
    """A docstring that states an energy per kg must not describe an output
    that units_data.h gives per mol, or the reverse (#1279: dh_dT, ds_dT and
    dcp_dT were documented as mass-specific, and return molar values). Only a
    docstring naming one basis is checked; one naming both may mix inputs and
    outputs."""
    import re

    kg = re.compile(r"J/\(?kg")
    mol = re.compile(r"J/\(?mol")
    mismatched = []
    for name in cb.__all__:
        if not cb.has_units(name):
            continue
        doc = getattr(getattr(cb, name, None), "__doc__", "") or ""
        out = cb.output_units(name) or ""
        doc_kg, doc_mol = bool(kg.search(doc)), bool(mol.search(doc))
        out_kg, out_mol = bool(kg.search(out)), bool(mol.search(out))
        if (doc_kg and not doc_mol and out_mol and not out_kg) or (
            doc_mol and not doc_kg and out_kg and not out_mol
        ):
            mismatched.append((name, out))
    assert not mismatched, mismatched
