# Development Workflow

This document outlines the development practices, testing procedures, and version synchronization rules for the CombAero project.

## Environment Setup

The root `.venv` is the single source of truth for all Python-based development tasks.

```bash
./scripts/bootstrap.sh
```

## Testing

### C++ Tests
C++ tests use GoogleTest and are integrated with CTest.

```bash
cd build
ctest --output-on-failure
```

### Python Tests

```bash
uv run pytest python/tests/
```

### Cantera Validation

```bash
CANTERA_DATA="$PWD/cantera_validation_tests" uv run --project cantera_validation_tests pytest -v
```

## Version Synchronization

When updating tool versions or Python targets, ensure the following files remain synchronized:

1. `.python-version`
2. `ruff.toml` (`target-version`)
3. `.github/workflows/ci.yml`
4. `.github/workflows/validation-tests.yml`
5. `.pre-commit-config.yaml`
6. `Makefile` and the `.sh` scripts in `scripts/`

When updating frontend tool versions, keep these synchronized: `gui/frontend/package.json`, `biome.json`
(`$schema` version), `.pre-commit-config.yaml` (biome rev and the `@biomejs/biome` pin) and `ci.yml`
(`setup-biome` version). `scripts/check_version_sync.py` enforces the biome set (CI and pre-commit).

Releases, Dependabot rounds and dependency floors: see [RELEASING.md](RELEASING.md).

## Pre-commit Hooks

The project uses `pre-commit` to ensure code style and basic checks pass before every commit.

```bash
# Install hooks (one-time)
pre-commit install

# Run manually on all files
pre-commit run --all-files
```

## Documentation

### Units Reference
`docs/UNITS.md` is auto-generated from `include/units_data.h`. To update it:

```bash
uv run scripts/generate_units_md.py
```
