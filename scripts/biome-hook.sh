#!/bin/bash
set -euo pipefail

# Pre-commit entry for Biome: runs the frontend's installed Biome, so the
# version is the one gui/frontend/package.json and pnpm-lock.yaml pin (the
# only files Dependabot bumps). scripts/check_version_sync.py guards that no
# second pin comes back.

cd "$(dirname "$0")/.."

BIOME="gui/frontend/node_modules/.bin/biome"
if [ ! -x "$BIOME" ]; then
    echo "ERROR: $BIOME not found. Run 'pnpm install' in gui/frontend first."
    exit 1
fi

exec "$BIOME" "$@"
