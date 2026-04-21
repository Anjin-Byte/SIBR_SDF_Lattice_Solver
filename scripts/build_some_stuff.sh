#!/usr/bin/env bash
# Build classic-MC versions of the four lattice deliverables, then remesh
# each to slicer-ready size via the xtask QEM tool.
#
# Geometry: cylinder r=17 mm, h=50 mm, aligned along +z.
# Mesh settings: grid-ratio 3 (printer-resolution-matched), Taubin smooth 15 iter.
# Remesh target: 0.0002 (~0.012 mm absolute) — brushes the printer resolution
# limit without crossing it.
#
# Outputs:
#   deliverables/classic/*.stl          — raw MC output (large)
#   deliverables/remesh/classic/*.stl   — decimated, slicer-ready
#
# Usage (from workspace root):
#   ./scripts/build_some_stuff.sh
#
# Estimated wall time on this machine: ~3 minutes of tool execution + cargo
# compile time (cached after first run).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p deliverables/classic deliverables/remesh/classic

# Shared flags — kept unquoted at the call site so shell word-splits them.
CYL="--primitive cylinder --cylinder-start 0,0,0 --cylinder-end 0,0,50 --cylinder-radius 17"
COMMON="--grid-ratio 3 --extraction-method classic --smooth-iterations 15"

echo "=== [1/4] Classic MC: Kelvin L=3, r*=0.1 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology kelvin --cell-length 3 --strut-radius 0.3 \
  $COMMON -o deliverables/classic/kelvin_L3_rstar0.1.stl

echo "=== [2/4] Classic MC: Kelvin L=2, r*=0.08 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology kelvin --cell-length 2 --strut-radius 0.16 \
  $COMMON -o deliverables/classic/kelvin_L2_rstar0.08.stl

echo "=== [3/4] Classic MC: BccXy L=3, r*=0.1 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology bccxy --cell-length 3 --strut-radius 0.3 \
  $COMMON -o deliverables/classic/bccxy_L3_rstar0.1.stl

echo "=== [4/4] Classic MC: BccXy L=2, r*=0.08 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology bccxy --cell-length 2 --strut-radius 0.16 \
  $COMMON -o deliverables/classic/bccxy_L2_rstar0.08.stl

echo
echo "=== Remesh pass (target_error = 0.0002, ~0.012 mm absolute) ==="
for name in \
  kelvin_L3_rstar0.1 \
  kelvin_L2_rstar0.08 \
  bccxy_L3_rstar0.1 \
  bccxy_L2_rstar0.08; do
  echo "--- remesh: $name ---"
  cargo run --release -q -p xtask -- remesh \
    -i "deliverables/classic/$name.stl" \
    -o "deliverables/remesh/classic/$name.stl"
done

echo
echo "=== Final sizes ==="
ls -lh deliverables/classic/ deliverables/remesh/classic/
