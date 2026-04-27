#!/usr/bin/env bash
# Build a "wafer" lattice sample for visual inspection in Blender +
# automated manifoldness analysis. **Pure MC33 output** — no Taubin
# smoothing, no QEM remesh, no other post-processing. This isolates the
# extraction stage from any pipeline transforms so the geometry you load
# in Blender is exactly what MC33 emitted (and then welded).
#
# Geometry: cylinder r=12 mm, h=16 mm with Kelvin L=4 cells. Several
# cells deep so cross-voxel symmetry features show up at scale.
#
# Resolution: `--grid-ratio 6` — twice the production-default density
# (vs. `build_some_stuff.sh`'s grid-ratio 3). Voxel size = strut_radius / 6
# = 0.3 / 6 = 0.05 mm. Roughly 8× the triangle count of grid-ratio=3 at
# the same geometry. Tune higher (8, 10, ...) for finer detail or lower
# (4, 5) if file size becomes unwieldy.
#
# Outputs in deliverables/wafer/:
#   kelvin_mc33_smooth0.stl   — MC33, no smoothing (raw extraction)
#
# (Comparison builds — classic MC and smoothed MC33 — are commented out
# below; uncomment if you want side-by-side artifacts.)
#
# Each STL is followed by a manifoldness summary via scripts/check_manifold.py.
#
# Usage (from workspace root):
#   ./scripts/build_wafer_diag.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

OUT=deliverables/wafer
mkdir -p "$OUT"

CYL="--primitive cylinder --cylinder-start 0,0,0 --cylinder-end 0,0,16 --cylinder-radius 12"
CELL="--cell-topology kelvin --cell-length 4 --strut-radius 0.3"
# grid-ratio 6 doubles the production density. Bump higher for finer
# detail; the analyzer scales fine well into the millions of triangles.
COMMON="--grid-ratio 6"

build() {
  local label="$1"; shift
  local stl="$OUT/$label.stl"
  echo "=== build: $label ==="
  RUST_LOG=warn cargo run --release -q -p sibr-lattice -- \
    $CYL $CELL $COMMON "$@" -o "$stl"
  python3 "$REPO_ROOT/scripts/check_manifold.py" "$stl"
  echo
}

build kelvin_mc33_smooth0   --extraction-method mc33    --smooth-iterations 0
#build kelvin_mc33_smooth15  --extraction-method mc33    --smooth-iterations 15
#build kelvin_classic_smooth0 --extraction-method classic --smooth-iterations 0

echo "=== file sizes ==="
ls -lh "$OUT/"
echo
echo "Tip: open the STLs in Blender, then Edit Mode → Select → All by Trait →"
echo "Non-Manifold. Note that Blender's selector includes BOUNDARY edges by"
echo "default — boundary edges are open edges (degree 1), not topology bugs."
echo "Use the analyzer above for the unambiguous breakdown."
