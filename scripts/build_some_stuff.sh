#!/usr/bin/env bash
# Build MC33 versions of the four lattice deliverables, then remesh each
# to slicer-ready size via the xtask QEM tool.
#
# Geometry: cylinder r=17 mm, h=50 mm, aligned along +z.
# Mesh settings: grid-ratio 3 (printer-resolution-matched), Taubin smooth 15 iter.
# Extraction: MC33 (Chernyaev 1995 / Lewiner 2003) — topologically correct
# at every non-degenerate configuration.
#
# Smoothing:
# - `--boundary-smoothness 0.001` — fillets the lattice-cylinder
#   intersection edge. Geometrically invisible at this scale (1 µm)
#   but eliminates the sub-voxel "noise islands" Marching Cubes can
#   produce where struts cross the curved cylinder wall at grazing
#   angles. See scripts/build_wafer_diag.sh for the diagnostic that
#   identified this and the SmoothIntersection combinator that
#   implements it.
# - `--strut-smoothness 0.3 * strut_radius` — fillets strut-to-strut
#   joints. Cosmetically rounds joint corners and can mitigate
#   stress concentrations in printed parts. Visually subtle at this
#   ratio.
#
# Remesh target: 0.0002 (~0.012 mm absolute) — brushes the printer resolution
# limit without crossing it.
#
# Outputs:
#   deliverables/mc33/*.stl          — raw MC33 output (large)
#   deliverables/remesh/mc33/*.stl   — decimated, slicer-ready
#
# Usage (from workspace root):
#   ./scripts/build_some_stuff.sh
#
# Estimated wall time on this machine: ~3 minutes of tool execution + cargo
# compile time (cached after first run).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

mkdir -p deliverables/mc33 deliverables/remesh/mc33

# Shared flags — kept unquoted at the call site so shell word-splits them.
# `--extraction-method mc33` is the CLI default as of the full MC33 port,
# but we pass it explicitly here so this script's behavior doesn't drift
# if the default ever changes again.
CYL="--primitive cylinder --cylinder-start 0,0,0 --cylinder-end 0,0,10 --cylinder-radius 17"
COMMON="--grid-ratio 6 --extraction-method mc33 --smooth-iterations 0 --boundary-smoothness 0.001"

echo "=== [1/4] MC33: Kelvin L=3, r*=0.1 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology kelvin --cell-length 3 --strut-radius 0.3 \
  --strut-smoothness 0.09 \
  $COMMON -o deliverables/mc33/kelvin_L3_rstar0.1.stl

echo "=== [2/4] MC33: Kelvin L=2, r*=0.08 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology kelvin --cell-length 2 --strut-radius 0.16 \
  --strut-smoothness 0.048 \
  $COMMON -o deliverables/mc33/kelvin_L2_rstar0.08.stl

echo "=== [3/4] MC33: BccXy L=3, r*=0.1 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology bccxy --cell-length 3 --strut-radius 0.3 \
  --strut-smoothness 0.09 \
  $COMMON -o deliverables/mc33/bccxy_L3_rstar0.1.stl

echo "=== [4/4] MC33: BccXy L=2, r*=0.08 ==="
cargo run --release -q -p sibr-lattice -- $CYL \
  --cell-topology bccxy --cell-length 2 --strut-radius 0.16 \
  --strut-smoothness 0.048 \
  $COMMON -o deliverables/mc33/bccxy_L2_rstar0.08.stl

#echo
#echo "=== Remesh pass (target_error = 0.0002, ~0.012 mm absolute) ==="
#for name in \
  #kelvin_L3_rstar0.1 \
  #kelvin_L2_rstar0.08 \
  #bccxy_L3_rstar0.1 \
  #bccxy_L2_rstar0.08; do
  #echo "--- remesh: $name ---"
  #cargo run --release -q -p xtask -- remesh \
    #-i "deliverables/mc33/$name.stl" \
    #-o "deliverables/remesh/mc33/$name.stl"
#done

#echo
#echo "=== Final sizes ==="
#ls -lh deliverables/mc33/ deliverables/remesh/mc33/
