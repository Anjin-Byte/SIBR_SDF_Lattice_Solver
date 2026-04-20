---
title: Meshing Complexity Analysis
tags: [domain, meshing, complexity, performance]
---

# Meshing Complexity Analysis

The memory-and-time budget for turning an SDF-represented lattice into a watertight triangle mesh. This is the "how much will it cost" reference that complements [[Lattice Generation Pipeline]] ("what are the stages").

## Correctness model

The meshing task: given `f: ℝ³ → ℝ` representing the lattice body (strut union, trimmed to primitive, with conformal skin), produce triangle mesh `M = (V, F)` that is:

- **Manifold** — every edge shared by exactly 2 faces
- **Closed** — no boundary edges
- **Consistently oriented** — face normals outward
- **Geometrically faithful** — Hausdorff distance from `M` to `{f=0}` bounded by the extraction resolution

"Watertight" = first three. "Performant" = tractable on a single workstation GPU for parts up to a few cm per side.

## Load-bearing assumptions

1. **Lattice SDF is expressible via [[SDF Primitives Catalog|repetition primitives]]** (`opRepetition` / `opLimitedRepetition`). Per-query cost is `O(k)` where `k` = struts per cell, not `O(N_struts)`. Without this, the approach dies.
2. **Strut SDFs are exact**, not bound. Bound contamination (from smooth combinators, twists, bends) still works for MC extraction but breaks sphere-tracing guarantees. See [[SDF Primitives Catalog]] "Bound contamination".
3. **Grid resolution `h ≲ r/3`** — finer is decorative, coarser is non-manifold.
4. **Primitive volume `V` and cell length `L_c` dominate**. Topology affects `k` (3–12) and surface area ratio (constant factor).
5. **GPU memory is the binding constraint.** Consumer cards: 8–24 GB. Swapping to host is catastrophic.

## Variables

| Symbol | Meaning |
|---|---|
| `V` | Primitive bounding volume |
| `L_c` | Unit cell length |
| `r` | Strut radius |
| `h` | Grid resolution (voxel edge length) |
| `k` | Struts per unit cell (topology-dependent, typically 3–12) |
| `G = V/h³` | Dense voxel count |
| `S ≈ 2π r k V / L_c²` | Total strut surface area |
| `G_s = S/h²` | Surface voxel count (sparse active) |
| `M ≈ 2 G_s` | Output triangle count (MC worst case) |

## Four candidate approaches

### Method A — Dense grid + Marching Cubes (MC33)

Evaluate `f` on a regular grid, run MC on every voxel.

- **Field memory:** `O(G) × 4 B`
- **Per-query work:** `O(k)`
- **Total eval time:** `O(k G)`
- **Output mesh:** `O(M)` triangles
- **Watertight + 2-manifold:** ✓ by construction (MC33 tables)
- **Parallelism:** embarrassing (per-voxel)

### Method B — Sparse / bricked grid + MC33

Activate only voxels with `|f| < ε · h`. Store in bricked structure (OpenVDB-style or wgpu-native).

- **Field memory:** `O(G_s) × 8 B`
- **Per-query work:** `O(k)`
- **Total eval time:** `O(k G_s × depth)`
- **Output mesh:** `O(M)` (same as dense)
- **Watertight + 2-manifold:** ✓ with brick-seam care (shared-sample aprons at brick boundaries)
- **Parallelism:** per-brick + per-voxel

### Method C — Octree + Dual Contouring

Adaptive octree refinement, dual contouring at leaves.

- **Field memory:** `O(G_s)`
- **Per-query work:** `O(k)` + octree navigation
- **Output mesh:** `O(M)` with fewer triangles per feature (~½ of MC)
- **Watertight:** ✓
- **2-manifold:** ✗ with standard DC; ✓ with Manifold DC (see [[Isosurface Extraction Methods]])
- **Parallelism:** harder (octree traversal on GPU is involved)

### Method D — Direct per-strut mesh + Boolean union

Tessellate each strut analytically (parameterized capsule), union the meshes.

- **Mesh memory:** `O(N_struts × T_strut)`
- **Per-strut work:** `O(1)`
- **Boolean merge:** `O(N_struts²)` naive; better with spatial partitioning; but **mesh Booleans are numerically fragile** (see [[SDF vs Mesh Booleans]])
- **Watertight:** **hard** — this is exactly what SDF methods exist to avoid

## Three reference scales

Fixed design: **cubic lattice, `L_c = 2.4 mm`, `r = 0.264 mm` (so `r* = 0.11`), `h = r/3 = 0.09 mm`, `k = 3`.** Only the primitive volume `V` varies across rows. Memory budgets assume f32 field values; sparse adds ~2× for brick addressing overhead. Woodward Fig. 4A is the Medium row.

| Scale | V (mm³) | Approx. size | Dense `G` | Sparse `G_s` | Dense memory | Sparse memory | Raw tris | Raw STL |
|---|---|---|---|---|---|---|---|---|
| Small | 10⁴ | ~22 mm cube | 13.7 M | 1.1 M | ~55 MB | ~9 MB | ~2 M | ~100 MB |
| Medium | 10⁵ | 50 mm cyl (Fig. 4A) | 137 M | 10.7 M | ~550 MB | ~85 MB | ~22 M | ~1.1 GB |
| Large | 10⁶ | 100 mm cube (Form 3 build vol) | 1.37 B | 107 M | ~5.5 GB | ~860 MB | ~210 M | ~10 GB |

**Observation: dense memory grows with volume; sparse memory grows with surface area.** The Large case is untenable dense (5.5 GB in the field buffer alone exceeds most consumer cards' free budget after the OS / driver overhead) but fits comfortably sparse.

## The dense-to-sparse ratio is scale-invariant

A closed form falls out of the formulas above:

```
G_dense / G_sparse = L_c² / (2π r k h)
```

With the design rule `h = r/3`:

```
G_dense / G_sparse = (3 / (2π k)) × (L_c / r)² = (3 / (2π k)) × (1/r*)²
```

For cubic (`k = 3`) this simplifies to **`(1/2π) × (1/r*)²`**. The ratio depends only on cell topology and dimensionless strut radius — **not on primitive size**. For typical Woodward `r*` values:

| `r*` | Cubic dense/sparse ratio |
|---|---|
| 0.05 (thin) | ~64× |
| 0.07 | ~32× |
| 0.11 (Fig. 4A) | ~13× |
| 0.15 | ~7× |
| 0.20 (thick) | ~4× |

Thinner struts give greater sparse savings — most of the volume is void, surface is a vanishing fraction. Thicker struts approach the solid case where dense and sparse converge.

## Effective sparse count (including halo)

The `S/h²` formula counts only voxels where the surface actually passes through. Real implementations include a halo of 2–3 voxels around the surface for:

- Gradient computation (central differences need ±1 voxel)
- Safety margin at brick boundaries (±1-voxel shared apron — see [[Isosurface Extraction Methods]])
- MC corner reads (a surface voxel needs 8 corners; some may fall in "inactive" neighbors)

**Effective active voxel count is typically 2–3× the `G_s` theoretical minimum.** For Medium: ~22–33 M actively processed voxels, not 11 M. Memory scales correspondingly (~170–260 MB for the Medium case); the argument for sparse over dense is unchanged but the numbers in the table above are lower bounds.

## Output mesh scale and the need for decimation

`M ≈ 2 × G_s` on average (MC emits 0–5 triangles per surface voxel, mean ~2 across a whole mesh). At 50 B/tri (binary STL, redundant vertices):

| Scale | Raw tris | Raw STL | After QEM to `h_edge = 0.1 mm` | Compressed 3MF |
|---|---|---|---|---|
| Small | ~2 M | ~100 MB | ~200 k | ~5 MB |
| Medium | ~22 M | ~1.1 GB | ~2.2 M | ~50 MB |
| Large | ~210 M | ~10 GB | ~21 M | ~500 MB |

**Raw STL is untenable for PreForm ingestion past the Small scale.** Decimation to printer-resolution edge length — see [[Mesh Decimation]] — collapses output to ~10% of the raw count with no visible loss at the physical printer resolution. This is why decimation is non-optional in the production pipeline.

## Per-strut sanity check (Medium scale)

Cross-validation of the aggregate formula against a direct strut count:

- Cells: `V/L_c³ = 10⁵ / 13.824 ≈ 7,231`
- Struts: `k × cells ≈ 21,700` (close to Woodward's 21,000)
- Per-strut lateral area: `2π × 0.264 × 2.4 ≈ 3.98 mm²`
- Total surface: `21,700 × 3.98 ≈ 86,400 mm²` — matches the aggregate `2π r k V / L_c²` formula within rounding ✓
- Surface voxels: `86,400 / 0.0081 ≈ 10.7 M` ✓

## Summary

| Method | Memory | Time | Watertight? | 2-manifold? | Mesh size |
|---|---|---|---|---|---|
| A. Dense MC33 | `O(V/r³)` ~100s MB | `O(kV/r³)` | ✓ | ✓ | large |
| B. Sparse MC33 | `O(S/r²)` ~10s MB | `O(kS/r²)` | ✓ | ✓ with brick care | large |
| C. Octree DC (standard) | `O(S/r²)` | `O(kS/r² log)` | ✓ | ✗ | medium |
| C'. Octree MDC | `O(S/r²)` | `O(kS/r² log)` | ✓ | ✓ | medium |
| D. Direct mesh | `O(N_struts)` ~MB | `O(N_struts²)` | hard | hard | small |

## Recommendation

**Sparse MC33 (Method B) + decimation + STL/3MF export** is the right first target:

- Memory scales with lattice surface area, not primitive volume.
- MC33 is the best-understood watertight + 2-manifold extractor; see [[Isosurface Extraction Methods]].
- Parallelizes cleanly on GPU compute per [[GPU Compute Pipeline]].
- Differential testing against a CPU reference is straightforward per [[CPU Reference Path]].

Octree MDC is a later optimization once profiling shows extraction is the bottleneck — unlikely given that SDF evaluation dominates.

Per-strut mesh Booleans are a dead end for watertightness.

## Caveats

1. **Boundary/skin stage cost ignored.** Stages D–F of the [[Lattice Generation Pipeline]] add `O(V^(2/3)/h²)` boundary cost, dominated by `S/h²` for high porosities. Asymptotically absorbed into the existing numbers.
2. **GPU memory pressure is real.** Medium fits on 16 GB cards with room to spare; Large (10⁶ mm³ primitive) needs sparse storage and careful brick streaming to stay under ~1 GB sparse field + mesh. Dense is infeasible at Large.
3. **Output mesh size dominates at Large scale.** For Medium, field is ~85 MB and raw output is ~1.1 GB — a ~13× asymmetry — and Large is a full order of magnitude worse. Decimation ([[Mesh Decimation]]) is required, not optional, for any pipeline past the Small scale.
4. **Numbers are upper bounds.** Real lattices share surface at node junctions and truncated boundary struts; expect 2–4× better in practice. The `~2–3×` halo factor (see "Effective sparse count" above) partially offsets this.
5. **Scale invariance of the ratio is not scale invariance of the absolute cost.** The dense/sparse ratio is ~13× regardless of scale at `r* = 0.11`, but the absolute cost grows linearly with V. A Large-scale lattice in sparse storage is still larger in absolute terms than a Medium-scale lattice in dense storage.

## Related

- [[Lattice Generation Pipeline]]
- [[Isosurface Extraction Methods]]
- [[SDF vs Mesh Booleans]]
- [[Mesh Decimation]]
- [[GPU Compute Pipeline]]
- [[AM Platforms]]
