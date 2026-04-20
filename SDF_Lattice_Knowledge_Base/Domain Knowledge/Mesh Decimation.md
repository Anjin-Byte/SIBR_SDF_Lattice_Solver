---
title: Mesh Decimation
tags: [domain, meshing, post-process]
---

# Mesh Decimation

The post-process that reduces triangle count while preserving geometric fidelity. Required after [[Isosurface Extraction Methods|MC33 extraction]] because the extractor over-tessellates uniformly regardless of local curvature.

## Why decimation is non-optional

Marching Cubes produces **~2 triangles per surface voxel**, everywhere. A dead-flat strut wall receives the same triangle density as a node junction. The printer can't resolve detail finer than its slice thickness + spot size anyway (see [[AM Platforms]]), so sub-printer-resolution triangles are decoration.

Concrete numbers for the [[Meshing Complexity Analysis|Woodward Fig. 4A cylinder]]:

| Stage | Triangle count | Uncompressed STL |
|---|---|---|
| Raw MC33 output, `h = 0.09 mm` | ~30 M | ~1.5 GB |
| After QEM to edge-length ≈ 0.1 mm | ~3–5 M | ~200 MB |
| After QEM + 3MF compression | same tri count | ~20–40 MB |

10× reduction with no visible difference at printer resolution. A shippable pipeline decimates; unoptimized STLs waste slicer time and disk.

## Quadric Edge Collapse (QEM)

Garland & Heckbert, 1997. The standard algorithm; every serious mesh library implements it.

**Idea.** Repeatedly pick an edge, collapse its two endpoints into one vertex, delete the two adjacent triangles. Pick edges in order of geometric error so the mesh shape is preserved.

**Algorithm:**

1. For each vertex `v`, compute a 4×4 symmetric *quadric matrix* `Q_v` encoding the sum of squared distances to the planes of its adjacent faces. Closed-form; no fitting.
2. For each edge `(v_a, v_b)`, `Q_edge = Q_a + Q_b`. The optimal collapsed position `v*` minimizes `v*ᵀ Q_edge v*` — one 3×3 linear solve.
3. Collapse *cost* is `v*ᵀ Q_edge v*` at the optimum. Interpretation: accumulated squared distance from the original surface that this collapse introduces.
4. Priority-queue all edges by cost. Pop the cheapest; collapse it; update `Q` for the new vertex; recompute costs of affected edges; reinsert.
5. Stop on any of: target triangle count, error threshold, shortest-edge-length threshold.

**Why it's good:**

- **Geometry-aware.** Flat regions (low quadric error) decimate aggressively — long, thin triangles survive where they're fine. High-curvature regions (strut junctions, skin features) retain triangles — quadric cost is high there. No parameters per region; it self-adapts.
- **Fast.** `O(E log E)` in mesh edge count.
- **Topology-preserving.** With standard guards, preserves the watertight + 2-manifold properties of the input.

## Manifoldness preservation — required guards

Naive QEM can create:

- Foldovers (triangle normals flip)
- Non-manifold edges (an edge shared by >2 faces after collapse)
- Self-intersections (geometrically distant parts become close)

Production QEM implementations refuse any collapse that would:

1. Create an edge with >2 incident faces.
2. Flip a triangle's normal relative to its original orientation (detects foldovers).
3. Produce a non-disk vertex one-ring.

These checks are O(1) per candidate collapse. Any serious library ships them.

**Rust option:** the [`meshopt`](https://crates.io/crates/meshopt) crate (bindings to Arseny Kapoulkine's `meshoptimizer` C++ library) exposes:

- `simplify(indices, vertices, target_count, target_error)` — topology-preserving, safe for printing.
- `simplify_sloppy(...)` — faster but can change topology; **do not use for printable output**.

## Stopping criterion for our pipeline

Use **shortest-edge-length** as the stopping criterion, not triangle count:

```
target_edge_length = max(slice_thickness, xy_spot_size × scale_factor)
```

Concrete for Formlabs SLA:
- Slice thickness: 0.05–0.1 mm (user-chosen in PreForm)
- Form 3/3+ spot size: 25 µm → ~0.05 mm
- `target_edge_length ≈ 0.1 mm` typically

This directly expresses "don't produce detail below the printer's ability to print" — the physical constraint, not an arbitrary triangle budget.

## Implementation notes for the project

- **Lives in [[lattice-gen crate]], not [[gpu crate]].** Operates on geometry (vertex + index arrays), no wgpu types. Pure CPU; no GPU dependency.
- **Single-threaded is fine.** For tens of millions of triangles, QEM is tens of seconds on one core — negligible compared to interactive use and tiny compared to the hours of print time. Parallelization is not first-pass work.
- **Runs after extraction, before export.** Pipeline shape: `SDF → MC33 (fine) → QEM (target printer resolution) → STL/3MF export`.
- **Don't try to skip the fine extraction.** Extracting directly at coarse resolution misses features smaller than `h`. Extract fine, decimate to target.

## Invariants the decimation stage must preserve

1. **2-manifold in, 2-manifold out.** QEM with the standard guards maintains this.
2. **Watertight in, watertight out.** Closed mesh stays closed.
3. **Self-intersection free in, self-intersection free out.** QEM with foldover detection maintains this; the "sloppy" variants do not.
4. **Orientation-consistent in, orientation-consistent out.** Face normals stay outward-pointing.

Any decimation stage that cannot document all four under its specific configuration is unacceptable for our pipeline. These are the conditions PreForm and similar slicers require of input ([[AM Platforms]]).

## Testing

- **Differential:** compare undecimated and decimated meshes on a library of reference lattices. Assert Hausdorff distance below the target tolerance.
- **Invariant:** on every decimated mesh, assert (a) closed, (b) 2-manifold, (c) triangle count reduced, (d) no flipped normals. These are sharp oracles per [[Sharp Oracles]].
- **Regression:** collect any mesh that ever produced non-manifold output from decimation; keep as a permanent failing-then-passing test per [[Regression Discipline]].
- **Adversarial:** test near-pathological inputs — slivers at grazing sign-changes, pinch points, near-zero-volume components. Per [[Tests Should Make Programs Fail]].

## Related

- [[Meshing Complexity Analysis]]
- [[Isosurface Extraction Methods]]
- [[GPU Compute Pipeline]]
- [[AM Platforms]]
- [[lattice-gen crate]]
- [[Sharp Oracles]]
