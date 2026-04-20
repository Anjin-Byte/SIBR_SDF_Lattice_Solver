---
title: Isosurface Extraction Methods
tags: [domain, meshing, extraction]
---

# Isosurface Extraction Methods

The algorithms for converting a sampled SDF `f: ℝ³ → ℝ` into a triangle mesh of the surface `{f = 0}`. This note compares the three algorithms relevant to the project — Marching Cubes, Dual Contouring, Manifold Dual Contouring — and justifies the recommendation for MC33.

## Watertight vs. 2-manifold — the crucial distinction

Two properties often conflated:

- **Watertight** = the mesh is closed; no boundary edges; bounds a volume.
- **2-manifold** = locally homeomorphic to a disk; every edge shared by exactly two faces; every vertex's one-ring is a topological disk.

A mesh can be watertight but non-manifold (e.g., two surfaces touching at an edge, sharing three or four faces per edge). Non-manifold meshes break most slicer software, including PreForm — see [[AM Platforms]]. **We need both properties.**

## Marching Cubes (MC33, Chernyaev 1995)

The canonical algorithm. For each grid cell:

1. Read 8 corner field values; compute 8-bit sign index (`Σ (f_i < 0) << i`, 256 cases).
2. Look up triangle count + edge configuration from a precomputed table.
3. Interpolate vertex positions along sign-change edges.
4. Emit triangles.

**Guarantees with MC33 (Chernyaev's disambiguated variant):**

- ✓ Watertight
- ✓ 2-manifold edges
- ✓ 2-manifold vertices
- ✓ No ambiguous-case holes (the original Lorensen-Cline 1987 MC had these; MC33 fixed them with additional face/interior tests)

**Drawbacks:**

- Doesn't preserve sharp features — interpolated vertices smooth over creases.
- No adaptive refinement; fixed `h` everywhere.
- Produces slivers at grazing sign-changes (harmless for printing, ugly for rendering).
- Raw output has duplicate vertices at voxel edges; needs a welding pass for compact indexed output.

**Reference:** Chernyaev, E. "Marching Cubes 33: Construction of Topologically Correct Isosurfaces." Technical Report CN/95-17, CERN, 1995.

## Standard Dual Contouring (DC, Ju et al. 2002)

Dual to MC — places vertices **inside cells**, not on edges.

For each sign-change cell:

1. Gather *Hermite data* — surface-intersection points on sign-change edges, plus surface normals at those points.
2. Place one vertex at the position minimizing the Quadratic Error Function (QEF) over the Hermite data.

For each sign-change edge: emit a quad connecting the four vertices of the four cells sharing that edge.

**Advantages:**

- ✓ Preserves sharp features (the QEF-minimizing position lands exactly on creases).
- ✓ Adaptive-octree-friendly.
- ✓ Fewer triangles per feature than MC.

**Guarantees:**

- ✓ Watertight
- ✗ **Not 2-manifold** in general

**Why not 2-manifold:**

1. **Multi-component cells.** A cell with two disjoint surface patches (e.g., two near-parallel struts) gets one vertex — both patches share it, spuriously joining surfaces.
2. **Shared-edge topology.** The dual quad may connect vertices whose inside/outside labels don't form a consistent manifold ring, producing non-manifold edges.

**Reference:** Ju, T., Losasso, F., Schaefer, S., Warren, J. "Dual Contouring of Hermite Data." SIGGRAPH 2002.

## Manifold Dual Contouring (MDC, Schaefer et al. 2007)

Fixes DC's manifold problem by allowing **more than one vertex per cell** — one per disconnected surface component in that cell.

**Algorithm sketch:**

1. For each cell, use MC-style face topology tables to determine which sign-change edges belong to the same surface component within the cell.
2. Place one vertex per component (each still QEF-minimizing).
3. Generate dual polygons with component-aware connectivity: connect vertices belonging to the same component at each sign-change edge.
4. Post-pass: walk each vertex's one-ring; if non-manifold, split the vertex and reassign faces.

**Guarantees:**

- ✓ Watertight
- ✓ 2-manifold edges
- ✓ 2-manifold vertices
- ✓ Surface-component-faithful

**Cost over standard DC:**

- Per-cell component classification (cheap — 256-entry tables like MC).
- Up to ~4 vertices per cell in pathological configurations (typically 1–2).
- Vertex-splitting pass is local (O(1) per non-manifold vertex).
- Net: ~1.5–2× slower than standard DC; 10–20% more vertices.

**Reference:** Schaefer, S., Ju, T., Warren, J. "Manifold Dual Contouring." *IEEE TVCG* 13(3), 2007, pp. 610–619.

## Comparison

| Property | MC33 | Standard DC | Manifold DC |
|---|---|---|---|
| Watertight | ✓ | ✓ | ✓ |
| 2-manifold edges | ✓ | ✗ | ✓ |
| 2-manifold vertices | ✓ | ✗ | ✓ |
| Sharp feature preservation | ✗ | ✓ | ✓ |
| Adaptive octree | — | ✓ | ✓ |
| GPU-compute friendly | ✓✓ | ✓ (octree nav) | ✓ (octree nav) |
| Table size | 256 entries | small | 256 entries + DC tables |
| Mesh quality (triangles) | skinny slivers possible | regular quads | regular quads |

## Recommendation for the lattice solver

**MC33 on a sparse bricked grid**, for these reasons:

1. **Lattice surfaces are intrinsically smooth.** Capsule/cylinder primitives have no sharp features. Smooth-union joints (where we use them) are deliberately smooth. DC's sharp-feature-preservation advantage is not needed.
2. **2-manifold output is non-negotiable** for PreForm / slicer ingestion. MC33 delivers it by construction with a well-understood 256-entry table. Standard DC would need the full MDC machinery.
3. **Parallelism is maximal.** MC33 is per-voxel independent; no octree navigation; no inter-cell dependencies. Maps cleanly to wgpu compute.
4. **Reference implementations are abundant.** CPU reference for differential testing ([[CPU Reference Path]]) is trivially available.
5. **The only cost — over-tessellation of flat regions — is handled by [[Mesh Decimation]]** as a post-process.

If a future feature set adds hybrid solid-lattice parts with sharp creases at the interface, revisit: extended MC variants handle that without needing octree adaptivity. MDC becomes relevant only if both (a) sharp features matter *and* (b) octree adaptivity is needed for memory — neither of which applies to our current application.

## Brick-seam care (a subtle watertightness concern)

MC33 is 2-manifold within a single contiguous grid. On a bricked grid (Method B in [[Meshing Complexity Analysis]]), adjacent bricks must **share the field values at their boundary planes** — typically via a one-voxel apron each brick holds of its neighbors' data. Without this, two bricks MC-extract boundary triangles at slightly different vertex positions, leaving sub-voxel cracks. The mesh is still topologically closed, but with non-manifold near-coincident edge pairs at the cracks.

The fix is textbook (every VDB-style library handles it) but it's a real implementation detail that belongs in the gpu crate's meshing module, not as an afterthought.

## Related

- [[Meshing Complexity Analysis]]
- [[Mesh Decimation]]
- [[GPU Compute Pipeline]]
- [[SDF Primitives Catalog]]
- [[CPU Reference Path]]
- [[Sharp Oracles]]
