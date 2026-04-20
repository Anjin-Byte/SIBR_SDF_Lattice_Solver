---
title: SDF vs Mesh Booleans
tags: [domain, sdf, meshing, terminology]
---

# SDF vs Mesh Booleans

Two operations share the name "Boolean" in 3D geometry. They have almost nothing in common beyond both producing a union/intersection/difference of two shapes. **Conflating them is a common source of performance misunderstandings**, so the distinction is canonicalized here.

## SDF Boolean

**Operation on distance functions.** Given two SDFs `a, b: ℝ³ → ℝ`, produce a new SDF representing the Boolean combination of their zero-level sets.

```text
opUnion(a, b)        = min(a, b)              // ~1 ALU op
opIntersection(a, b) = max(a, b)              // ~1 ALU op
opSubtraction(a, b)  = max(-a, b)             // ~2 ALU ops
opXor(a, b)          = max(min(a,b), -max(a,b))  // ~4 ALU ops
```

**Smooth variants** (produce bounds, not exact SDFs — see [[SDF Primitives Catalog]] "Bound contamination"):
- `opSmoothUnion(a, b, k)` — blended min via `smin`
- `opSmoothIntersection`, `opSmoothSubtraction` — symmetric forms

**Cost profile:**

- **Per-query cost:** sum of the operand evaluation costs + a handful of scalar ops.
- **Parallelism:** point-wise; no cross-query dependencies; no memory traffic beyond what the operands need.
- **Hardware fit:** uniform across SIMD lanes; no branch divergence (`min`/`max` are hardware instructions); parallelizes trivially per voxel / per thread / per lane.
- **Failure modes:** essentially none — closed-form, deterministic, floating-point-stable for well-behaved operands.

**This is the operation every function in the [[SDF Primitives Catalog]] uses to compose primitives into the lattice body.**

## Mesh Boolean

**Operation on triangle meshes (also called CSG, constructive solid geometry).** Given two triangulated surfaces `A, B`, produce a new triangulated surface representing their Boolean combination.

**Algorithm structure:**

1. **Compute intersection curves** between faces of `A` and faces of `B`. Each pair of intersecting triangles produces a line segment; globally these form closed loops.
2. **Refine both meshes** by splitting faces along the intersection curves.
3. **Classify** each resulting sub-face as *inside A*, *outside A*, *inside B*, *outside B*. Use point-in-polyhedron tests.
4. **Select** sub-faces according to the Boolean operation (union keeps *outside A and inside B* complement logic, etc.).
5. **Stitch** the selected faces into a consistent manifold, handling edge pairings and orientations.

**Cost profile:**

- **Per-operation cost:** `O(N_A × N_B)` for naive intersection testing; `O((N_A + N_B) log N)` with spatial partitioning; intersection-curve computation is expensive per triangle pair.
- **Parallelism:** pairwise intersection computation parallelizes; **stitching is serial-ish** (has data dependencies on edge equivalence classes).
- **Numerical robustness:** infamously fragile. Coplanar faces, near-grazing intersections, vertices exactly on faces — all produce either wrong topology or outright crashes. Industrial CGAL-style implementations use **exact arithmetic** (rational or constructive reals) to survive, at substantial cost.
- **Failure modes:** self-intersections, non-manifold edges, lost faces, disconnected components, inside-out orientations. The literature on robust mesh Booleans is a 40-year exercise in dodging these.

## Quick comparison

| Property | SDF Boolean | Mesh Boolean |
|---|---|---|
| Operates on | scalar functions | triangle meshes |
| Per-query cost | 1–10 ALU ops | ~1k–10k ALU ops + memory access |
| Global cost per op | N/A (point-wise) | `O(N_A N_B)` naive; `O(N log N)` with BVH |
| Parallelism | trivially embarrassing | harder (stitching phase has dependencies) |
| Numerical robustness | closed-form, stable | notoriously fragile |
| Failure modes | essentially none | self-intersection, non-manifold, lost faces |
| Production implementations | ~50 lines of GLSL/WGSL | CGAL, libigl, OpenSCAD (all non-trivial) |

## Why this matters for the project

The entire reason for choosing an SDF-based approach over a CAD-style mesh-CSG approach is to **avoid mesh Booleans**. The pipeline shape is:

1. **Compose in SDF space** — primitives combined via SDF Booleans ([[SDF Primitives Catalog]]).
2. **Extract once, at the end** — a single marching-cubes pass produces a watertight manifold mesh from the final SDF ([[Isosurface Extraction Methods]]).
3. **No mesh Booleans anywhere in the pipeline.**

This is what makes the approach both performant and watertight — the combinatorial explosion of mesh-Boolean robustness issues never happens because we never do them. See also [[Lattice Generation Pipeline]] and [[Meshing Complexity Analysis]] Method D (the "direct per-strut mesh Boolean" anti-pattern we explicitly reject).

## Confusion to avoid

When someone (including me) says "Boolean" in a lattice / SDF / CAD context, **always disambiguate**:

- In a voxel grid or sphere-tracing loop → **SDF Boolean**, cheap.
- In a CAD model tree (Rhino, Grasshopper, OpenSCAD, CGAL) → **mesh Boolean**, expensive and fragile.
- In the Woodward paper's generation pipeline ([[Lattice Generation Pipeline]]) → *mesh Boolean* (they use Grasshopper + Dendro). Our implementation replaces this with SDF Booleans throughout.

When reviewing code or discussing performance, if the answer hinges on "is this Boolean cheap or expensive," pause and ask which kind of Boolean.

## Related

- [[SDF Primitives Catalog]]
- [[Meshing Complexity Analysis]]
- [[Isosurface Extraction Methods]]
- [[GPU Compute Pipeline]]
- [[Explicit Not Magical]]
