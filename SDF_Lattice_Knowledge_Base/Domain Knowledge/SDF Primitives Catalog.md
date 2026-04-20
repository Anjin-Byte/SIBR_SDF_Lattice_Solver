---
title: SDF Primitives Catalog
tags: [domain, sdf, primitives, gpu]
source: reference/iq_distfunctions_exact_from_uploaded_html.md
source-url: https://iquilezles.org/articles/distfunctions/
---

# SDF Primitives Catalog

The kernel library of signed-distance functions that the [[gpu crate]] WGSL shaders will port. Source is Inigo Quilez's distance-functions article (HTML export at `reference/iq_distfunctions_exact_from_uploaded_html.md`). The original GLSL is the authoritative form; ported WGSL is verified by the [[CPU Reference Path]].

## Definitions and contract

A **signed distance function** `f: R³ → R` for a closed surface S satisfies:

- `f(p) = 0` for `p ∈ S`
- `f(p) > 0` for `p` outside S
- `f(p) < 0` for `p` inside S
- `|f(p)| = distance(p, S)` (exact case)
- `|f(p)| ≤ distance(p, S)` (bound case — useful for sphere tracing but with reduced step size)

**Critical invariant for sphere tracing:** the function must be 1-Lipschitz: `|f(p) - f(q)| ≤ |p - q|`. Operations that violate this (most deformations) require either a discount factor on the step size or restriction to non-tracing uses.

## Exact primitives

Solid shapes whose SDF is the true Euclidean distance:

| Primitive | Inputs | Notes |
|---|---|---|
| Sphere | radius `r` | The trivial example: `length(p) - r` |
| Box | half-extents `b` | Axis-aligned |
| Round Box | `b`, radius `r` | Box minus inflated negative-r |
| Box Frame | `b`, edge thickness `e` | Three-axis combination |
| Torus | `(R_major, r_minor)` | XZ-plane oriented |
| Capped Torus | sin/cos of cap angle, `ra`, `rb` | Slice of a torus |
| Link | length `le`, `r1`, `r2` | Two linked tori |
| Cylinder (capped) | `r`, `h` | Vertical, finite |
| Arbitrary Capped Cylinder | endpoints `a`, `b`, `r` | Aligned along arbitrary axis |
| Rounded Cylinder | `ra`, `rb` (round), `h` | Smoothed edges |
| Cone (finite) | sin/cos `c`, height `h` | Apex up |
| Capped Cone | `h`, `r1`, `r2` | Two radii (frustum) |
| Capped Cone (arbitrary) | endpoints `a`, `b`, radii `ra`, `rb` | Aligned along arbitrary axis |
| Round Cone | `r1`, `r2`, `h` | Spherically capped frustum |
| Round Cone (arbitrary) | endpoints `a`, `b`, `r1`, `r2` | Aligned along arbitrary axis |
| Capsule (line) | endpoints `a`, `b`, `r` | Cylinder with hemisphere caps |
| Vertical Capsule | `h`, `r` | Y-axis capsule |
| Plane | normal `n` (unit), offset `h` | Half-space |
| Hexagonal Prism | `h = (radius, height)` | Six-sided |
| Cut Sphere | `r`, cut height `h` | Sphere with planar cap |
| Cut Hollow Sphere | `r`, `h`, thickness `t` | Shell with cap |
| Death Star | `ra`, `rb`, separation `d` | Sphere with spherical cavity |
| Solid Angle | sin/cos `c`, `ra` | Cone-bounded sphere wedge |
| Vesica Segment | endpoints `a`, `b`, width `w` | Lens between endpoints |
| Rhombus | `la`, `lb`, height `h`, round `ra` | Diamond prism |
| Octahedron | size `s` | Eight-faced |
| Pyramid | height `h` | Square base, unit-sized |
| Triangle (UD) | vertices `a`, `b`, `c` | **Unsigned** distance — `udTriangle` |
| Quad (UD) | vertices `a`, `b`, `c`, `d` | **Unsigned** distance — `udQuad` |

## Bound (non-exact) primitives

These functions return a value `≤` true distance. They are still useful for sphere tracing if the marcher uses a step discount, but they are **not** safe inputs to operations that depend on exactness.

- **Ellipsoid** — `k0*(k0-1)/k1` form, returns a bound only
- **Octahedron (bound variant)** — cheap form: `(p.x+p.y+p.z-s)*0.57735027`
- **Triangular Prism** — bound only

**Implication for the codebase:** the type system should distinguish exact from bound. A possible shape: `pub enum SdfKind { Exact, Bound }` carried in the SDF function's metadata, or two separate trait families. See [[Make Invalid States Unrepresentable]] — feeding a `Bound` SDF into an operation that requires `Exact` should be a compile error, not a silent correctness loss.

## 2D → 3D constructors

Take a 2D SDF primitive and lift to 3D:

- `opRevolution(p, primitive2d, offset)` — rotate around Y axis at radial offset `o`
- `opExtrusion(p, primitive2d, h)` — extrude along Z to half-height `h`

## Boolean combinators

Exact (over exact inputs):

- `opUnion(a, b) = min(a, b)`
- `opSubtraction(a, b) = max(-a, b)` — subtracts `a` from `b`
- `opIntersection(a, b) = max(a, b)`
- `opXor(a, b) = max(min(a, b), -max(a, b))`

**Smooth** combinators (return bounds, not exact distances; controlled by a smoothing parameter `k`):

- `opSmoothUnion(a, b, k)` — `min(a,b) - h²/(4k)` with `h = max(k - |a-b|, 0)` (after `k *= 4`)
- `opSmoothSubtraction(a, b, k) = -opSmoothUnion(a, -b, k)`
- `opSmoothIntersection(a, b, k) = -opSmoothUnion(-a, -b, k)`

The smooth variants are critical for blended joints in lattices but are bounds — see the **bound contamination** caveat below.

## Positioning

- `opTx(p, t, primitive)` — apply inverse transform to the query point: `primitive(invert(t)*p)`
- `opScale(p, s, primitive)` — `primitive(p/s) * s` — preserves Lipschitz; non-uniform scale would break it

## Symmetry

- `opSymX(p)` — `p.x = abs(p.x)` then evaluate primitive
- `opSymXZ(p)` — `p.xz = abs(p.xz)` then evaluate primitive
- (Generalize as needed for Y, XY, etc.)

## Repetition

- **Infinite repetition:** `opRepetition(p, s, primitive) = primitive(p - s*round(p/s))` — tiles the primitive on a grid of cell size `s`
- **Limited repetition:** `opLimitedRepetition(p, s, l, primitive) = primitive(p - s*clamp(round(p/s), -l, l))` — same but bounded to `±l` cells

These are the central operations for lattice generation. Combined with intersection against the primitive shape (per [[Lattice Generation Pipeline]]), repetition produces a tiled lattice trimmed to a bounding shape.

## Deformations

These break exact SDF properties and produce only bounds. **They are not safe to compose with Boolean ops that require exactness, and they require step discounts when sphere-traced.**

- `opDisplace(p, primitive)` — adds a per-point displacement: `primitive(p) + displacement(p)`
- `opTwist(p, primitive)` — rotates by angle proportional to `p.y`
- `opCheapBend(p, primitive)` — rotates by angle proportional to `p.x`
- `opElongate(p, primitive, h)` — stretches the primitive along chosen axes (two forms; the second is the well-behaved one)
- `opRound(primitive, rad)` — `primitive(p) - rad` — strictly safe for exact inputs (just inflates), but degrades to bound when applied to a bound primitive
- `opOnion(sdf, t)` — `abs(sdf) - t` — turns a solid into a thick shell

## Change of metric

Replace `length()` with non-Euclidean norms (`length2`, `length6`, `length8`). Notes from the source: bound only; the author "does not recommend" them. We will likely not implement these.

## Bound contamination

A function that combines an exact primitive with a bound primitive (or with a bound combinator) yields a **bound**. This propagates through the SDF graph. The type system can prevent silent contamination:

- Exact ∘ Exact via `opUnion`, `opIntersection`, `opSubtraction`, `opXor` → Exact
- Anything ∘ Smooth combinator → Bound
- Anything via deformation (twist, bend, displace, onion) → Bound
- `opScale` (uniform) → preserves Exact/Bound classification of input
- `opTx` → preserves Exact/Bound classification of input
- `opSymX`, `opSymXZ` → preserves classification (the abs is on the input, not the output)

This is the kind of invariant [[Reason About Correctness]] wants externalized at the type level.

## Project implications

1. **WGSL ports** of these functions live in `crates/gpu/shaders/`, with naming aligned to the source (`sdSphere`, `sdBox`, `opUnion`, etc. — kept verbatim where possible to preserve traceability to the original article).
2. **CPU reference implementations** of all exact primitives live in `crates/sdf` (see [[CPU Reference Path]]). The differential test between WGSL and CPU is the [[Sharp Oracles|sharp oracle]] for SDF correctness.
3. **Exactness tracking** must be preserved through the type system, not just by convention. See [[Make Invalid States Unrepresentable]].
4. **Adversarial tests** (per [[Tests Should Make Programs Fail]]) for SDFs should include: degenerate inputs (`r=0`, `h=0`, coincident endpoints), points exactly on the surface, points at infinity, points inside corners and edges where derivatives are discontinuous.

## Related

- [[Three-Layer Stack]]
- [[Lattice Generation Pipeline]]
- [[gpu crate]]
- [[CPU Reference Path]]
- [[Make Invalid States Unrepresentable]]
