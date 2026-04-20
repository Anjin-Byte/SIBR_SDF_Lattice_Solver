---
title: sdf crate
tags: [crate, architecture, sdf]
crate-path: crates/sdf
---

# `sdf` crate

The **math library.** Exact signed-distance functions, combinators, operators, and deformations. Pure Rust, no wgpu, no wgsl runtime. This is the crate `lattice-gen` builds on and `gpu` ports to WGSL.

## Responsibilities (what goes IN)

- The SDF trait hierarchy — precision tracked at the type level
- Primitive types (Sphere, Box, Capsule, Cylinder, Torus, …) — the [[SDF Primitives Catalog]] contents
- Combinators: Union, Intersection, Subtraction, Xor (exact-preserving); SmoothUnion and friends (bound-producing)
- Operators: Translate, UniformScale, SymX/Y/Z, Repeat, LimitedRepeat (precision-preserving)
- Deformations: Twist, Bend, Displace, Onion, Round (always bound-producing)
- CPU evaluation for all of the above — this is the [[CPU Reference Path]] oracle
- Typed errors at the API boundary
- Deterministic test fixtures and property tests

## What stays OUT

- `wgpu` — no GPU types, ever
- `winit` — no window code
- Lattice-specific logic (unit cell topologies, population, skin generation — all lives in [[lattice-gen crate]])
- WGSL shader source — lives in [[gpu crate]] under `crates/gpu/shaders/sdf/`. Differential tests live in cross-crate tests.
- CLI / file I/O

## Why this matters

Three clean properties fall out of this split:

1. **Unit-testable in full isolation.** No GPU, no lattice context. Every primitive has a closed-form correctness oracle — `|p| - r` for a sphere, etc.
2. **Reusable.** A future non-lattice SDF application (e.g., a rendering tool, a collision system) can consume this crate without touching lattice code.
3. **Refactor-safe.** When the [[lattice-gen crate]] pipeline changes, the SDF math doesn't. When the SDF math changes, the pipeline's public API doesn't (only its internal body of SDF ops).

This is the math half of the workspace's pure-CPU domain layer — companion to [[lattice-gen crate]], which holds the pipeline half. Earlier iterations of the architecture collapsed these into a single `core` crate; the split was made when it became clear the SDF math is generally reusable and the lattice pipeline is specific to this project.

## Precision as a type-level concern

The [[SDF Primitives Catalog]] identifies a safety-critical distinction between:

- **Exact SDFs** — return the true Euclidean distance; 1-Lipschitz; safe inside Boolean ops
- **Bound SDFs** — return `≤` the true distance; usable for sphere tracing with step discount; **unsafe** inside Boolean ops (produce wrong distances silently)

This crate enforces the distinction via marker traits so that contamination is a compile error, not a runtime bug. Tentative shape:

```rust
pub trait Sdf {
    fn eval(&self, p: glam::Vec3) -> f32;
}

/// Marker: `eval` returns the true Euclidean distance.
pub trait ExactSdf: Sdf {}

/// Marker: `eval` returns a value `≤` true distance (bound only).
pub trait BoundSdf: Sdf {}

// Combinators preserve precision when inputs are exact
impl<A: ExactSdf, B: ExactSdf> Sdf for Union<A, B> { /* ... */ }
impl<A: ExactSdf, B: ExactSdf> ExactSdf for Union<A, B> {}

// Smooth combinators always produce bounds
impl<A: Sdf, B: Sdf> Sdf for SmoothUnion<A, B> { /* ... */ }
impl<A: Sdf, B: Sdf> BoundSdf for SmoothUnion<A, B> {}
```

A lattice generator that requires exact SDFs for correctness can be generic over `T: ExactSdf` — the compiler rejects bound inputs. See [[Make Invalid States Unrepresentable]].

## Public surface guidance

- Types are plain data plus trait impls; avoid hidden state.
- No `async` — pure math, sync evaluation.
- Errors via `thiserror`, typed per-module.
- `glam` for vector math (small, fast, widely adopted in the Rust graphics ecosystem).

## What this crate depends on

- `glam` for Vec3/Mat4 types
- `thiserror` for typed errors
- `serde` (with `derive`) for SDF-tree serialization (for CLI inputs eventually)
- Small test-only deps: `proptest` or similar for property tests

**No** wgpu, winit, tokio, or anything async / platform-specific.

## Testing strategy

Per [[Three Levels of Tests]]:

- **Unit.** Each primitive against its closed-form definition on sampled points.
- **Property.** For every exact primitive: `|f(p)| = dist(p, surface)` for sampled surface points; 1-Lipschitz property `|f(p) - f(q)| ≤ |p - q|`.
- **Differential (cross-crate).** CPU eval vs. WGSL eval (the [[CPU Reference Path]] oracle for [[gpu crate]] shader correctness). Lives in `tests/` at workspace root.
- **Adversarial.** Degenerate inputs (zero radii, coincident endpoints, points on the surface, at infinity, inside discontinuities). Per [[Tests Should Make Programs Fail]].

## Related

- [[Core Principles]]
- [[GPU Boundary]]
- [[CPU Reference Path]]
- [[SDF Primitives Catalog]]
- [[lattice-gen crate]]
- [[gpu crate]]
- [[Make Invalid States Unrepresentable]]
