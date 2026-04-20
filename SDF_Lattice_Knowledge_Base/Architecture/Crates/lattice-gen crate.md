---
title: lattice-gen crate
tags: [crate, architecture, lattice]
crate-path: crates/lattice-gen
---

# `lattice-gen` crate

The **pipeline library.** Unit cell topologies, lattice body construction, and the Woodward seven-stage generation process. Built on [[sdf crate]]. Pure Rust, no wgpu.

## Responsibilities (what goes IN)

- Unit cell topologies — Cubic, Kelvin, BCCxy, RD, WP (from [[Unit Cell Topologies]])
- Lattice body composition: cell template + primitive + strut radius (possibly graded) → an [[sdf crate|ExactSdf]] representing the full lattice body
- The seven-stage [[Lattice Generation Pipeline]] as [[State Transition Types]]:
  - `PrimitiveSpec → VoxelGrid<FullVox> → VoxelGrid<BoundaryVox> → LatticeSpec<Populated> → LatticeSpec<BoundaryPopulated> → LatticeSpec<Trimmed> → LatticeSpec<Skinned> → Mesh`
- CPU-side meshing for validation (slow but correct — marching cubes reference implementation)
- Geometric property queries — open porosity ε_o, window diameter d_w, specific surface area (for the [[Pressure Drop Correlation]] downstream consumer)
- Printability gating (feature sizes vs. platform from [[AM Platforms]])
- Optional dimensional compensation per [[Printing Artifacts and Compensation]]
- Typed errors at the API boundary

## What stays OUT

- `wgpu` — no GPU types. The [[gpu crate]] accelerates this crate's pipeline, but this crate does not know about it.
- `winit` — no window code.
- The SDF primitives themselves (they live in [[sdf crate]])
- WGSL source — shaders live in [[gpu crate]]'s `shaders/` directory
- CLI argument parsing / file I/O

## Why this matters

`lattice-gen` is the **specification of what is computed**. [[gpu crate]] is the specification of **how it is accelerated**. Keeping them separate means:

- The pipeline is testable on CPU without any GPU code existing.
- The CPU reference implementation *is* the spec the GPU implementation must match (see [[CPU Reference Path]]).
- Unit cell topology additions don't touch GPU code (they're new SDF compositions).
- GPU backend refactors don't touch pipeline logic (same composed SDF, different executor).

This is the pipeline half of the workspace's pure-CPU domain layer — companion to [[sdf crate]], which holds the math half. Earlier iterations of the architecture collapsed these into a single `core` crate; the split was made to keep the SDF math independently reusable and the lattice pipeline narrowly focused.

## Public API shape

```rust
pub struct LatticeJob {
    pub primitive: PrimitiveShape,
    pub cell: UnitCell,
    pub strut: StrutSpec,          // radius, optionally a grading field
    pub grid: GridSpec,            // voxel resolution
}

pub enum PrimitiveShape {
    Cylinder { radius: Millimeters, height: Millimeters },
    Cube { side: Millimeters },
    Sphere { radius: Millimeters },
    // ... or a user-provided boxed ExactSdf
}

pub enum UnitCell {
    Cubic { length: Millimeters },
    Kelvin { length: Millimeters },
    // ... per Unit Cell Topologies
}

/// The lattice body expressed as an SDF — before meshing.
pub fn lattice_body(job: &LatticeJob) -> impl sdf::ExactSdf { /* ... */ }

/// Full CPU pipeline: produces a watertight mesh. Slow but correct.
pub fn generate_mesh(job: &LatticeJob) -> Result<Mesh, GenError> { /* ... */ }
```

Note: `lattice_body(...)` is the seam the [[gpu crate]] uses. It returns the composed SDF; the GPU backend then runs its own fast evaluation + [[Isosurface Extraction Methods|MC33 extraction]] over that SDF without re-implementing the lattice composition.

## What this crate consumes

- [[sdf crate]] — primitives, combinators, repetition, precision types
- `glam` for vector math (inherited via sdf)
- `thiserror` for typed errors
- `meshopt` (or equivalent) for CPU-side decimation reference
- `serde` for job serialization

## What this crate exposes

- `LatticeJob` and related input types
- `Mesh` output type (vertex + index buffers, plain data)
- The composable `lattice_body()` function — the seam for GPU acceleration
- Geometric query functions (porosity, window diameter, surface area)
- Typed errors — no `sdf` or internal errors leak directly

## Phasing

- **Phase 1 (first implementation):** CPU-only. Cubic topology only. Small-scale meshes (Woodward Fig. 4B-scale). No compensation. Establishes the pipeline spec and test suite.
- **Phase 2:** [[gpu crate]] arrives; GPU backend uses `lattice_body()` as input, runs field + MC33 on wgpu. Differential test vs. CPU reference.
- **Phase 3:** Additional topologies (Kelvin, BCCxy, RD, WP). Compensation pipeline. Larger scales.

## Testing strategy

Per [[Three Levels of Tests]]:

- **Unit.** Each pipeline stage tested in isolation on known inputs.
- **Integration.** Full CPU pipeline on small lattices; vertex count, Euler characteristic, manifold invariants asserted.
- **Property.** For any input job: output mesh is watertight + 2-manifold + orientation-consistent. Per [[Sharp Oracles]].
- **Differential.** CPU pipeline on large job vs. [[gpu crate]] GPU pipeline on same job (Phase 2+).
- **Adversarial** per [[Tests Should Make Programs Fail]]: tiny primitives, cells larger than primitive, strut radius approaching cell size, degenerate primitives (flat, thin, concave), extreme porosity.

## Related

- [[Core Principles]]
- [[CPU Reference Path]]
- [[Lattice Generation Pipeline]]
- [[Unit Cell Topologies]]
- [[Printing Artifacts and Compensation]]
- [[sdf crate]]
- [[gpu crate]]
- [[State Transition Types]]
