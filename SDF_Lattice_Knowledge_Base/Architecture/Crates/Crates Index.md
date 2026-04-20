---
title: Crates Index
tags: [index, crates]
---

# Crates Index

Per-crate responsibility notes. Each note answers: what goes in, what stays out, and why.

## Phase 1 (starting state — CPU-only, no GPU yet)

- [[sdf crate]] — SDF math library; primitives, combinators, operators; precision tracked at the type level
- [[lattice-gen crate]] — unit cell topologies; Woodward seven-stage pipeline; CPU-side meshing
- [[xtask crate]] — build/automation tasks

## Phase 2 (GPU acceleration + CLI)

- [[gpu crate]] — wgpu backend that accelerates what `lattice-gen` specifies
- [[cli crate]] — argument parsing, file I/O, orchestration

## Phase 3 (optional)

- [[preview crate]] — winit + surface + interactive viewer

## Dependency direction

```
          ┌────────────────────────┐
          │         cli            │
          │ (argv, file I/O)       │
          └──┬─────────────────┬───┘
             │                 │
             ▼                 ▼
          ┌──────┐      ┌───────────────┐
          │ gpu  │◀────▶│  lattice-gen  │
          │(wgpu)│      │  (pipeline)   │
          └──┬───┘      └───────┬───────┘
             │                  │
             └────────▶┌────────┘
                       ▼
                    ┌──────┐
                    │ sdf  │ (math, no deps on above)
                    └──────┘
```

- **`sdf`** depends on nothing project-local. Pure math. Reusable outside this project.
- **`lattice-gen`** depends on `sdf`. Composes SDFs into lattice bodies; runs the pipeline.
- **`gpu`** depends on `sdf` (for shader parity against CPU reference) and `lattice-gen` (to execute its pipeline on GPU).
- **`cli`** depends on `gpu` and `lattice-gen`.
- **`preview`** (when added) depends on `gpu` and `lattice-gen`, never on wgpu directly.
- **`xtask`** is independent of `cli`, `preview`. May depend on `sdf` / `lattice-gen` for shared fixture types.

No crate depends on `wgpu` or `winit` except those where the vault explicitly permits it. See [[GPU Boundary]] and [[Headless First]].

## Why this split (vs. a single `core` crate)

The project's domain decomposes into **math** (SDF primitives and combinators) and **pipeline** (unit cells, Woodward's seven stages). These evolve at different rates and serve different audiences:

- `sdf` is a general-purpose library that could be lifted out and reused for any SDF application.
- `lattice-gen` is specific to this project's lattice-generation task.

Splitting them keeps `sdf` portable, keeps `lattice-gen` focused on one problem, and matches [[Workspace Layout]]'s rule: *"if something can be built, tested, or shipped independently, it deserves its own crate."*

## Related

- [[Workspace Architecture MOC]]
- [[Workspace Layout]]
