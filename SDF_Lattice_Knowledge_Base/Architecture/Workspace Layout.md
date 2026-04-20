---
title: Workspace Layout
tags: [architecture, layout]
---

# Workspace Layout

The recommended directory shape for the SDF Lattice Solver workspace.

## Full layout (all phases)

```
SIBR_SDF_Lattice_Solver/
  Cargo.toml              # virtual workspace manifest
  .cargo/
    config.toml
  crates/
    sdf/                  # SDF math library (no GPU, no lattice-specific logic)
      Cargo.toml
      src/
    lattice-gen/          # Unit cells + Woodward pipeline; CPU reference
      Cargo.toml
      src/
    gpu/                  # wgpu backend (Phase 2+)
      Cargo.toml
      src/
      shaders/            # *.wgsl lives next to pipeline code
    cli/                  # clap, file I/O, orchestration (Phase 2+)
      Cargo.toml
      src/
    preview/              # OPTIONAL: winit + surface (Phase 3)
      Cargo.toml
      src/
  tests/                  # cross-crate integration tests
  reference/              # primary-source reference documents (papers, etc.)
  xtask/                  # build/automation tasks
    Cargo.toml
    src/
```

## Why this shape

- **Virtual workspace root** — see [[Cargo Workspace Configuration]] for why and how.
- **`crates/*` glob** — every member sits in one place; trivial to add or move.
- **Shaders under `gpu/`** — see [[gpu crate]]. Pipeline code and shader assets evolve together. The [[sdf crate]] does *not* hold WGSL — instead, cross-crate differential tests in `tests/` verify CPU-Rust and WGSL implementations agree on the same inputs.
- **`xtask` at root** — see [[xtask crate]]. Sibling to `crates/` because it is build infrastructure, not a library.
- **`tests/` at root** — for cross-crate integration tests; per-crate unit tests live in their own `tests/` or inline `#[cfg(test)]`.
- **`reference/` at root** — primary-source PDFs, HTML exports, datasets. Not a crate. The [[Domain Knowledge MOC|Domain Knowledge notes]] cite from here.

## Phase 1 — starting state

```
workspace root
├── crates/sdf
├── crates/lattice-gen
├── xtask
└── reference/
```

This is what we build first. CPU-only. No GPU, no CLI, no windowing.

**Rationale:** the SDF math and the lattice pipeline are the *specification* of what the system computes. Establishing them correctly — with tests, a CPU reference path that runs anywhere, and a clean precision-tracked type system — gives the [[gpu crate]] a known-good oracle to target in Phase 2. See [[Push Correctness Left]] and [[CPU Reference Path]].

The cost: Phase 1 can only produce small meshes (minutes to hours at Woodward-scale). That is acceptable for establishing correctness before optimizing.

## Phase 2 — GPU acceleration + CLI

```
workspace root
├── crates/sdf
├── crates/lattice-gen
├── crates/gpu         ← added
├── crates/cli         ← added
├── xtask
└── reference/
```

The [[gpu crate]] accelerates `lattice-gen`'s pipeline; the [[cli crate]] exposes the result as a headless command-line tool. At this point the system is a shippable CLI.

## Phase 3 — optional interactive preview

```
workspace root
├── crates/sdf
├── crates/lattice-gen
├── crates/gpu
├── crates/cli
├── crates/preview     ← added when needed
├── xtask
└── reference/
```

See [[preview crate]] and [[Headless First]]. Only added if an interactive viewer is actually needed.

See also: [[Core Principles]], [[Binary Strategy]], [[Crates Index]]
