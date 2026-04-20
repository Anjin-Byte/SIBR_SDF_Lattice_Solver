---
title: Core Principles
tags: [architecture, principles]
---

# Core Principles

Four principles drive every workspace decision. If a future change violates one of these, treat that as a smell and revisit the design.

## 1. Keep the GPU boundary narrow

Core algorithms should not know about wgpu types unless absolutely necessary.

**Bad shape:** application logic directly manipulates `Device`, `Queue`, `Buffer`, `Texture`, surface config, window state, CLI args, and file I/O all in one crate.

**Better shape:** split by responsibility:

- **sdf crate** — SDF primitives, combinators, operators; precision tracked at the type level
- **lattice-gen crate** — unit cell topologies; Woodward seven-stage pipeline; CPU reference meshing
- **gpu crate** — all wgpu device, pipeline, shader, buffer, and command encoding work
- **cli crate** — argument parsing, file I/O, logging, exit codes
- **optional preview crate** — window + event loop + surface presentation

That way, most logic can be tested without a GPU at all.

See: [[GPU Boundary]], [[sdf crate]], [[lattice-gen crate]], [[gpu crate]]

## 2. Treat headless as the default, preview as an add-on

A CLI app should be able to run with no window. If on-screen preview is needed later, bolt it on as a separate concern. This aligns with the ecosystem split: [[winit]] is a cross-platform window/event-loop library, distinct from wgpu itself.

See: [[Headless First]], [[preview crate]]

## 3. Make the workspace reflect deployable units

Rule of thumb: **if something can be built, tested, or shipped independently, it deserves its own crate.**

Do not split crates just to look sophisticated. A workspace with five meaningful crates is better than one with twelve decorative gourds.

See: [[Workspace Layout]]

## 4. Always keep a CPU reference path

For GPU-heavy CLI work, this is almost a law of nature.

Reasons:
- correctness oracle
- testability in CI without a GPU
- easier debugging
- deterministic comparisons

Workspace consequence:
- **core**: canonical behavior
- **gpu**: accelerated implementation

Saves you from the classic GPU debugging experience: *"it is very fast and very wrong."*

See: [[CPU Reference Path]]
