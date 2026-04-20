---
title: Workspace Architecture MOC
tags: [moc, architecture]
---

# Workspace Architecture — Map of Content

Entry point for understanding **why** the workspace is shaped the way it is, and **where** to look for any specific architectural decision.

## The thesis

> A wgpu CLI app should be modeled as a **native systems app with a portable GPU layer**, not as a web-app crate layout. Concerns evolve at different speeds — CLI UX changes often, GPU code changes carefully, domain logic must stay testable without a GPU, and preview/windowing must be optional.

## Foundational reading order

1. [[Core Principles]] — the four rules
2. [[Workspace Layout]] — directory structure
3. [[Cargo Workspace Configuration]] — root manifest
4. [[Architecture/Crates/]] — per-crate responsibilities
5. [[Architecture/Patterns/]] — recurring design patterns

## Crates

- [[sdf crate]] — SDF math library with type-level precision tracking
- [[lattice-gen crate]] — unit cell topologies + Woodward pipeline, CPU reference
- [[gpu crate]] — wgpu backend, narrow public API (Phase 2)
- [[cli crate]] — argument parsing, orchestration, I/O (Phase 2)
- [[preview crate]] — optional winit + surface (Phase 3)
- [[xtask crate]] — build/automation tasks

## Patterns

- [[GPU Boundary]] — keep wgpu types out of upstream crates
- [[Headless First]] — windowing is opt-in, never required
- [[CPU Reference Path]] — correctness oracle for GPU output
- [[Binary Strategy]] — one binary vs. two
- [[Feature Gating]] — when to use Cargo features
- [[Shared Dependencies]] — what belongs in `workspace.dependencies`
- [[Workspace Lints and Profiles]] — root-only configuration
- [[winit Lifecycle]] — modern event-loop expectations

## The condensed rules

1. Core logic must survive without a GPU.
2. wgpu lives behind a narrow crate boundary.
3. Windowing is optional and isolated.
4. Workspace root owns versions, profiles, and lints.
5. Default workspace commands target the CLI, not every crate.
6. Crates reflect responsibilities, not aesthetics.
7. Always keep a CPU reference path.
