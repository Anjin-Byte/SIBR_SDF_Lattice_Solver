---
title: Shared Dependencies
tags: [pattern, cargo, dependencies]
---

# Shared Dependencies

What belongs in `[workspace.dependencies]` — and what doesn't.

## Belongs at workspace level

Pin these once in the root `Cargo.toml` so all crates use the same version:

- `wgpu`
- `winit`
- `serde` (with `derive`)
- `tracing`, `tracing-subscriber`
- `anyhow`, `thiserror`
- `clap` (with `derive`)
- `pollster` (or whichever async runtime is chosen)
- `bytemuck` if used for GPU upload structs

A version drift between crates on `wgpu` or `winit` would be painful — they ship breaking changes meaningfully often.

## Member crate inheritance

```toml
# crates/gpu/Cargo.toml
[dependencies]
wgpu = { workspace = true }
serde = { workspace = true, features = ["derive"] }
```

Note: features can be **added** at the member level on top of the workspace baseline.

## Does NOT belong everywhere

Just because something is in `[workspace.dependencies]` does not mean every crate should pull it in.

Critical examples:
- [[sdf crate]] and [[lattice-gen crate]] should **not** depend on `wgpu` or `winit`. That's the whole point of [[GPU Boundary]] and [[Headless First]].
- [[cli crate]] should **not** depend on `winit`. That's [[preview crate]]'s job.
- `xtask` should **not** depend on `wgpu` or `winit` unless a specific automation task genuinely needs them.

## Why discipline here matters

Each unnecessary dependency:
- slows clean builds
- expands the attack/maintenance surface
- couples crates that should evolve independently
- makes "what does this crate actually do" harder to read at a glance

A dependency that doesn't appear in a crate's `Cargo.toml` is one you cannot accidentally use.

## Related

- [[Cargo Workspace Configuration]]
- [[Core Principles]]
- [[GPU Boundary]]
