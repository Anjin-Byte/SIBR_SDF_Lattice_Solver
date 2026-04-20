---
title: gpu crate
tags: [crate, architecture, gpu]
crate-path: crates/gpu
---

# `gpu` crate

All wgpu lives here. The rest of the workspace sees a small, app-level API; the backend ugliness stays internal.

## Responsibilities (what goes IN)

- Adapter and device acquisition
- Pipeline creation
- Bind groups
- Shader module loading (WGSL files under `crates/gpu/shaders/`)
- Upload/download utilities
- Command encoders
- Debug labeling
- Error-scope helpers
- A **narrow public API** like `GpuEngine::run(job) -> Output`

## Shaders live here

```
crates/gpu/shaders/*.wgsl
```

Pipeline code and shader assets evolve together — keeping them adjacent is the right default. Move shaders to a dedicated crate **only** if real reasons emerge (preprocessing, reflection, validation helpers, codegen across multiple consumers).

## Public API shape

**Do this:**

```rust
pub struct RunOptions { /* ... */ }
pub struct GpuEngine { /* ... */ }
pub struct RunResult { /* ... */ }

impl GpuEngine {
    pub async fn new(/* ... */) -> Result<Self>;
    pub async fn process_field(
        &self,
        input: &Field,
        opts: &RunOptions,
    ) -> Result<RunResult>;
}
```

**Not this:**

```rust
pub fn do_everything(device: &wgpu::Device, queue: &wgpu::Queue, /* ... */)
```

The latter leaks too much machinery upward and makes future refactors miserable. See [[GPU Boundary]].

## What this crate consumes

- SDF primitives and combinators from [[sdf crate]] — used to generate WGSL shaders whose behavior matches the CPU reference
- Pipeline types (`LatticeJob`, `Mesh`, state-transition markers) from [[lattice-gen crate]]
- `wgpu`, `pollster` (or another async runtime), `bytemuck`, etc.

## What this crate exposes

- High-level engine types
- Typed errors at the boundary (no `wgpu::Error` leaking)
- Optional capability flags (e.g. backend choice, debug mode)

## Related

- [[Core Principles]]
- [[GPU Boundary]]
- [[wgpu]]
- [[Feature Gating]]
