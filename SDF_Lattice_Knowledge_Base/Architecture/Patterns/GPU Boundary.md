---
title: GPU Boundary
tags: [pattern, architecture, gpu]
---

# GPU Boundary

> **Rule:** wgpu types must not appear in the public API of any crate other than [[gpu crate]].

## What this means in practice

Crates that depend on `gpu` see structures like:

```rust
pub struct GpuEngine { /* ... */ }
pub struct RunOptions { /* ... */ }
pub struct RunResult { /* ... */ }

impl GpuEngine {
    pub async fn process_field(
        &self,
        input: &Field,            // from `core`
        opts: &RunOptions,        // from `gpu`
    ) -> Result<RunResult>;       // from `gpu`
}
```

They do **not** see:

```rust
pub fn do_everything(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    /* ... */
)
```

## Why the boundary matters

1. **Refactor safety.** If wgpu changes its API (it has, repeatedly), only the gpu crate breaks.
2. **Backend swap-ability.** The codebase can swap to a different backend (or add a CPU-only fallback) without rewriting [[cli crate]].
3. **Testability.** Crates above the boundary can be tested without a GPU adapter.
4. **Mental load.** Readers of `cli` don't have to know wgpu to understand orchestration.

## Errors cross the boundary too

Don't leak `wgpu::Error` upward. Wrap in a typed error (`thiserror`) at the gpu crate edge:

```rust
#[derive(thiserror::Error, Debug)]
pub enum GpuError {
    #[error("no compatible adapter found")]
    NoAdapter,
    #[error("validation: {0}")]
    Validation(String),
    /* ... */
}
```

## Practical heuristic

When writing code in `cli` or `core`, if you reach for `use wgpu::*` — stop. The thing you need probably belongs in [[gpu crate]] as a higher-level operation.

## Related

- [[Core Principles]]
- [[gpu crate]]
- [[sdf crate]]
- [[lattice-gen crate]]
