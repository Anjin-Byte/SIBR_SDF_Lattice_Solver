---
title: Domain Errors at Boundaries
tags: [rust, error-handling, idiom]
---

# Domain Errors at Boundaries

> **Rule:** typed `enum` errors (via `thiserror`) at library boundaries; `anyhow`-style wrappers only at the *application* boundary.

## The pattern

```rust
#[derive(Debug, thiserror::Error)]
pub enum ParseConfigError {
    #[error("missing field: {0}")]
    MissingField(&'static str),
    #[error("invalid port: {0}")]
    InvalidPort(u16),
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
}
```

This gets you:
- **Testability** — assert on specific variants without string matching.
- **API clarity** — callers see exactly what can go wrong.
- **Future refactoring** — adding a new failure variant is a deliberate API change, not a silent string drift.
- **Pattern-matching at the call site** — callers can handle `MissingField` differently from `Io` cleanly.

## anyhow vs thiserror

- **`thiserror`** — define typed errors at *library* boundaries: per-crate, per-module, per-public-API. Use this for the `gpu` crate's public surface, the `core` crate's domain operations, etc.
- **`anyhow`** — wrap heterogeneous errors at the *application* boundary, typically in [[cli crate]]'s `main` or top-level orchestration. Convenient for "I just want to bubble a context-rich message to the user and exit."

The mistake to avoid: using `anyhow` deep inside library code. That throws away the structure the rest of the system needs to handle errors precisely.

## Conversion plumbing

`thiserror`'s `#[from]` lets you compose errors without manual `From` impls:

```rust
#[derive(Debug, thiserror::Error)]
pub enum LoadError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse: {0}")]
    Parse(#[from] ParseConfigError),
}
```

Now `?` can convert `std::io::Error` or `ParseConfigError` into `LoadError` without ceremony.

## Adding context

When converting between error types, add context the caller actually needs:

```rust
fn load(path: &Path) -> Result<Config, LoadError> {
    let bytes = std::fs::read(path)?;             // Io variant
    let config = parse(&bytes)?;                   // Parse variant
    Ok(config)
}
```

If the bare `std::io::Error` doesn't tell the caller *which path failed*, wrap it with `.map_err(|e| LoadError::IoAt(path.to_owned(), e))` or attach context via `anyhow::Context` (only at the application boundary).

## Variant naming

- Name variants for the **failure category**, not the underlying cause. `InvalidPort` is more useful than `IntegerError`.
- Keep variant payloads small and structural: the offending field name, the bad value, the path. Resist the urge to embed full backtraces (that's `anyhow`'s job).
- Reserve `Other(String)` as a last resort, not a first reach.

## Boundary, repeated

Each crate boundary is its own opportunity for typed errors. The `gpu` crate's public functions return `gpu::Error` enums — they do *not* leak `wgpu::Error` (see [[GPU Boundary]]). The `core` crate's public functions return `core::DomainError`. The CLI orchestrates these and finally collapses everything to `anyhow::Error` for the user.

## Related

- [[Result vs Panic]]
- [[Predictable APIs]]
- [[GPU Boundary]]
