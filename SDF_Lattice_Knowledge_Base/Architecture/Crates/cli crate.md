---
title: cli crate
tags: [crate, architecture, cli]
crate-path: crates/cli
---

# `cli` crate

The user-facing entry point. Thin orchestration over [[sdf crate]], [[lattice-gen crate]], and [[gpu crate]].

## Responsibilities (what goes IN)

- `clap` (or equivalent) argument parsing
- File loading and saving (STL / 3MF export, input job descriptions)
- Selecting backend or mode
- Diagnostics output (`tracing`, structured logs)
- Orchestration — wiring `lattice-gen` → `gpu` → output
- Exit codes

## What stays OUT

- Actual compute logic — that lives in [[sdf crate]], [[lattice-gen crate]], and [[gpu crate]]
- Window/event-loop code — that lives in [[preview crate]]

The CLI must not become the place where domain logic accidentally accumulates. If you find yourself reaching for math or pipeline code here, push it down.

## Default workspace target

The root `Cargo.toml` sets:

```toml
default-members = ["crates/cli"]
```

So `cargo run`, `cargo build`, and `cargo test` at the repo root operate on this crate by default. See [[Cargo Workspace Configuration]].

## Binary strategy

Two viable shapes — see [[Binary Strategy]]:

- **Option A:** one binary `myapp` with `--preview` flag (single distributable, preview is convenience)
- **Option B:** separate `myapp-cli` and `myapp-viewer` binaries (cleaner separation, headless build has zero windowing dependencies)

Option B ages better for many projects.

## Optional `preview` feature

If using Option A, gate the preview path behind a Cargo feature:

```toml
[features]
default = []
preview = ["dep:preview-crate"]
```

See [[Feature Gating]].

## Related

- [[Core Principles]]
- [[Headless First]]
- [[Binary Strategy]]
