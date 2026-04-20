---
title: xtask crate
tags: [crate, architecture, tooling]
crate-path: xtask
---

# `xtask` crate

The standard "grown-up Rust workspace" move. A binary crate at the workspace root that owns ad-hoc automation, written in Rust instead of shell.

## Why have one early

Without `xtask`, automation tends to crystallize as a layered geology of shell scripts, makefiles, and CI YAML — none of which run the same way locally. Cursed archaeological layers.

With `xtask`, automation is:
- typed
- cross-platform
- discoverable (`cargo xtask --help`)
- testable
- versioned with the code

## What goes in `xtask`

- Shader validation / build steps
- Snapshot generation for golden tests
- Benchmark orchestration
- Packaging and release artifacts
- Test data generation
- Backend smoke tests (e.g. spin up adapter, log capabilities)
- Anything you'd otherwise put in a `Makefile` or `scripts/foo.sh`

## What does NOT go in `xtask`

- Library code reused by [[sdf crate]] / [[lattice-gen crate]] / [[gpu crate]] / [[cli crate]] — that's a normal crate.
- Production CLI commands the user runs — that's [[cli crate]].

## Convention

Invoke via `cargo xtask <subcommand>`. Use `clap` for subcommand structure so adding tasks stays cheap.

`xtask` typically depends on [[sdf crate]] or [[lattice-gen crate]] for any shared types (e.g. fixture descriptions) but should **not** depend on `cli` or `preview`.

## Workspace membership

In the root `Cargo.toml`:

```toml
[workspace]
members = ["crates/*", "xtask"]
```

`xtask` sits as a sibling to `crates/` because it is build infrastructure, not a library shipped with the application.

## Related

- [[Workspace Layout]]
- [[Cargo Workspace Configuration]]
