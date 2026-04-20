---
title: Workspace Lints and Profiles
tags: [pattern, cargo, configuration, quality]
---

# Workspace Lints and Profiles

> **Rule:** profiles and lints live in the **root** manifest only. Member crates inherit.

## Profiles

`[profile.*]` settings are **only** recognized in the root manifest of a workspace. If you set them in a member crate's `Cargo.toml`, Cargo silently ignores them.

GPU-aware defaults:

```toml
[profile.dev]
opt-level = 1
debug = 2

[profile.release]
lto = "thin"
codegen-units = 1
debug = 1
```

### Why these settings

- **`profile.dev` opt-level = 1**: GPU code at `opt-level = 0` is misleading slow. Math-heavy CPU paths (and the [[CPU Reference Path]]) feel broken when really they're just unoptimized. `opt-level = 1` keeps debug-info friendliness while not lying about performance.
- **`profile.release` lto = "thin"**: meaningful inlining across crate boundaries, much faster to link than full LTO.
- **`profile.release` codegen-units = 1**: best optimization at the cost of compile time. Acceptable for shipping builds.
- **`profile.release` debug = 1**: keeps line tables. When something crashes in production or you want to attach a profiler, you'll be glad. Strip later if binary size matters.

## Lints

Cargo supports `[workspace.lints]` and member crates inherit with `lints.workspace = true`.

```toml
# root Cargo.toml
[workspace.lints.rust]
unsafe_op_in_unsafe_fn = "warn"
missing_docs = "warn"

[workspace.lints.clippy]
unwrap_used = "warn"
expect_used = "warn"
pedantic = { level = "warn", priority = -1 }
```

```toml
# each member Cargo.toml
[lints]
workspace = true
```

### Why centralize

For a GPU CLI app, consistent rules across crates keep the workspace feeling like one system instead of four roommates sharing a fridge. In particular:

- **`unsafe` usage** — relevant for any FFI or raw-pointer GPU code
- **`unwrap` / `expect` in non-test code** — surfaces footguns early
- **`missing_docs` on public APIs** — especially valuable for [[gpu crate]]'s narrow public surface
- **Clippy pedantry** in [[sdf crate]], [[lattice-gen crate]], and [[gpu crate]] — these are the slow-changing crates where rigor pays off

### Per-crate overrides

A member crate can downgrade or allow specific lints when needed:

```toml
[lints.clippy]
unwrap_used = "allow"  # e.g. in xtask, where pragmatism beats purity
```

But such overrides should be rare and intentional.

## Related

- [[Cargo Workspace Configuration]]
- [[Workspace Layout]]
