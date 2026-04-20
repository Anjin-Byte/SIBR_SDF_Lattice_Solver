---
title: Clippy as Discipline
tags: [rust, tooling, clippy, quality]
---

# Clippy as Discipline

Clippy is a large collection of lints organized by category. Use it **aggressively**, configure it **intentionally**, and run it on **every build**.

## What Clippy is for

It catches:

- Idiomatic mistakes (clippy::style)
- Likely bugs (clippy::correctness)
- Suspicious patterns (clippy::suspicious)
- Performance hazards (clippy::perf)
- Complexity smells (clippy::complexity)
- Pedantic improvements you might want (clippy::pedantic, off by default)

Each lint is opt-in or opt-out. The defaults are well-tuned; the project should layer specific choices on top.

## Configuration in this workspace

Clippy is configured in the **root** `Cargo.toml` via `[workspace.lints.clippy]` and inherited by member crates. See [[Workspace Lints and Profiles]].

A reasonable baseline:

```toml
[workspace.lints.clippy]
# Catch real footguns
unwrap_used = "warn"
expect_used = "warn"
panic = "warn"

# Pedantry — opt in with low priority so individual lints can override
pedantic = { level = "warn", priority = -1 }

# Common pedantic lints to relax (often noisy without value)
module_name_repetitions = "allow"
must_use_candidate = "allow"
missing_errors_doc = "allow"  # consider promoting to "warn" for the sdf/lattice-gen/gpu crates
```

Tune over time as the project develops opinions.

## "Run clippy" is not enough

Make Clippy part of the build discipline, not a thing you run when bored. See [[CI Quality Bar]] for how it fits in CI.

Useful command:

```sh
cargo clippy --workspace --all-targets --all-features -- -D warnings
```

- `--workspace` — all crates.
- `--all-targets` — tests, examples, benches too. Catches issues that production-only checks miss.
- `--all-features` — all combinations build cleanly. (For larger workspaces with many features, scope this to meaningful sets.)
- `-D warnings` — turn lint warnings into errors for the run. Forces them to be addressed.

## Per-crate overrides, when justified

Some crates have legitimately different needs:

```toml
# xtask/Cargo.toml — pragmatism beats purity here
[lints.clippy]
unwrap_used = "allow"
```

`xtask` is build automation, not production code. Loosening lints there is fine. Loosening them in `core` or `gpu` should require a written reason.

## What Clippy doesn't catch

Clippy is not a substitute for:

- **Tests** — see [[Three Levels of Tests]].
- **Type-driven design** — Clippy can't redesign your types for you. See [[Type System as Design Tool]].
- **Code review** — humans still catch architectural drift, missing test cases, and bad naming.

Clippy is the floor of code quality, not the ceiling.

## Related

- [[CI Quality Bar]]
- [[Workspace Lints and Profiles]]
- [[Documentation as Truth]]
