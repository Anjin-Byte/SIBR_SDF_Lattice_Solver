---
title: Feature Gating
tags: [pattern, cargo, configuration]
---

# Feature Gating

> **Rule:** Cargo features represent **real capability slices**, not random internal toggles.

## Good features

- `preview` — pulls in [[preview crate]] and any windowing deps
- `trace` — wires up `tracing` exporters / spans for production diagnostics
- `profiling` — enables backend profilers (e.g. `wgpu-profiler`, GPU timestamps)
- `shader-hot-reload` — dev-only, watches WGSL files and rebuilds pipelines

Each one corresponds to a thing a user or developer would consciously turn on.

## Bad features (avoid)

- `use-fast-path` — internal toggles that should just be conditional code
- `experimental` — too vague; what does enabling it actually do?
- `legacy` — usually a sign the old code should just be deleted
- one feature per minor switch — feature soup that nobody understands

## Default features

Default to **off** for anything that:
- adds heavy dependencies (windowing, GUI)
- changes runtime cost meaningfully
- is only useful in a specific environment (dev, CI, profiling)

Default to **on** for things that are part of the normal expected experience.

## Cross-crate feature plumbing

Feature in [[cli crate]] that turns on preview:

```toml
# crates/cli/Cargo.toml
[features]
default = []
preview = ["dep:preview"]

[dependencies]
preview = { path = "../preview", optional = true }
```

In `cli` source:

```rust
#[cfg(feature = "preview")]
mod preview_mode;
```

If you go with [[Binary Strategy]] Option B, you may not need a `preview` feature at all — the viewer is just a separate binary.

## Anti-pattern: testing every feature combination

If `N` features compose freely, you have `2^N` build configurations. Keep features orthogonal, well-scoped, and few. CI should test the meaningful combinations, not the cartesian product.

## Related

- [[Cargo Workspace Configuration]]
- [[Binary Strategy]]
- [[Headless First]]
