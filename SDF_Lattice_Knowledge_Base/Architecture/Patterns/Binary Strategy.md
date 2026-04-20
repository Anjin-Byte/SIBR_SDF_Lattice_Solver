---
title: Binary Strategy
tags: [pattern, architecture, deployment]
---

# Binary Strategy

How to ship the application: one binary, or two?

## Option A — single CLI binary

```
myapp           # supports --preview
```

**Good when:**
- Preview is a convenience mode the same user occasionally wants.
- A single distributable simplifies install/docs.
- The dependency cost of dragging windowing into headless invocations is acceptable.

**Implementation:** [[cli crate]] gains a `preview` Cargo feature that pulls in [[preview crate]]. See [[Feature Gating]].

## Option B — separate binaries

```
myapp-cli       # headless, no winit dependency
myapp-viewer    # window, surface, interaction
```

**Good when:**
- Headless builds (server, CI) should have **zero** windowing dependencies.
- Deployment for batch/server use should be lean.
- The viewer grows richer (panels, inspectors, tools) and benefits from independent evolution.

**For this project,** Option B will likely age better. Server and CI use of the SDF solver should not require windowing libraries to be present.

## Recommendation for SDF Lattice Solver

Start with **headless-only** (no `preview` crate at all). When/if interactive preview becomes useful, prefer **Option B** — a separate `myapp-viewer` (or named appropriately) binary in its own crate that depends on [[gpu crate]], [[lattice-gen crate]], and [[sdf crate]].

## What both options share

- [[sdf crate]] and [[lattice-gen crate]] are consumed by both.
- [[gpu crate]] is consumed by both.
- Neither binary contains compute logic itself — both orchestrate.

## Related

- [[Headless First]]
- [[cli crate]]
- [[preview crate]]
- [[Feature Gating]]
