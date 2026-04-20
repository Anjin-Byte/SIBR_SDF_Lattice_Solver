---
title: preview crate
tags: [crate, architecture, optional]
crate-path: crates/preview
status: optional
---

# `preview` crate

**Add only when an interactive surface is actually needed.** The CLI must work without it.

## Responsibilities (what goes IN)

- `winit` window and event loop
- Surface lifecycle (create / configure / resize / destroy)
- Redraw loop
- Viewport logic
- Interactive inspection UI (if any)

## What stays OUT

- Compute logic (lives in [[gpu crate]])
- SDF math and pipeline types (live in [[sdf crate]] and [[lattice-gen crate]])
- File I/O / argument parsing (lives in [[cli crate]])

## Why this is its own crate (or feature)

Windowing is platform-y, drags in a heavy event-loop dependency, and changes on a different cadence than the GPU compute path. Isolating it means:

- Server/CI builds of the CLI can omit it entirely.
- The redraw loop and event-handling code don't pollute the headless code path.
- If the viewer grows (panels, inspectors, tools), it can evolve independently.

## winit lifecycle expectations

Modern winit guidance — **follow it from day one**, not after porting from old patterns:

- Windows are created from inside the active event loop, typically via `ActiveEventLoop::create_window` in `resumed()`.
- Render in response to **redraw events**, not by driving rendering from `about_to_wait`.

Do not architect this crate like an old immediate-mode game loop from 2017 blog posts. Those patterns are nostalgic but not portable truth.

See [[winit Lifecycle]].

## Integration shape

`preview` depends on [[gpu crate]] (and likely [[sdf crate]] / [[lattice-gen crate]] for shared types). It does **not** reach into wgpu directly — it asks the gpu crate for what it needs (surface formats, texture handles to present, etc.) through a small additional API surface.

## Related

- [[Headless First]]
- [[Binary Strategy]]
- [[winit]]
- [[Feature Gating]]
