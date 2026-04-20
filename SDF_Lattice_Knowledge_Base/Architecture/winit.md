---
title: winit
tags: [reference, dependency, windowing]
type: external-dependency
---

# winit

Cross-platform Rust window-creation and event-loop library. Companion to [[wgpu]] but **deliberately separate** from it.

## Mental model

> Window creation, OS event delivery, input — and **nothing GPU-specific**.

winit gives you a `Window` (which has a raw handle that wgpu can build a `Surface` from) and an `ApplicationHandler`-driven event loop. It is the platform layer, not a renderer.

## Where it lives in the workspace

**Only** in [[preview crate]]. The headless CLI ([[cli crate]]) does not depend on it. See [[Headless First]].

## Modern lifecycle (must follow)

See the dedicated [[winit Lifecycle]] note. Two key rules in summary:

1. **Create windows from inside the active event loop** — typically in `resumed()`, via `ActiveEventLoop::create_window`.
2. **Render in response to `RedrawRequested`** — not from `about_to_wait`.

Old blog-post patterns (pre-`ApplicationHandler` style, render-every-tick loops) will misbehave on platforms with strict suspend/resume semantics. Don't port from them.

## Why it matters that wgpu and winit are separate

The split is a feature, not friction. It lets:

- Headless GPU work happen with no windowing dependency (the SDF solver's primary mode).
- Preview/viewer code evolve independently of the compute backend.
- Test infrastructure run anywhere.

Fusing them in a single crate would undermine [[Headless First]].

## Related

- [[preview crate]]
- [[winit Lifecycle]]
- [[wgpu]]
- [[Headless First]]
