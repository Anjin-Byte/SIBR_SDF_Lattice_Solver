---
title: winit Lifecycle
tags: [pattern, winit, preview]
---

# winit Lifecycle

The modern winit event-loop expectations. **Build [[preview crate]] to these from day one** — don't port from old patterns.

## The two key rules

### 1. Create windows from inside the active event loop

Don't construct a `Window` ahead of time and hand it to the event loop. Instead, create it from `ActiveEventLoop::create_window`, typically in the `resumed()` handler.

This matters because window/surface creation needs to happen in the right platform context (especially on mobile and on suspended/resumed application lifecycle), and the active event loop is the only place that holds it.

### 2. Render in response to redraw events

Drive rendering from `WindowEvent::RedrawRequested`, not by stuffing a render call at the bottom of `about_to_wait`. Request a redraw (`window.request_redraw()`) when you actually need one — animation, input, or external state change.

## Why old blog posts mislead

Pre-`ApplicationHandler` winit code (and many tutorials still floating around) shows patterns like:

- creating the window before the event loop runs
- rendering every iteration in `EventsCleared` / `MainEventsCleared` / `about_to_wait`
- treating the event loop as a tight spin

These patterns do not portably represent how the platforms actually want apps to behave. They will misbehave on suspend/resume, on Wayland, and on platforms with stricter lifecycle semantics.

## Architectural consequence for [[preview crate]]

- The preview crate's main type implements winit's `ApplicationHandler` (or equivalent).
- Window and surface creation happen in `resumed()`.
- The render path is in `window_event(.., RedrawRequested)`.
- Compute work is **not** initiated from the event loop directly — it's requested from [[gpu crate]]'s `GpuEngine` and the result is presented when ready.

## Don't do

- Don't drive a 60Hz render loop from `about_to_wait`.
- Don't create the surface before the event loop is active.
- Don't block the event loop on long GPU work — submit and continue.

## Related

- [[preview crate]]
- [[winit]]
- [[Headless First]]
