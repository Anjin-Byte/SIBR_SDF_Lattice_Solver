---
title: Headless First
tags: [pattern, architecture]
---

# Headless First

> **Rule:** the CLI runs to completion with no window, no surface, and no event loop. Preview is opt-in.

## Why this is non-negotiable

A CLI tool gets used in places where windows can't or shouldn't open:
- CI pipelines
- Headless servers
- Containerized batch jobs
- Automated regression suites
- Remote SSH sessions

If the headless path is an afterthought, all of the above become awkward — extra flags, X server stubs, suppressed exceptions, dummy event loops, etc.

## Architectural consequence

- [[gpu crate]] must support **off-screen rendering / compute** without a surface. Acquire an adapter and device without ever requesting a `Surface`.
- [[cli crate]] never imports `winit`.
- [[preview crate]] is either a separate binary (preferred) or behind a Cargo feature gate.

See [[Binary Strategy]] and [[Feature Gating]].

## Ecosystem alignment

The Rust GPU ecosystem already mirrors this split: `wgpu` handles the GPU API, `winit` handles windowing. They are deliberately separate libraries. Treating them that way in the workspace — instead of fusing them — keeps you on the well-trodden path.

See [[wgpu]] and [[winit]].

## What "preview" means here

- A window opens.
- A surface is created and presented to.
- The user can watch results visually, possibly interact (rotate, inspect, scrub a parameter).

It is **not**:
- Required for any compute output.
- A debugging substitute for proper logging or snapshot tests.

## Related

- [[Core Principles]]
- [[preview crate]]
- [[winit Lifecycle]]
