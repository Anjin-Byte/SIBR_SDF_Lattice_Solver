---
title: Pure Core, Effectful Edges
tags: [rust, principle, effects, testability]
---

# Pure Core, Effectful Edges

Isolate file I/O, env vars, clocks, randomness, networking, and GPU/device access at the **edges** of the system. Keep the **center** pure and deterministic.

## The shape

```
[ filesystem ] ─┐
[   network  ] ─┤
[    clock   ] ─┼─▶ adapters ─▶  PURE CORE  ─▶ adapters ─▶ [ filesystem ]
[    GPU     ] ─┤                                          [    stdout    ]
[ env / args ] ─┘
```

- `main` parses args, sets up logging, wires dependencies.
- `lib` (or core modules) holds the actual behavior — pure functions, deterministic transformations.
- Adapter modules at the edge talk to the OS, GPU, network, etc.

## Why this is the highest-value habit in Rust

Tests that depend on the filesystem, the network, the clock, or a GPU adapter are slow, flaky, hard to set up, and hard to reason about. Tests that depend on **pure functions** are fast, cheap, and boring — the best kind.

By moving logic out of `main` (and out of effect-laden modules) into pure code, you don't just make it testable. You make refactors safer, parallel execution trivial, and reasoning local.

## How this maps to this project's architecture

This principle is the **per-crate** expression of [[GPU Boundary]] and [[Headless First]]:

- [[sdf crate]] and [[lattice-gen crate]] are the pure center.
- [[gpu crate]] is the adapter that wraps wgpu.
- [[cli crate]] is the adapter that wraps argv, stdin/stdout, and the filesystem.
- [[preview crate]] is the adapter that wraps winit.

The same shape recurs **inside** each crate too. A crate is allowed to talk to one effect (the gpu crate talks to wgpu) but should still keep its internal logic — pipeline assembly, layout planning — as pure as possible, with the actual `Device::*` calls in a thin layer.

## Diagnostic questions

Use these in code review:

- Does this function read the clock or env var directly? Could it take the value as a parameter?
- Does this function open a file path, or could it accept a `Read` / `Write`?
- Does this function call `wgpu::*` directly, or could it produce a description that an adapter executes?
- Could I test this function with no setup at all?

The closer the answer is to "yes, no setup," the further left along [[Push Correctness Left]] you've pushed it.

## Related

- [[Traits as Seams]]
- [[Tests Without Heroics]]
- [[GPU Boundary]]
- [[Headless First]]
