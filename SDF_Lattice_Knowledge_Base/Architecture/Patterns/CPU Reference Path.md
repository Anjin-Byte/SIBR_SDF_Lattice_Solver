---
title: CPU Reference Path
tags: [pattern, architecture, testing]
---

# CPU Reference Path

> **Rule:** for every algorithm implemented on the GPU, a CPU implementation exists in [[sdf crate]] or [[lattice-gen crate]] and is the reference oracle.

## Why

GPU debugging without an oracle has a notorious failure mode: *"it is very fast and very wrong."* You get plausible-looking output, you have no idea whether it's correct, and bisecting is miserable because the GPU pipeline is opaque.

A CPU reference solves this by giving you:

1. **Correctness oracle.** Compare GPU output to CPU output on the same input. Diffs above a small numerical tolerance are bugs.
2. **CI testability.** Tests of the algorithm itself run anywhere, including on machines without a GPU adapter.
3. **Easier debugging.** Standard debuggers and prints work. You can step through.
4. **Deterministic comparisons.** Floating-point order of operations is controlled.

## Workspace consequence

- **`core`** holds the canonical, deterministic, testable implementation.
- **`gpu`** holds the accelerated implementation.
- A test crate (or `tests/` at the workspace root) exercises both with shared fixtures.

This is the structural payoff for [[Core Principles]] #4.

## Tolerances

Don't compare GPU and CPU outputs with bitwise equality — different summation orders and FMA fusion will diverge harmlessly. Use:

- absolute tolerance for small magnitudes
- relative tolerance for larger ones
- structural equality for non-numeric outputs (indices, topology)

## When the CPU path can lag behind

Pragmatically, the CPU reference doesn't have to ship in the binary. It can be:

- behind a `#[cfg(test)]` feature
- in a sibling crate consumed only by tests
- a slower, more obviously correct version of the algorithm

What matters is that **for any GPU output, there exists a CPU computation that produces the expected answer**.

## Related

- [[Core Principles]]
- [[sdf crate]]
- [[lattice-gen crate]]
- [[gpu crate]]
