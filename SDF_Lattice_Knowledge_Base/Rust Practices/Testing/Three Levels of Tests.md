---
title: Three Levels of Tests
tags: [rust, testing, principle]
---

# Three Levels of Tests

Rust distinguishes **unit tests**, **integration tests**, and **doc tests**. Strong projects use all three because they catch different classes of mistake.

## The split

### Unit tests
- Live alongside the code they test, in `#[cfg(test)] mod tests { use super::*; }`.
- Have access to private items in the same module.
- Target small pure functions, type invariants, edge cases.
- Should be fast and independent.

### Integration tests
- Live in the crate's `tests/` directory.
- See only the **public API** of the crate — exactly what external users see.
- Target end-to-end flows through the public surface.
- Catch API drift and accidental breakage of public behavior.

### Doc tests
- Live inside `///` doc comments as fenced ` ```rust ... ``` ` blocks.
- Compiled and run by `cargo test` automatically.
- Keep documentation **honest** — examples can't lie if they're executed.

## Why all three

| Test kind | Catches |
|---|---|
| Unit | Logic bugs, off-by-one, invariant violations, branch coverage |
| Integration | API regressions, wiring mistakes, real-world flow bugs |
| Doc | Documentation drifting from reality |

Each level catches things the others can't. Skipping a level means a class of bug has nowhere to be caught.

## Where to put what

A pure function on a domain type? **Unit test** it.

A scenario like "load this config file, run the solver, write this output"? **Integration test** in `tests/`.

A worked example in a public type's docs? **Doc test**, so it stays compiling.

## Doc tests in practice

```rust
/// Constructs a [`Port`] from a non-zero `u16`.
///
/// # Examples
///
/// ```
/// use mycrate::Port;
/// let p = Port::new(8080).expect("nonzero");
/// assert_eq!(p.get(), 8080);
/// ```
///
/// Returns `None` for the zero value:
///
/// ```
/// use mycrate::Port;
/// assert!(Port::new(0).is_none());
/// ```
pub fn new(value: u16) -> Option<Self> { /* ... */ }
```

Both blocks compile and run as part of `cargo test`. If the API changes and the example breaks, CI catches it before the docs become a liar.

## Doc tests that intentionally don't compile

Use code-block annotations:

- ` ```ignore ` — present in docs, not compiled.
- ` ```no_run ` — compiled but not executed (for examples that need a real GPU, network, etc.).
- ` ```compile_fail ` — must fail to compile (useful for documenting type-system safety).

Use these sparingly — the whole point is that the example matches reality.

## How this layers with the project's architecture

- **`core` crate**: lots of unit tests + doc tests. Easy because the crate is pure.
- **`gpu` crate**: unit tests for layout/planning logic; integration tests for end-to-end pipeline runs (some `no_run` doc examples for snippets that need an adapter).
- **`cli` crate**: integration tests in `tests/` exercising the binary's argument parsing and orchestration.
- **Cross-crate** integration tests: at the workspace root in `/tests/` for flows that span crates.

## Related

- [[Edge Cases and Properties]]
- [[Tests Without Heroics]]
- [[Documentation as Truth]]
- [[CI Quality Bar]]
