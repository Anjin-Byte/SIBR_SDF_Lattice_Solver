---
title: Small Functions, Narrow Contracts
tags: [rust, functions, principle, design]
---

# Small Functions, Narrow Contracts

> Testability improves when each function has **one job and a crisp contract**. Stability improves because fewer things can go wrong in one place.

## Smell test

For any function you write, ask:

1. **Can it be named with a verb phrase?** (`load_config`, `validate_mesh`, `submit_compute_pass`)
2. **Can I describe its preconditions in one sentence?** ("Takes a non-empty buffer of vertices and returns the bounding box.")
3. **Can I test it with 3–6 representative cases?** (typical, empty, boundary, malformed)

If any answer is "no," the function is probably doing too much — split it.

## Why Rust rewards this style especially

Ownership and borrowing become much easier to reason about when scopes are small. A 200-line function juggling multiple `&mut` references will fight the borrow checker constantly. Three 60-line functions, each with a clear input and output, almost never do.

The borrow checker's frustration is, again, [[Type System as Design Tool|design information]] — a sign that responsibilities are tangled.

## Contract anatomy

A "narrow contract" means stating, for each function:

- **What goes in** — types, valid ranges (often via [[Newtypes and Domain Types|newtypes]]), required state.
- **What comes out** — type, possible variants for `Result` errors.
- **What the function promises** — observable effects on inputs, side effects on the world.
- **What the function does NOT promise** — order of operations, performance characteristics if irrelevant.

Most of this is encoded in the signature. The doc comment fills in the rest. See [[Documentation as Truth]].

## Splitting heuristics

When a function is too big, useful axes to split along:

- **Validation vs operation.** Split "check the inputs" from "do the work."
- **Pure vs effectful.** Pull the pure transformation out of the I/O wrapper. See [[Pure Core Effectful Edges]].
- **Per-phase.** If the function has clear phases ("parse → analyze → emit"), each phase can be its own function.
- **Error-prone subroutine.** If one inner step fails in interesting ways, lift it out so its tests can target it directly.

## What "small" means

Not a hard line count, but a useful heuristic: if you can hold the whole function in your head while reading it, it's small enough. If you have to scroll back up to remember what a variable was, it isn't.

For Rust, ~20–60 lines is a common sweet spot for non-trivial logic. Trivial helpers can be much smaller. Some legitimately complex routines (parsers, schedulers) are larger — but those should still decompose into named, contracted helpers internally.

## Anti-pattern: the "kitchen sink" function

Signs:
- More than ~5 parameters, especially if some are bool flags.
- A name like `do_thing_with_options` or `handle_request`.
- Comments inside dividing the body into sections.
- A `match` on a "mode" parameter that switches between fundamentally different operations.

Each of these is a hint to extract.

## Related

- [[Push Correctness Left]]
- [[Pure Core Effectful Edges]]
- [[Tests Without Heroics]]
- [[Boring Data Layouts]]
