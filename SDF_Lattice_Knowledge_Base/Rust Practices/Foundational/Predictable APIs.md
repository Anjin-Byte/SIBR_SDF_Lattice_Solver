---
title: Predictable APIs
tags: [rust, principle, foundational, api-design]
---

# Predictable APIs

Optimize public APIs for **predictability**, not cleverness. Aligns with the Rust API Guidelines on consistency, naming, and dependability.

## What predictable means

A predictable API is one where a competent caller can guess the right thing without reading the implementation.

Concretely:

- **Names imply behavior.** `len()` is cheap. `count_unique()` might not be. If a method called `id()` does a database query, the name is lying.
- **Ownership behavior does not surprise.** A function called `set_name(&mut self, name: String)` taking ownership is fine; if it instead clones internally and the caller didn't expect that, performance goes sideways quietly.
- **Constructors and conversions are unsurprising.** `From<A> for B` should be lossless and obvious. `TryFrom` when failure is possible. `parse` for string-shaped input. Don't invent novel conversion verbs.
- **Methods don't hide expensive work.** If a getter triggers I/O or allocates, the name should hint at it (`fetch_*`, `load_*`, `compute_*`).
- **Defaults are explicit when ambiguity matters.** `Default::default()` is fine for true zero values. For domain types, prefer named constructors so the default isn't invisible.

## Why this is a stability concern

Stability is partly technical and partly social: callers should not need telepathy. APIs that surprise generate bug reports, defensive wrapper code, and cargo-culted workarounds — all of which become load-bearing.

A surprising API is also harder to refactor, because callers have built mental models around the surprises. "Fixing" the surprise becomes a breaking change.

## Naming heuristics

- Use the verb form the standard library uses: `iter`, `into_iter`, `as_str`, `to_string`, `clone`, `with_capacity`. Match those patterns.
- Reserve `unwrap_*`, `expect_*` for explicit panic-on-failure conversions.
- `is_*` returns `bool`. `has_*` returns `bool`. Don't use them for fallible getters.
- Avoid abbreviations unless they're domain-standard (`buf`, `ctx`, `cfg` are fine; `prc`, `mng`, `xfr` are not).

## Related

- [[Push Correctness Left]]
- [[Type System as Design Tool]]
- [[Documentation as Truth]]
- [[Domain Errors at Boundaries]]
