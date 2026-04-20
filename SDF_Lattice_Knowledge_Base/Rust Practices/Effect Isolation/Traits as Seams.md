---
title: Traits as Seams
tags: [rust, idiom, traits, testability]
---

# Traits as Seams

Traits are great for **injecting behavior** in tests and **separating interfaces from implementations**. They are not great for sport.

## When a trait is the right move

Use a trait when you genuinely need a substitution point:

- **Clock abstraction** — production uses system time, tests inject a fixed instant.
- **Filesystem abstraction** — production uses `std::fs`, tests use an in-memory implementation.
- **Repository / service boundaries** — for swapping persistence backends.
- **CPU vs GPU engine boundary** — see [[CPU Reference Path]].
- **Parser / serializer boundaries** — for plugging in different formats.

In each case there is a real second implementation (test double, alternate backend, future replacement) that justifies the trait.

## When a trait is *not* the right move

- "I might want to swap this someday" — without a concrete second implementation in mind, this leads to abstractions that fit nothing well. Add the trait when the second implementation actually arrives.
- "Six layers of generic indirection to avoid naming a concrete type" — generic soup makes errors longer, compile times slower, and readers more lost. If only one type ever implements the trait, prefer the concrete type.
- "I feel architecturally guilty using a struct" — guilt is not a design constraint.

## Concrete vs trait object

Two ways to consume a trait:

```rust
// Generic — monomorphized, zero runtime cost, can balloon binary size
fn run<C: Clock>(clock: &C) { /* ... */ }

// Trait object — single instantiation, runtime dispatch, smaller binary
fn run(clock: &dyn Clock) { /* ... */ }
```

Pick **generic** when:
- Performance matters in a hot path.
- The trait has associated types or generic methods.
- You need different concrete types in different call sites.

Pick **trait object** when:
- The collection holds a heterogeneous set (`Vec<Box<dyn Renderer>>`).
- You don't want monomorphization explosion.
- The dynamic dispatch cost is in the noise.

For test seams in non-hot paths, trait objects are usually cleaner.

## Keep trait surface small

A trait should be the smallest surface a substitute needs to provide. Big "kitchen sink" traits make test doubles painful (you implement 12 methods to test one thing).

Split into:
- A minimal trait with what callers actually need.
- Free functions or extension traits for derived helpers.

## Related

- [[Pure Core Effectful Edges]]
- [[Tests Without Heroics]]
- [[CPU Reference Path]]
