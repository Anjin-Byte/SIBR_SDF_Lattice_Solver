---
title: Push Correctness Left
tags: [rust, principle, foundational]
---

# Push Correctness Left

The meta-principle behind every other Rust practice in this knowledge base.

> **Push correctness as far left as possible** — into types, construction, ownership, module boundaries, and tests — so fewer bugs are left to runtime chance.

## "Left" means earlier in the bug's lifecycle

Imagine the timeline of how a defect becomes user-visible:

1. **Type system** rejects it (best — bug never compiles)
2. **Constructor / validation** rejects it at object birth
3. **Module boundary** rejects it at API entry
4. **Unit test** catches it before merge
5. **Integration test** catches it before release
6. **Production** catches it (worst — bug becomes incident)

The further left you push the catch, the cheaper, faster, and louder the failure. Rust gives unusually good tools at the leftmost three slots — the type system, construction patterns, and module privacy — and we should use them.

## Why this principle leads everything

Every other practice in [[Rust Practices MOC]] is a corollary:

- [[Make Invalid States Unrepresentable]] — left-shift to slot 1
- [[Checked Constructors and Builders]] — left-shift to slot 2
- [[Domain Errors at Boundaries]] — left-shift to slot 3
- [[Three Levels of Tests]] — slots 4 and 5
- [[Clippy as Discipline]] — automated nagging that pulls slot 6 issues back to slot 4

## When this principle is uncomfortable

Pushing left often costs *more* code up front: a newtype instead of a `u32`, a builder instead of a struct literal, an enum instead of two booleans. That cost is real but local. The payoff is global — fewer surprises across every caller, every refactor, every future contributor.

The wrong reaction to "this newtype is verbose" is to delete the newtype. The right reaction is to ask whether the *use sites* could be cleaner.

## What this principle is **not**

It is not "ban runtime checks." Some invariants can only be checked at runtime (network input, file contents, GPU validation). For those, push the check **once**, at the boundary, and produce a validated type — then never re-check downstream. See [[State Transition Types]].

## Related

- [[Type System as Design Tool]]
- [[Predictable APIs]]
- [[Rust Practices MOC]]
