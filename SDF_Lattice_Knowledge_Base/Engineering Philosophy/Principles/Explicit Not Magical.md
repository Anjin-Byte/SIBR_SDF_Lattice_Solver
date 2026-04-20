---
title: Explicit, Not Magical
tags: [philosophy, principle, design]
---

# Explicit, Not Magical

> Prefer **clear control flow**, **stable interfaces**, and **obvious data movement**. Avoid cleverness that obscures reasoning. Make debugging and inspection straightforward.

## What "explicit" means

- **Control flow you can follow with your eyes.** A reader should be able to trace a value from input to output by reading top-to-bottom, not by chasing through a web of macros, traits, and reflection.
- **Stable interfaces.** A function's signature should imply how it's called, what it requires, and what it returns. Surprising side effects, hidden global lookups, and "do something different based on context" behavior are anti-features.
- **Obvious data movement.** Where data is constructed, where it's transformed, where it's consumed. If a struct's contents change, it should be clear *which function* changed them.
- **Inspection-friendly.** When something goes wrong, a `Debug` print, a logged span, or a stepped debugger should reveal what's happening. State that requires a custom tool to inspect is state that hides bugs.

## Why "magical" is a problem

Magical code — clever metaprogramming, deep generic chains, action-at-a-distance via globals or interior mutability, behavior that depends on attributes the compiler interprets — fails on three axes:

1. **Reading.** The reasoning is hidden. The reader cannot see what the program does without already understanding the magic.
2. **Debugging.** When the magic misbehaves, you cannot step through it; you can only stare at the output and guess.
3. **Refactoring.** Touching anything risks invalidating an invariant nobody documented because the magic "obviously" maintained it.

"Magic" is fine to *use* (libraries are full of it). It is **not** fine to *introduce* in your own code unless the alternative is meaningfully worse.

## Concrete heuristics

- **Prefer concrete types over generics** when only one type is ever used. See [[Traits as Seams]].
- **Prefer functions over macros** when the macro doesn't earn its keep. Macros that exist to save typing are debt.
- **Prefer explicit collaborators over global lookups.** A function that takes a `&Clock` is reasonable; one that consults a `static CLOCK` is a hidden coupling.
- **Prefer one obvious way over three clever ones.** A codebase with three different "configuration" patterns has none.
- **Prefer Debug-printable state.** If a type is opaque, derive `Debug` thoughtfully or implement a useful manual one. A field whose value cannot be inspected is one you can't debug.
- **Prefer `match` over deep `if let`/`Option` chains** when the intent is "handle each case." The match makes the case set visible.

## When "explicit" loses

A few legitimate cases where the magical option wins:

- **Derives.** `#[derive(Debug, Clone, PartialEq)]` is magic, but it is **shared, well-known magic** that every Rust reader understands. Don't reimplement these by hand for purity points.
- **Well-established framework conventions.** Serde, clap, thiserror — these are magical *and* universal. Reinventing them by hand is worse for readers, not better.
- **Genuine performance-critical hotspots** where the explicit version measurably loses. Measure first.

The principle is "explicit, not magical" — not "no abstraction at all." Use abstractions whose magic is well-documented and widely understood.

## Composes with [[Human Understanding First]]

This principle is the *defensive* sibling of "human understanding first." That principle says: structure code so a reader can follow it. This one says: do not undermine that with cleverness in the interest of saving lines.

The two together: **build code a reader can read, and don't sneak in things the reader can't see**.

## Related

- [[Human Understanding First]]
- [[Predictable APIs]]
- [[Boring Data Layouts]]
- [[Traits as Seams]]
