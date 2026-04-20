---
title: Human Understanding First
tags: [philosophy, principle]
---

# Human Understanding First

> Structure code so that a reader can understand the **reasoning** behind it. Make the design **narratable**.

## What this means in practice

- **Expose invariants, assumptions, preconditions, postconditions, state transitions.** A reader should not have to reverse-engineer them.
- **Prefer organization that supports explanation and correctness over organization that is merely convenient for the machine.** A clever flat layout that's hard to walk through is a worse choice than a more verbose one a reader can follow.
- **Each major component should have a paragraph explaining its job before its code.** Not what it does (the code shows that) — *why it exists*, *what contract it satisfies*, *what invariants it maintains*.
- **Names earn their length.** A six-character abbreviation that needs a comment to decode is worse than a twenty-character name that doesn't.

## Why this is the first principle

If a reader cannot understand the reasoning, every other principle fails:

- They cannot **audit the correctness argument** ([[Reason About Correctness]]).
- They cannot **identify the neglected paths** to attack with adversarial tests ([[Normal Usage Is Not Testing]]).
- They cannot **debug a failure** when one occurs ([[Explicit Not Magical]]).
- They cannot **strengthen the suite** when a bug is found ([[Regression Discipline]]).

Understanding is the substrate. Without it the rest is decorative.

## Practical tests

Read the code as if you were a teammate joining the project six months from now:

- Can you state, after reading the file, what **invariants** the type or function maintains?
- Can you point to **where** each invariant is established?
- Can you say **why** the design chose this shape rather than an alternative?
- Can you find the **state transitions** without grepping for `mut`?

If any answer is "no," the structure is hiding the reasoning, not expressing it.

## What this is NOT

- **It is not "more comments."** Most comments restate code or rot into lies. Comments are useful for the *why* the code can't show — invariants, non-obvious constraints, references to bugs. See [[Documentation as Truth]] for the discipline.
- **It is not novel formatting or decorative diagrams.** A clear struct definition with one paragraph above it usually beats a UML diagram in a separate file.
- **It is not exhaustive prose.** A trivial helper does not need a treatise; surface the reasoning in proportion to the stakes.

## How this composes with [[Rust Practices MOC]]

- [[Newtypes and Domain Types]] makes invariants legible at the type level, not just in prose.
- [[State Transition Types]] makes phases visible at the type level — the reader sees the lifecycle in the signatures.
- [[Boring Data Layouts]] is this principle expressed at the data-modeling level.

## Related

- [[Reason About Correctness]]
- [[Explicit Not Magical]]
- [[Reasoning Beside Code]]
- [[Predictable APIs]]
