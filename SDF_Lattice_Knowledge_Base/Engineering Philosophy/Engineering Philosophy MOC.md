---
title: Engineering Philosophy MOC
tags: [moc, philosophy]
---

# Engineering Philosophy — Map of Content

How we **reason about, design, and test** software in this project. The working method by which the workspace architecture and Rust practices are applied.

## The thesis (one sentence)

> Write programs so their correctness is **explainable to a human**, then **attack them with adversarial tests** aimed at the strange, brittle, neglected corners that ordinary usage will not exercise.

This is the project's working method. It is in the spirit of Donald Knuth's approach to correctness, testability, and long-term reliability.

## The eight principles

1. [[Human Understanding First]] — narratable design; expose invariants
2. [[Reason About Correctness]] — argue for it; don't hope for it
3. [[Normal Usage Is Not Testing]] — typical use only exercises the obvious paths
4. [[Tests Should Make Programs Fail]] — adversarial testing finds important bugs
5. [[Sharp Oracles]] — exact, decisive correctness criteria; never "looks right"
6. [[Regression Discipline]] — every bug strengthens the suite, permanently
7. [[Explicit Not Magical]] — clarity over cleverness; debuggable by inspection
8. [[Reasoning Beside Code]] — produce the surrounding argument, not just code

## Output format

When designing or implementing non-trivial work, follow the structure documented in [[Output Format]]:

1. **Correctness model** — what the system must guarantee
2. **Invariants and assumptions** — the load-bearing claims
3. **Design and implementation** — the actual code
4. **Test strategy** — split into ordinary, edge, adversarial, and regression tests, each annotated with the failure class it targets

## How this composes with the other sections

| Section | Answers |
|---|---|
| [[Workspace Architecture MOC]] | **What does the system look like?** Crate layout, boundaries, deployable units. |
| [[Rust Practices MOC]] | **What does good Rust code look like?** Idioms, type-driven design, error handling, test levels. |
| **Engineering Philosophy** (this section) | **How do we *reason and test* toward correctness?** The attitude with which the other two are applied. |

The architecture defines *what good looks like in the large*. The Rust practices define *what good looks like at the line level*. The engineering philosophy defines *how we know we've achieved either*.
