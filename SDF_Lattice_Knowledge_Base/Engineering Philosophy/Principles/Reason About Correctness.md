---
title: Reason About Correctness
tags: [philosophy, principle, correctness]
---

# Reason About Correctness

> Treat correctness as something to be **reasoned about**, not merely **hoped for**.

## What this means

- **Use types, interfaces, constructors, and module boundaries to make invalid states difficult or impossible.** Reasoning is easier when the set of representable states is closer to the set of legal states.
- **State the key invariants explicitly.** A property that holds because "the code happens to do that" is not maintained — it's lucky.
- **For each major component, explain *why* it should work before discussing *how* it is implemented.** The "why" is the correctness argument; the "how" is the artifact that satisfies it.
- **Where formal proof is impractical, provide an informal but disciplined argument.** "I tested it" is not an argument. "Here is the invariant and here is why every operation preserves it" is.

## The two halves of correctness work

### Halve 1 — make wrong states unrepresentable

This is the constructive half. See [[Make Invalid States Unrepresentable]] in [[Rust Practices MOC]] for the Rust-specific techniques: newtypes, validated constructors, type-state, exhaustive enums.

When the type system rejects an illegal state, you have not "added a check." You have **eliminated the failure mode from the design**. There is nothing left to reason about, because the broken case cannot exist.

### Halve 2 — argue for the cases the type system can't capture

Some invariants are not type-shaped:

- "After `flush()`, the on-disk state matches the in-memory state."
- "The neighborhood graph constructed from this lattice is connected."
- "GPU output equals CPU output within tolerance ε."

For these, write the argument. Even informally, even in a comment near the function, even in the design doc — but write it. The shape of a good argument:

1. **State the invariant.** Be specific.
2. **Identify the operations that could violate it.** Which functions touch the relevant state?
3. **Show that each one preserves the invariant.** Either by construction (the operation cannot do otherwise) or by maintenance (the operation actively re-establishes it).
4. **Note any assumptions** the argument depends on (clock monotonicity, single-writer, GPU validation enabled, etc.).

This is what [[Output Format]]'s "Invariants and Assumptions" section is for.

## Why "reasoning" beats "testing alone"

Tests are samples. Reasoning is universal.

A test that passes shows the program worked **for those inputs**. An invariant argument shows the program works **for every input the invariant covers**. The two are complementary — see [[Tests Should Make Programs Fail]] — but reasoning is what gets you general statements; tests probe the boundary.

The strongest defense uses both: invariants narrow the space of possible bugs; adversarial tests probe what slips through.

## Anti-patterns

- **"It works in the tests we have"** as a stand-in for an argument. The tests bound a region of input space; they don't cover the rest.
- **"The compiler accepts it"** as a stand-in for an argument. The compiler enforces type rules, not domain invariants.
- **Implicit invariants.** A property only one person knows is one resignation away from being lost. Externalize it.

## Related

- [[Human Understanding First]]
- [[Sharp Oracles]]
- [[Make Invalid States Unrepresentable]]
- [[Push Correctness Left]]
