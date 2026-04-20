---
title: Tests Should Make Programs Fail
tags: [philosophy, principle, testing]
---

# Tests Should Make Programs Fail

> The purpose of testing is to make the program fail until the important bugs are found.

This is the offensive posture. Tests are not validation theatre — they are an adversary that wants the program to break, so that the bugs it finds are removed before users find them.

## What this looks like in practice

Create tests that are **intentionally unlike ordinary usage**. Specifically:

- **Combine edge cases in unpleasant ways.** Empty input AND maximum-size flag. First-byte-corrupt AND otherwise valid. Zero AND negative AND NaN.
- **Exercise limits.** Off-by-one is the eternal classic, but also: integer overflow, capacity boundaries, recursion depth, allocation failure, file size near filesystem limits, GPU buffer size near hardware limits.
- **Feed malformed inputs.** Truncated headers, null bytes in the middle of strings, unicode normalization differences, mixed line endings, BOM in unexpected places.
- **Test boundary values.** `0`, `1`, `MAX-1`, `MAX`, `MAX+1` (where applicable). The interesting bugs cluster here.
- **Provoke degenerate states.** A graph with one node. A "polygon" with two vertices. A texture that is 1×1 or 0×N.
- **Pass contradictory options.** `--quiet --verbose`. `--input file --stdin`. `--output -` with `--format binary`. Either reject coherently or pick deterministically — but never silently choose one and discard the other.
- **Use unusual orderings.** Operations the API doesn't forbid but didn't anticipate. Calling `flush` before `write`. Iterating concurrently with mutation (in languages where that's even a question).
- **Repeat operations.** Initialize twice. Close twice. Submit the same job back-to-back. Idempotence violations are common and silent.
- **Exercise empty cases.** Most code paths assume nonempty input somewhere. Find the assumption.
- **Exercise maximal cases.** Most code paths assume "reasonable" size somewhere. Find that assumption.
- **Attempt invalid transitions.** Use the type system if you can ([[State Transition Types]]); test them at runtime if you can't.

## The slogan

> **Prefer tests that would embarrass a weak design.**

A test passes either because:

(a) the design is robust enough to handle that case, or
(b) the test was too gentle to expose the weakness.

Aim for (a) by writing tests that would expose (b). Tests that no one would feel embarrassed by failing are not earning their place in the suite.

## The relationship between this principle and design

Adversarial tests **inform design**. When a test reveals an awkward corner, the answer is often not just "add a guard." The answer is sometimes:

- This case should be made unrepresentable by the type system.
- This combination of options should not exist as a configuration.
- This API should not expose the operation in this state.

In other words: failed adversarial tests sometimes point toward [[Make Invalid States Unrepresentable]] rather than toward more runtime validation. Use them as design feedback, not just as bug-finding tools.

## The relationship between this principle and oracles

A test that "exercises" a case but cannot say whether the result is right is a **smoke test**, not an adversarial one. The adversarial value comes from a **sharp verdict**: this output is wrong, this state is invalid, this invariant was broken. See [[Sharp Oracles]].

## When the suite has caught up

A useful sign that the suite is doing its job: **adding new adversarial tests stops finding new bugs** for a while. That doesn't mean stop. It means the suite has tightened around the current design. The next round of bugs will arrive with the next round of changes.

## Related

- [[Normal Usage Is Not Testing]]
- [[Sharp Oracles]]
- [[Edge Cases and Properties]]
- [[Regression Discipline]]
