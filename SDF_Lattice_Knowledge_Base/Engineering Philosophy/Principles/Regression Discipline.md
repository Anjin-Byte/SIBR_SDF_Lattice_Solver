---
title: Regression Discipline
tags: [philosophy, principle, testing, regressions]
---

# Regression Discipline

> Every discovered bug should produce or strengthen a test. Tricky past failures become **permanent members** of the suite. Treat bug history as engineering knowledge.

## The principle

When a bug is found:

1. **Reproduce it as a failing test first.** Before fixing the code. The failing test is the proof you understand the bug; without it, the "fix" is speculative.
2. **Fix the code.** Make the test pass.
3. **Leave the test in the suite forever.** Even — especially — if the bug is "embarrassing" or "obviously won't recur."
4. **Generalize when you can.** If the bug arose from a class of input (e.g., empty inputs, max-sized values, repeated operations), add tests for nearby members of that class. The bug you found is rarely alone.

## Why this matters more than it looks

Three reasons.

### 1. Bugs come back

Code is rewritten. Refactors happen. Optimizations replace working code with cleverer working code, except the cleverer version handles the empty case differently. A regression test catches this **the first time it happens**, not the third time it ships.

### 2. Bugs are sample points of weakness

Each bug found is evidence of a class of failure the design did not protect against. Recording the bug as a test makes the weakness **legible**. Future maintainers — and future tools — can see the shape of past failures and reason about whether new changes are likely to violate the same invariants.

### 3. Tests preserve hard-won knowledge

A bug found takes hours, sometimes days, to root-cause. The fix is often one line. The reasoning behind the fix lives only in the engineer's head, the commit message, and (if you're lucky) a comment.

A regression test makes that reasoning **executable**. Anyone who breaks the invariant in the future learns immediately, without needing to re-discover what the original engineer learned.

## What a good regression test looks like

- **Named for what it defends.** `truncated_header_returns_typed_error` is better than `bug_42` (even if the comment links to issue #42).
- **Annotated with the failure mode.** A short doc comment naming what would happen if the regression returned. This is what makes the test legible to a future reader.
- **Targets the property, not the byte-for-byte behavior.** A regression test for "off-by-one on empty input" should test the empty case generally, not just the specific input that revealed the bug.
- **Not deleted "because it duplicates an existing test."** It might. The duplication is the point — it documents that this case is load-bearing.

## What NOT to do

- **Don't fix the code without a failing test first.** You'll never know if the fix actually fixes the reported bug, vs. some adjacent thing.
- **Don't delete regression tests during refactors** because they "no longer match the structure." If they no longer compile, port them. The behavior they defend is still required.
- **Don't replace adversarial regression tests with weaker ones.** "We replaced the brittle test with a smoke test" usually means "we removed the part that caught the regression."
- **Don't bury regressions in a `#[ignore]` block** because they're flaky. Flakiness is its own bug; investigate root causes. An ignored test is a deleted test with extra paperwork.

## Bug history as a knowledge source

The accumulation of regression tests, taken together, is a **map of the system's historical weaknesses**. It tells future contributors:

- Which types of inputs have caused trouble.
- Which APIs have been misused.
- Which invariants were once silently violated.
- Which platform/environment edge cases bit us.

Reading the regression suite is one of the fastest ways to learn what's actually fragile in a codebase. Curate it like the asset it is.

## Composing with adversarial testing

When [[Tests Should Make Programs Fail|adversarial tests]] reveal a bug, they too produce a regression test — not just the fix, but the originating adversarial case becomes part of the permanent record. Over time, the regression suite becomes a library of attacks the system has been hardened against.

## Related

- [[Tests Should Make Programs Fail]]
- [[Sharp Oracles]]
- [[Three Levels of Tests]]
- [[Reasoning Beside Code]]
