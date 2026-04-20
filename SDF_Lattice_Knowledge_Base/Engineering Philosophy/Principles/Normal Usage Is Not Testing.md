---
title: Normal Usage Is Not Testing
tags: [philosophy, principle, testing]
---

# Normal Usage Is Not Testing

> Do not confuse normal usage with real testing.

## The claim

- **Typical user behavior exercises only the obvious paths.** Most users do not stress-test; they accomplish their goal and leave.
- **The hardest and least-traveled code paths are the ones most likely to hide defects.** They have not been worn smooth by use.
- Therefore, **tests must reach the rare, awkward, pathological, composite, and adversarial cases** — because nothing else will.

## Why this matters

A program that "works in normal use" has been validated against approximately the easiest 10–20% of the input space. The remaining 80–90% — empty inputs, max-sized inputs, simultaneous edge conditions, malformed-yet-superficially-valid data, contradictory option combinations, repeated invalid transitions — has been **observed not at all**.

If the only tests in the suite reflect the same shape as normal use, the suite is just **automated normal use**. It catches regressions in the well-traveled paths. It does not catch lurking defects.

The well-traveled paths are not where bugs hide. They are where bugs have been *found and fixed*. The unvisited paths are where bugs **still are**.

## Where the lurking defects live

Categories of code that historically conceal the most bugs:

- **Error paths.** They run rarely, are tested less, and often have their own subtle assumptions.
- **Resource cleanup paths.** Drop, close, abort, cancel. Frequently exercised in tests only on the happy path.
- **Empty cases.** Empty slice, zero-length string, empty iterator, no input at all. Often a separate code path that nobody walked through carefully.
- **Maximum cases.** Max-int, max-capacity, max-recursion-depth. Frequently fail in subtle ways (overflow, allocation failure, stack growth).
- **Concurrent operations.** Two operations interleaved in an order no one anticipated.
- **Repeated operations.** Calling `init` twice, `close` twice, `flush` after `close`.
- **State-machine edge transitions.** Going from a rare state to another rare state.
- **Cross-feature interactions.** Two features that work alone but interact badly.

Each of these is a place ordinary usage will not visit, and where bugs accumulate quietly.

## Implication for the testing strategy

Normal-shaped tests (the "ordinary behavioral" class in [[Output Format]]) are a baseline, not the goal. The goal is to **also** systematically attack the categories above. See:

- [[Tests Should Make Programs Fail]] — the offensive posture
- [[Sharp Oracles]] — the precision of the verdict
- [[Edge Cases and Properties]] — concrete tactical patterns
- [[Tests Without Heroics]] — keep the design testable so the attacks can be cheap

## Honest framing

You will not eliminate all bugs this way. But you will eliminate the **important** bugs that hide in the neglected half of the system — the ones that, left alone, become production incidents. The cost is upfront test design effort. The payoff is the bugs your users do not file.

## Related

- [[Tests Should Make Programs Fail]]
- [[Sharp Oracles]]
- [[Three Levels of Tests]]
- [[Edge Cases and Properties]]
