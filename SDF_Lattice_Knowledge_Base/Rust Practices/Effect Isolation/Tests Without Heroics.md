---
title: Tests Without Heroics
tags: [rust, testing, design, smell]
---

# Tests Without Heroics

> **Heuristic:** if testing requires heroics, the *design* is fighting testability — not the test framework.

## Symptoms of a fighting design

- Temp directories created and torn down for every case.
- Thread sleeps to "wait for" something to happen.
- Network mocks or local server processes spun up per test.
- Giant fixture files duplicated across cases.
- Half the crate marked `pub` purely so tests can reach in.
- `#[cfg(test)]` shims in production modules to "make this testable."

When you see these, the cause is usually upstream: a function that mixes effects with logic, or a type that hides behavior the test needs to observe.

## Healthier alternatives

- **Dependency injection through traits or function parameters.** Pass the clock, filesystem, or GPU engine in. See [[Traits as Seams]].
- **Constructors that accept collaborators.** A type that takes its dependencies in `new()` is testable; one that constructs them internally is not.
- **Explicit config objects.** Beats long parameter lists, and tests can build a known-good config once.
- **Pure helpers factored out of orchestrators.** The orchestrator stays effectful and is integration-tested; the helpers become unit-testable.

## The "moved-out-of-main" pattern

The Rust Book's classic CLI example demonstrates this: code stuck in `main` is hard to test, so move it into a library with explicit inputs. The same move applies recursively — code stuck in any effectful function is hard to test, so move the pure part out.

## Avoid `pub(crate)` as a workaround

Marking things `pub(crate)` "so tests can see them" is a smell. The right moves are:

1. Test through the public API instead. If that's hard, the public API needs adjusting.
2. Move the pure logic to a place it can be tested at the unit level without exposing internals.
3. Use `#[cfg(test)] mod tests { use super::*; }` for tests that legitimately need private access — that's idiomatic and doesn't widen the public surface.

## When heroics are unavoidable

Some integration tests genuinely need a real filesystem, a real GPU adapter, or a real subprocess. That's fine — but those should be **few**, **clearly labeled** (consider a separate `tests/integration_*.rs` file or a feature gate), and not the only way the system is tested. Most tests should still be cheap.

## Related

- [[Pure Core Effectful Edges]]
- [[Traits as Seams]]
- [[Three Levels of Tests]]
- [[Small Functions Narrow Contracts]]
