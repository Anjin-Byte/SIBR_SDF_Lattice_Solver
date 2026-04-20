---
title: Result vs Panic
tags: [rust, error-handling, principle]
---

# Result vs Panic

> **Rule:** `Result` for **recoverable failure**. `panic!` for **bugs and impossible invariants**.

## The distinction

- **`Result<T, E>`** — the failure is something a *caller* might reasonably handle. Examples: file not found, network unreachable, malformed input, GPU adapter not available.
- **`panic!` / `expect(...)`** — the failure means *the program itself is wrong*. Examples: an internal index map missing a key that was just inserted, an assertion about a private invariant, an unreachable branch in a closed enum.

The test for which to use: **"could a sensible caller do anything other than crash?"** If yes → `Result`. If no → panic.

## Where panics belong

- `unreachable!()` for branches the type system can't prove dead but logic guarantees.
- `assert!`/`debug_assert!` for invariants you want spelled out.
- `expect("message")` over `unwrap()` *only* when the message documents the invariant being relied on (e.g. `.expect("port was checked nonzero in constructor")`).
- Inside private code where invariants are local and the code that maintains them is nearby.

## Where panics do NOT belong

- Library-public functions that handle external input (parsing, deserialization, network).
- Functions whose failure mode the user might want to handle differently (e.g. retry, prompt, fall back).
- Anywhere the failure depends on **runtime data** rather than **programmer mistake**.

For library-shaped code, panic paths should be **rare and loud** — when they fire, they're surfacing a real bug, not a runtime situation.

## On `unwrap()`

`unwrap()` on a `Result` or `Option` says "I'm certain this can't fail." Often that certainty is wrong, and the panic message gives no hint what was assumed.

Prefer:
- `expect("explanation")` — at least the message documents the assumption.
- `?` — when propagation is the right move.
- A real `match` — when both arms have meaningful behavior.

Allow `unwrap()` freely in:
- Tests
- `xtask` and one-off scripts
- Examples in docs (where brevity matters more than robustness)

In production code, prefer `expect` over `unwrap` purely for the error-message hygiene.

## A worked example

```rust
// File missing — caller could retry, prompt user, fall back. Result.
fn load_config(path: &Path) -> Result<Config, LoadError> { /* ... */ }

// Index map corrupted despite private construction — programmer error. Panic.
fn lookup(&self, id: VertexId) -> &Vertex {
    self.vertices.get(&id)
        .expect("VertexId can only be obtained from this map")
}
```

The first failure has a sensible caller-side response. The second one shouldn't be possible without a bug, and silently returning `None` would just push the bug downstream.

## Why this matters for stability

Mixing the two roles erodes both. If `Result` is used for "this should never happen," callers get verbose error handling for impossible cases and stop reading it. If `panic` is used for ordinary failures, the program crashes when it should recover.

Keeping the categories sharp keeps both behaviors meaningful.

## Related

- [[Domain Errors at Boundaries]]
- [[Push Correctness Left]]
- [[Predictable APIs]]
