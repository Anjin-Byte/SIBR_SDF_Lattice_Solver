---
title: State Transition Types
tags: [rust, idiom, types, type-driven]
---

# State Transition Types

When a value moves through phases, **model the phases as distinct types**. Use the type system to track where in its lifecycle a value currently is.

## The pattern

```rust
pub struct ConfigDraft { /* mutable, partial */ }
pub struct Config { /* immutable, validated */ }

impl ConfigDraft {
    pub fn new() -> Self { /* ... */ }
    pub fn port(mut self, p: Port) -> Self { /* ... */ }
    pub fn build(self) -> Result<Config, ConfigError> { /* validate */ }
}

impl Config {
    // Methods only valid on a fully-built config
    pub fn run(&self) -> Result<(), RunError> { /* ... */ }
}
```

A function that takes `&Config` knows — without re-checking — that the value has been validated. A function that takes `&mut ConfigDraft` knows the value is still mutable and incomplete.

## Common pairs

| Before | After |
|---|---|
| `ConfigDraft` | `Config` |
| `UncheckedMesh` | `ValidatedMesh` |
| `RawInput` | `ParsedInput` |
| `DisconnectedClient` | `ConnectedClient` |
| `OpenFile` | `ClosedFile` |
| `EncoderRecording` | `EncoderFinished` |

This is just [[Make Invalid States Unrepresentable]] extended over time.

## Why it prevents whole bug categories

- **No "half-built object" bugs.** A consumer cannot accept a partially-initialized `Config` because the type doesn't exist in that shape.
- **No re-validation drift.** Once validated, the validated form carries the proof; nothing downstream needs to revalidate.
- **No invalid operation calls.** You can't call `.send()` on a `DisconnectedClient` because that method only exists on `ConnectedClient`.

## Type-state vs runtime state

Two ways to model "what phase am I in":

1. **Type-state** (this note): the phase is part of the type. Transitions consume the old value and return the new one. Compile-time enforcement.
2. **Runtime state** (an enum field): the phase is data inside the type. Methods check it and return errors. Runtime enforcement.

Prefer type-state when:
- Phase transitions are rare and explicit
- Calling a wrong-phase method should be a compile error, not a runtime one
- The number of phases is small (2–4)

Prefer runtime state when:
- Phases change frequently, dynamically, or based on external events
- The type needs to be stored uniformly (e.g. all clients in one `Vec`, regardless of state)

You can also combine: a runtime enum for the dynamic case, but distinct constructors that produce a value already known to be in a particular state.

## Cost

Type-state can make container types awkward (heterogeneous collections need a wrapping enum). Use it where the safety win is worth the API cost — usually at boundaries of important state machines, not on every flag.

## Related

- [[Make Invalid States Unrepresentable]]
- [[Checked Constructors and Builders]]
- [[Push Correctness Left]]
