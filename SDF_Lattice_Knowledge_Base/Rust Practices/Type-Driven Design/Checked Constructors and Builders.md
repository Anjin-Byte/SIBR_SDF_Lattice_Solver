---
title: Checked Constructors and Builders
tags: [rust, idiom, types, type-driven, construction]
---

# Checked Constructors and Builders

Prefer **explicit construction** over post-init mutation. Push validation into birth, not into every later access.

## The shift

```rust
// Worse — every field accessible, every combination possible,
// invariants checked nowhere
let mut config = Config::default();
config.port = 8080;
config.fast = true;
config.other_flag = false;

// Better — invariants enforced by build(), invalid combos rejected once
let config = Config::builder()
    .port(Port::new(8080)?)
    .mode(Mode::Fast)
    .build()?;
```

## Why this is more stable

- Invariants are checked **once**, at construction. Downstream code trusts the type.
- The set of legal configurations is enumerated by the builder API. Combinations the builder can't produce are de facto banned.
- Adding a required field is a compile-time refactor (callers fail to build), not a "good luck finding all the call sites" hunt.
- The constructed type can have all-private fields, eliminating "someone mutated this in the middle of an operation" bugs.

## When a constructor is enough

For 1–3 fields with simple validation, a plain `new` is fine:

```rust
impl Config {
    pub fn new(port: Port, mode: Mode) -> Self { /* ... */ }
}
```

Reach for a builder when:
- there are many optional parameters
- some parameters depend on others
- defaults are non-trivial
- adding a field shouldn't break callers

## Builder shape

Two common Rust builder styles:

**Consuming (movable) builder** — typical:
```rust
pub struct ConfigBuilder { /* ... */ }
impl ConfigBuilder {
    pub fn port(mut self, p: Port) -> Self { self.port = Some(p); self }
    pub fn build(self) -> Result<Config, ConfigError> { /* ... */ }
}
```

**Type-state builder** — for required fields the compiler must enforce:
```rust
pub struct ConfigBuilder<HasPort> { /* phantom marker */ }
impl ConfigBuilder<NoPort> {
    pub fn port(self, p: Port) -> ConfigBuilder<HasPort> { /* ... */ }
}
impl ConfigBuilder<HasPort> {
    pub fn build(self) -> Config { /* ... */ }
}
```

The type-state form costs more code but turns "forgot to set port" into a compile error.

## Constructor naming

- `new(...)` — the canonical constructor for the type
- `with_*` — alternative constructors emphasizing a chosen parameter
- `from_*` — conversion-style construction (often pairs with `From`/`TryFrom`)
- `builder()` — returns a builder
- `default()` — zero/identity value (via `Default`); only when "default" is unambiguous

## Anti-pattern

```rust
pub struct Config { pub port: u16, pub fast: bool, pub mode: Mode, /* ... */ }
```

All-public fields look ergonomic but mean every caller can construct any combination, including invalid ones. Save this shape for **plain data records** with no invariants — config DTOs read from disk that you immediately validate into a private-fields version.

## Related

- [[Make Invalid States Unrepresentable]]
- [[State Transition Types]]
- [[Newtypes and Domain Types]]
- [[Predictable APIs]]
