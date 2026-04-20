---
title: Newtypes and Domain Types
tags: [rust, idiom, types, type-driven]
---

# Newtypes and Domain Types

When a value has **semantic meaning** beyond its primitive representation, give it a type.

## The pattern

```rust
pub struct UserId(u64);
pub struct LatticeIndex(u32);
pub struct Meters(f32);
pub struct VertexCount(u32);
```

These are zero-cost at runtime — same memory layout as the primitive — but distinct at the type level. You cannot accidentally pass a `LatticeIndex` where a `VertexCount` was expected.

## Why this matters more than it looks

- **Function signatures self-document.** `fn get(id: UserId)` is unambiguous; `fn get(id: u64)` could be anything.
- **Refactoring is safe.** Changing `LatticeIndex` from `u32` to `u64` is one definition site, not a grep across the codebase.
- **Mistakes that compile in primitives don't compile in newtypes.** You cannot subtract `Meters` from `VertexCount` by accident.
- **Conversions become explicit.** Need to multiply `Meters` by a scaling factor? You write the operation deliberately.

## Implementation guidance

### Tuple struct vs named field

```rust
pub struct Port(u16);             // tuple: fine for one field
pub struct Port { value: u16 }    // named: fine when extra fields might join
```

Tuple form is conventional for single-field newtypes.

### Constructor shape

- **Infallible:** plain `pub fn new(v: T) -> Self`.
- **Fallible:** `pub fn new(v: T) -> Result<Self, E>` or `Option<Self>` if the failure carries no info.
- **From / TryFrom:** implement when conversion is the natural verb.

```rust
impl TryFrom<u16> for Port {
    type Error = PortError;
    fn try_from(v: u16) -> Result<Self, Self::Error> {
        if v == 0 { Err(PortError::Zero) } else { Ok(Self(v)) }
    }
}
```

### Accessors

- `pub fn get(self) -> T` for `Copy` inner types.
- `pub fn as_inner(&self) -> &T` for non-`Copy`.
- Implement `Deref` only when the newtype is genuinely a "smart pointer" — overusing `Deref` undermines the whole point.

### Derives

Match the inner type's traits where appropriate: `Debug`, `Clone`, `Copy`, `Eq`, `PartialEq`, `Hash`, `Ord`, `PartialOrd`. Add `serde::Serialize` / `Deserialize` when the type crosses a serialization boundary.

## When NOT to bother

- One-shot internal scalar that never crosses a function boundary.
- Trivial enum-like value where `enum` is clearer.
- Cases where the inner primitive's full operation set really is the desired API (`Meters + Meters` makes sense; `Meters * Meters` doesn't — that's `SquareMeters`).

## Related

- [[Make Invalid States Unrepresentable]]
- [[Checked Constructors and Builders]]
- [[Predictable APIs]]
