---
title: Edge Cases and Properties
tags: [rust, testing, principle]
---

# Edge Cases and Properties

Test the places where humans improvise badly. Happy-path tests are necessary but rarely sufficient.

## What to actually target

The high-value test categories for stability:

- **Boundary values** — zero, one, max, max+1, empty, full.
- **Empty inputs** — empty slice, empty string, empty iterator. Many bugs live here.
- **Malformed inputs** — wrong types, wrong order, garbage bytes, partial reads.
- **Repeated application / idempotence** — does `f(f(x)) == f(x)` when it should?
- **Ordering assumptions** — does the function depend on input order? Does it preserve order?
- **State transition mistakes** — calling methods in the wrong sequence; double-init, double-close.
- **Concurrency boundaries** — two operations racing; cancellation mid-flight (where applicable).

## A useful question

For every public function, ask:

> **What inputs would an annoyed future version of me accidentally pass here?**

Then write *that* test now, while current-you is still cooperative. The annoyed-future-you scenarios are exactly where bugs land in production.

## Specific Rust patterns

### Boundary tests for integer types

```rust
#[test]
fn handles_zero() { assert!(Port::new(0).is_none()); }
#[test]
fn handles_one() { assert!(Port::new(1).is_some()); }
#[test]
fn handles_max() { assert!(Port::new(u16::MAX).is_some()); }
```

### Empty-collection tests

```rust
#[test]
fn empty_input_returns_empty_output() {
    assert_eq!(transform(&[]), vec![]);
}
```

### Malformed-input tests

```rust
#[test]
fn rejects_truncated_header() {
    let bad = &[0x00, 0x01]; // valid headers are 8 bytes
    assert!(matches!(parse(bad), Err(ParseError::TruncatedHeader)));
}
```

### Idempotence

```rust
#[test]
fn validate_is_idempotent() {
    let v1 = validate(input.clone()).unwrap();
    let v2 = validate(input).unwrap();
    assert_eq!(v1, v2);
}
```

## Property-based testing

For algorithms with rich input spaces, consider `proptest` or `quickcheck`:

```rust
proptest! {
    #[test]
    fn round_trip_serialization(input in any::<Config>()) {
        let bytes = serialize(&input).unwrap();
        let back = deserialize(&bytes).unwrap();
        prop_assert_eq!(input, back);
    }
}
```

Properties to look for in this project:
- GPU result ≈ CPU reference within tolerance (see [[CPU Reference Path]])
- Serialization round-trips
- Idempotent transformations
- Invariants preserved by every operation on a type

## What NOT to over-test

- Trivial getters/setters — testing them tests the language, not your code.
- Implementation details — tests that break on every refactor are noise.
- Behavior already enforced by the type system — if `Port` can't be 0, no test needs to check `Port != 0`.

The goal is tests that catch real defects without rotting on every change. Test the **behavior**, not the **implementation**.

## Related

- [[Three Levels of Tests]]
- [[Tests Without Heroics]]
- [[Make Invalid States Unrepresentable]]
- [[CPU Reference Path]]
