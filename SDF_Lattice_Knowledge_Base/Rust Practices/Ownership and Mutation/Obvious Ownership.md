---
title: Obvious Ownership
tags: [rust, ownership, principle]
---

# Obvious Ownership

> A lot of brittle Rust comes from **fighting the borrow checker** instead of using it as a design signal.

## Healthy ownership tendencies

- **Pass `&T` when borrowing is natural** — read-only access, function doesn't outlive the caller's value.
- **Take owned `T` when consuming is natural** — function will store, transform, or destructure the value.
- **Clone deliberately, not as a reflex** — every `.clone()` should answer "why is this an owned copy needed here?"
- **Avoid long-lived shared mutable state** — see [[Minimize Shared Mutability]].
- **Build values immutably, then hand them off** — avoid mutating-while-borrowing patterns.

## Reading the borrow checker as a signal

When the borrow checker pushes back, it's worth asking which kind of pushback it is:

| Symptom | Likely design issue |
|---|---|
| "cannot borrow as mutable because also borrowed as immutable" | Two responsibilities tangled in one scope; split them |
| "lifetime parameter required" on many functions | A reference is escaping a scope it shouldn't; consider owning |
| Frequent `clone()` to satisfy the checker | The boundaries between owned and borrowed data are wrong |
| Frequent `Rc<RefCell<_>>` | Shared mutability is being modeled where ownership would be cleaner |
| Long generic lifetime annotations propagating everywhere | A struct is holding a reference when it should hold the value |

The fix is rarely "add another lifetime parameter." It's usually "change which thing owns which other thing."

## Common shape that works

Most well-shaped Rust functions look like one of these:

```rust
// Read: borrow immutably, return owned result
fn render(scene: &Scene) -> Image { /* ... */ }

// Modify: borrow mutably, return nothing or a status
fn normalize(field: &mut Field) { /* ... */ }

// Consume: take owned, return owned (transformation)
fn into_validated(self) -> Validated { /* ... */ }

// Build: take pieces, return owned
fn assemble(parts: Parts) -> Whole { /* ... */ }
```

These four cover most cases. Reach for `Rc`, `Arc`, `RefCell`, or interior mutability only when the actual sharing/mutation pattern genuinely requires it.

## When `Rc<RefCell<_>>` IS the right answer

- Shared graph structures with cycles (rare in this project's domain).
- GUI/event-loop nodes that need observer-style callbacks.
- Test doubles that need to record calls from inside an injected dependency.

These exist. The issue is using them as a *first* reach instead of a deliberate one.

## "Build immutably, hand off" pattern

```rust
// Build (mutable scope, local)
let mut buffer = Vec::with_capacity(n);
for item in input { buffer.push(transform(item)); }

// Hand off (immutable from here on)
process(&buffer);
```

Mutation is contained inside the construction phase. Once `buffer` is handed to `process`, no one mutates it. This is much easier to reason about than a buffer that's modified throughout its lifetime by multiple functions.

## Related

- [[Type System as Design Tool]]
- [[Minimize Shared Mutability]]
- [[Small Functions Narrow Contracts]]
