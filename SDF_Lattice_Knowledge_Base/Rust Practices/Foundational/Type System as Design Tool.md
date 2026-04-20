---
title: Type System as Design Tool
tags: [rust, principle, foundational, types]
---

# Type System as Design Tool

The deepest mindset shift in writing strong Rust: **the type system is a design tool, not a compiler obstacle.**

## The reframe

When the borrow checker or type checker pushes back on a design, the instinct from other languages is to **placate** the compiler — `clone`, `Rc`, `RefCell`, `Box<dyn Any>`, `unsafe`, whatever makes the red squiggle go away.

In Rust, that pushback is usually **information**. The compiler is saying "the shape you've drawn doesn't actually fit the constraints you've claimed." The right move is often not to placate but to **redesign**.

## Diagnostic questions

When a Rust design feels like a fight, run through these:

- Should this be **two types** instead of one? (e.g. `Draft` vs `Validated`)
- Should this state be an **enum** instead of two booleans?
- Should this field be **private** with a checked constructor in front?
- Should this dependency be **passed in** rather than looked up globally?
- Should this owned value really be **borrowed**, or this borrow really be **owned**?
- Should this trait object be a **generic**, or vice versa?

Most "Rust is hard" moments resolve when the answer to one of these flips.

## What this looks like in practice

A signal that the type system is doing useful design work:

- The compiler refuses an operation, you change the *types* (not the operation), and the operation becomes unnecessary.
- A bug class disappears from your test plan because it can no longer be expressed.
- Refactors propagate through the codebase by following compile errors instead of running every test.

A signal you're fighting the type system instead of using it:

- Lots of `.clone()` calls that exist only to dodge a borrow.
- Frequent `Rc<RefCell<_>>` in code that isn't actually shared between owners.
- `unsafe` blocks added to express something the type system "doesn't get."
- Long chains of generic parameters whose only job is to stay generic.

## Related

- [[Push Correctness Left]]
- [[Make Invalid States Unrepresentable]]
- [[Obvious Ownership]]
- [[Predictable APIs]]
