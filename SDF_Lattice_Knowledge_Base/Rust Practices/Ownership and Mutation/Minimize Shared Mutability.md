---
title: Minimize Shared Mutability
tags: [rust, ownership, mutation, principle]
---

# Minimize Shared Mutability

Rust makes mutation **safe**. That is not the same as **simple**. Stability improves when mutation is **localized and explicit**.

## The preference order

When in doubt, prefer earlier options:

1. **Immutable data flow** — functions take input, return output, no internal state changes.
2. **Returning new values** — instead of mutating an argument, return a transformed version.
3. **Short mutable scopes** — local mutation inside a function, immediately followed by handoff to immutable code.
4. **One owner of state transitions** — a single type or function is responsible for advancing the state; everything else reads.
5. **Multi-owner shared mutability** (`Arc<Mutex<_>>`, `Rc<RefCell<_>>`, atomics) — only when the problem genuinely requires it.

Each step adds complexity. Stop at the level that fits the problem.

## What to be cautious with

**Global mutable state.** `static mut`, `Mutex`-wrapped globals, lazily-initialized singletons. These are bug farms with nice landscaping. Symptoms:

- Tests interfere with each other and have to run single-threaded.
- Behavior depends on prior test order.
- "Reset the world between tests" helpers proliferate.

If global state is genuinely needed (tracing subscriber, the GPU adapter cache), confine it to **one** initialization point and **read-only** access afterward.

**Cross-thread shared mutation.** Each `Arc<Mutex<_>>` is a place where contention, deadlock, or a forgotten lock-order can hide. Prefer message-passing (channels) where it fits.

**Caches that silently change behavior.** A cache that returns stale data, or a cache whose miss path has different observable behavior from its hit path, is one of the hardest categories of bug to reproduce. Either:
- Make the cache pure (same input → same output, always).
- Make the cache an explicit collaborator that tests can swap.

## Why this matters more in Rust

Other languages let you mutate shared state without the compiler complaining; you discover the cost in production. Rust makes the cost visible (locks, atomics, runtime borrow checks) but it does not make the cost go away.

The fact that `Arc<Mutex<T>>` compiles does not mean it's the right design. It means Rust will refuse to make it *unsound* — but it can still make it slow, deadlock-prone, or hard to reason about.

## A short test

When reaching for shared mutability, ask:

- Could this be **one owner** with explicit handoff?
- Could this be **immutable, with a new value produced** on each transition?
- Could the mutation be **scoped to a single function**?
- Is the shared mutability **actually shared at runtime**, or just structurally?

If the honest answer is "yes" to any of these, take that path instead.

## Related

- [[Obvious Ownership]]
- [[Pure Core Effectful Edges]]
- [[Type System as Design Tool]]
