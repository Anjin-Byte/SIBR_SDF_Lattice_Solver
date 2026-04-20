---
title: Unsafe Quarantine
tags: [rust, unsafe, principle]
---

# Unsafe Quarantine

> Unsafe should be treated like plutonium: useful in some reactors, but not something to smear over the kitchen counters.

## The principle

If `unsafe` is necessary, **isolate it behind a safe API with a clearly documented invariant**. Never let raw `unsafe` operations leak into general code paths.

## Anatomy of a healthy unsafe module

A healthy `unsafe` module has:

1. **One narrow purpose.** "FFI to the foo library" or "in-place reinterpret of byte buffers." Not a grab-bag.
2. **A safe wrapper** that callers actually use. The `unsafe` keyword should not appear in the public API of higher-level code.
3. **`# Safety` doc comments** on every public `unsafe fn`, listing the invariants the caller must uphold.
4. **`// SAFETY: ...` comments** above every `unsafe { ... }` block, explaining *why* the operation is sound at this call site.
5. **Tests around the safe wrapper**, including miri runs where applicable.
6. **As little surface area as possible** — minimize the lines that are technically `unsafe`, even within the module.

## The shape

```rust
mod raw {
    /// # Safety
    /// `ptr` must point to a valid `Foo` allocated by `foo_create`,
    /// and must not have been freed.
    pub unsafe fn destroy(ptr: *mut Foo) { /* ... */ }
}

pub struct OwnedFoo {
    ptr: *mut Foo,
}

impl Drop for OwnedFoo {
    fn drop(&mut self) {
        // SAFETY: ptr was obtained from foo_create in OwnedFoo::new
        // and is not aliased (we own it).
        unsafe { raw::destroy(self.ptr) }
    }
}
```

Outside this module, no one writes `unsafe`. The dangerous operations are wrapped in a type that maintains the invariant via ownership.

## When to reach for unsafe

Legitimate reasons:
- FFI to a C library (e.g. some parts of GPU bindings, OS APIs).
- Performance-critical low-level operations (bytemuck-style reinterpretation, SIMD intrinsics) where the safe alternative measurably matters.
- Implementing a fundamental data structure that the safe API can't express (rare in application code).

Bad reasons:
- "The borrow checker is annoying." See [[Type System as Design Tool]] — this is information, not adversity.
- "I want to skip a check." If the check matters, you'll regret skipping it.
- "Performance" without measurements showing the safe version is the bottleneck.

## What this project will likely use it for

For an `wgpu` CLI app, `unsafe` will most likely show up in:

- `bytemuck` casts of plain-old-data structs to byte slices for GPU upload (already abstracted by the crate).
- Raw window handles bridging winit ↔ wgpu (already abstracted by the libraries).
- Possibly small custom math intrinsics if profiling demands it.

In all of these, the goal is to keep the `unsafe` usage **inside the gpu crate**, behind a safe API the rest of the workspace consumes.

## Related

- [[GPU Boundary]]
- [[Obvious Ownership]]
- [[Type System as Design Tool]]
