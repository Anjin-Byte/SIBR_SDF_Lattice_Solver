---
title: Make Invalid States Unrepresentable
tags: [rust, principle, types, type-driven]
---

# Make Invalid States Unrepresentable

The Rustiest answer to a whole class of bugs: don't write a runtime check, write a type that the bad value cannot inhabit.

## The shift

This pattern turns many tests from *"check this can't happen"* into *"the compiler won't let it happen."* That is the cleanest possible left-shift along [[Push Correctness Left]].

## Techniques

### 1. Prefer `enum` over booleans with hidden meaning

```rust
// Bad
fn render(stereo: bool, hdr: bool) { /* ... */ }
// Caller: render(true, false)  — what does that mean?

// Better
enum AudioMode { Mono, Stereo }
enum ColorRange { Sdr, Hdr }
fn render(audio: AudioMode, color: ColorRange) { /* ... */ }
```

Boolean parameters are a 50/50 chance of being remembered correctly. Enums document themselves at the call site.

### 2. Prefer newtypes over raw primitives when values have semantics

See [[Newtypes and Domain Types]] for the full pattern.

### 3. Keep fields private, force construction through checked constructors

See [[Checked Constructors and Builders]].

### 4. Separate validated from unvalidated data

See [[State Transition Types]].

## A worked example

```rust
pub struct Port(u16);

impl Port {
    pub fn new(value: u16) -> Option<Self> {
        if value == 0 { None } else { Some(Self(value)) }
    }

    pub fn get(self) -> u16 { self.0 }
}
```

Once a function takes `Port`, no test inside that function ever needs to assert `port != 0` again — it cannot be 0. The check moved left, exactly once, to the constructor.

Compare the alternative: every function downstream that uses a `u16` "port" either re-checks (defensive, wasteful) or trusts (fragile). Both are worse than just having a `Port`.

## When this is overkill

For a one-shot internal helper that takes a `u16` and immediately uses it once, a newtype is ceremony. The principle scales with **how many functions touch the value** and **how bad the consequence of a bad value is**.

Rule of thumb: if more than ~3 functions transmit the value across module boundaries, give it a type.

## Related

- [[Push Correctness Left]]
- [[Newtypes and Domain Types]]
- [[State Transition Types]]
- [[Checked Constructors and Builders]]
