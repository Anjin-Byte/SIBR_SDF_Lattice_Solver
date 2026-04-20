---
title: Boring Data Layouts
tags: [rust, data, principle, design]
---

# Boring Data Layouts

For long-term stability, prefer **boring** data shapes over magical ones. A little verbosity in modeling buys a lot of reliability.

## The preferences

- **Obvious structs beat clever tuple soup.**
  ```rust
  // Worse
  fn move_to(pos: (f32, f32, f32)) { /* which one is x? */ }
  // Better
  pub struct Point3 { pub x: f32, pub y: f32, pub z: f32 }
  fn move_to(pos: Point3) { /* unambiguous */ }
  ```

- **Explicit config types beat long parameter lists.**
  ```rust
  // Worse
  fn run(input: &Path, output: &Path, threads: usize, verbose: bool, dry_run: bool, /* ... */) { /* ... */ }
  // Better
  pub struct RunOptions { /* named fields */ }
  fn run(opts: &RunOptions) { /* ... */ }
  ```

- **Named states beat magic combinations of flags.**
  ```rust
  // Worse — what does (true, false, true) mean?
  fn render(stereo: bool, hdr: bool, denoise: bool) { /* ... */ }
  // Better
  pub enum AudioMode { Mono, Stereo }
  pub enum ColorPipeline { Sdr, Hdr }
  pub enum Denoising { Off, On }
  ```

## Where this matters most

These habits pay off disproportionately in:

- **Systems code** — long-lived state structures whose layout future contributors must read.
- **GPU code** — buffer layouts that need to match WGSL exactly; a misnamed field corrupts a pipeline silently.
- **Parsing** — token / AST types where structural clarity prevents whole bug classes.
- **Concurrency work** — shared structures where readers must understand the layout at a glance.

In each, "compact" often becomes "cryptic." The compactness saves you keystrokes once and costs you understanding for years.

## When tuple structs are fine

- One-field newtypes — see [[Newtypes and Domain Types]].
- Genuinely positional pairs where the position carries the meaning (e.g. `(usize, usize)` for "row, column" if it's local enough that no one will forget).
- Returning a small handful of values from a private helper.

For anything that crosses a module boundary or persists in a struct field, prefer a named struct.

## The "long parameter list" trap

A function with 6+ parameters is almost always bundling several concepts. Common bundlers:

- A `RunOptions` struct grouping "how to run."
- A `Context` struct grouping shared collaborators (clock, logger, allocator).
- A separate function per mode, instead of a `mode` flag.

The fix is rarely "use named arguments" (Rust doesn't have them). The fix is usually "introduce a struct or split the function."

## Default values, explicitly

When a config has many fields, providing defaults is helpful — but not via `Default::default()` if the defaults are non-obvious. Prefer:

```rust
impl RunOptions {
    pub fn balanced() -> Self { /* explicit, named preset */ }
    pub fn aggressive() -> Self { /* explicit, named preset */ }
}
```

Named presets document what "default" actually means at the call site.

## Related

- [[Make Invalid States Unrepresentable]]
- [[Newtypes and Domain Types]]
- [[Predictable APIs]]
- [[Small Functions Narrow Contracts]]
