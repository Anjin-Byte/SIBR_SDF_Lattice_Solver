---
title: Documentation as Truth
tags: [rust, docs, principle, quality]
---

# Documentation as Truth

Document **invariants** and **examples** where they live. Because Rust doc examples are compiled and tested, documentation can act as **executable truth** — not decoration.

## What good docs explain

For a public type or function, useful docs cover:

- **What the type guarantees.** "A `Port` is always nonzero." The reader should not need to read the implementation to learn this.
- **What the caller must guarantee.** Preconditions, valid input ranges, expected ordering of method calls.
- **Failure modes.** Which errors are returned when, and what they mean operationally.
- **Examples of intended use.** Short, executable, written at the level of a caller.

This is not the same as describing the implementation. Docs should describe the **contract**, not the body.

## Doc tests are the secret weapon

```rust
/// Loads a [`Config`] from a TOML file.
///
/// # Errors
///
/// Returns [`LoadError::Io`] if the file cannot be read,
/// or [`LoadError::Parse`] if the contents are not valid TOML.
///
/// # Examples
///
/// ```no_run
/// use mycrate::Config;
/// let config = Config::load("app.toml")?;
/// # Ok::<(), mycrate::LoadError>(())
/// ```
pub fn load(path: impl AsRef<Path>) -> Result<Config, LoadError> { /* ... */ }
```

Because the example compiles as part of `cargo test --doc`, it can't drift. Rename `load` and the doc test breaks. Change the error type and the doc test breaks. The docs and the code are *married* instead of "separated but technically still in contact."

## Conventions that make docs more useful

- **Front-load the summary.** First line is a single-sentence description, used in API listings.
- **Use sections.** `# Errors`, `# Panics`, `# Safety`, `# Examples` are conventional and parsed by tools.
- **Link types with `[` `]`.** Rustdoc resolves `[`Type`]` to a clickable link.
- **Document `# Panics` for any function that can panic.** Including for `expect` calls in non-obvious places.
- **Document `# Safety` for every `unsafe fn`.** Non-negotiable. See [[Unsafe Quarantine]].

## Where docs are most valuable

- **Public APIs of the [[gpu crate]]** — the narrow surface external code consumes. Worth spending time on.
- **[[Newtypes and Domain Types]]** — they encode invariants; the docs should make those invariants legible.
- **Error enums** — what each variant means, when it's returned.
- **Builders and constructors** — required vs optional, validation rules.

Internal helper functions need less doc weight; their callers are nearby and can read the code.

## When to skip docs

- Trivial getters whose name fully describes them (`get`, `len`, `is_empty`).
- Internal `pub(crate)` items used in one place.
- Generated/derived code where the trait's own docs apply.

This is not "skip docs everywhere convenient." It's "don't write empty boilerplate that says what the signature already says."

## CI for docs

Useful additions in CI:

```sh
# Catch broken intra-doc links
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps

# Doc tests already run via `cargo test --doc` — see [[CI Quality Bar]]
```

These keep documentation honest the same way Clippy keeps code honest.

## Related

- [[Predictable APIs]]
- [[Three Levels of Tests]]
- [[Clippy as Discipline]]
- [[CI Quality Bar]]
