---
title: Rust Practices MOC
tags: [moc, rust, practices]
---

# Rust Practices — Map of Content

Principles and idioms the project optimizes for. Entry point for understanding **how we write Rust here** and **why those choices favor stability and testability**.

## The thesis (one sentence)

> Push correctness as far left as possible — into types, construction, ownership, module boundaries, and tests — so fewer bugs are left to runtime chance.

See [[Push Correctness Left]].

## Foundational

- [[Push Correctness Left]] — the meta-principle behind everything else
- [[Type System as Design Tool]] — types are not just a compiler obstacle
- [[Predictable APIs]] — predictability over cleverness

## Type-Driven Design

- [[Make Invalid States Unrepresentable]]
- [[Newtypes and Domain Types]]
- [[State Transition Types]]
- [[Checked Constructors and Builders]]

## Effect Isolation

- [[Pure Core Effectful Edges]]
- [[Traits as Seams]]
- [[Tests Without Heroics]]

## Error Handling

- [[Result vs Panic]]
- [[Domain Errors at Boundaries]]

## Ownership and Mutation

- [[Obvious Ownership]]
- [[Minimize Shared Mutability]]
- [[Unsafe Quarantine]]

## Functions and Data

- [[Small Functions Narrow Contracts]]
- [[Boring Data Layouts]]

## Testing

- [[Three Levels of Tests]] — unit, integration, doc
- [[Edge Cases and Properties]]

## Tooling and Quality

- [[Clippy as Discipline]]
- [[CI Quality Bar]]
- [[Documentation as Truth]]

## The condensed checklist

See [[Rust Practices Checklist]] for the short-form version.

## Relationship to architecture

These practices are how individual crates are written. The workspace shape that hosts them is documented in [[Workspace Architecture MOC]]. The two should be read together — for example, [[Pure Core Effectful Edges]] is the per-crate expression of the [[GPU Boundary]] architectural rule.
