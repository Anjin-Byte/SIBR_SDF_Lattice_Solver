---
title: SDF Lattice Knowledge Base
tags: [index, vault-root]
---

# SDF Lattice Knowledge Base

Working knowledge base for the **SIBR SDF Lattice Solver** — a Rust CLI application built on [[wgpu]] for GPU-accelerated SDF (signed distance field) lattice computation.

## Project shape (one line)

A native, headless-by-default CLI tool with an optional preview window, structured as a Cargo workspace where the GPU backend is firewalled behind a narrow API and a CPU reference path keeps everything testable without a GPU.

## Start here

- [[Workspace Architecture MOC]] — map of content for the workspace design
- [[Rust Practices MOC]] — map of content for code-level principles
- [[Engineering Philosophy MOC]] — map of content for how we reason and test
- [[Domain Knowledge MOC]] — map of content for the math and physics the system computes
- [[Core Principles]] — the four architectural rules
- [[Push Correctness Left]] — the meta-principle behind the Rust practices
- [[Output Format]] — the four-section template for non-trivial work
- [[Three-Layer Stack]] — the SDF → lattice pipeline → pressure-drop decomposition

## Folder map

- **Architecture/** — workspace structure, crate boundaries, and architectural patterns
  - [[Workspace Architecture MOC]] — the entry-point map
  - [[Crates Index]] — per-crate responsibility notes
  - [[Patterns Index]] — cross-cutting design patterns
  - External dependency references: [[wgpu]], [[winit]]
- **Rust Practices/** — principles for writing the Rust code inside crates
  - [[Rust Practices MOC]] — the entry-point map
  - [[Rust Practices Checklist]] — condensed code-review form
  - Subfolders: Foundational, Type-Driven Design, Effect Isolation, Error Handling, Ownership and Mutation, Functions and Data, Testing, Tooling and Quality
- **Engineering Philosophy/** — how we reason about correctness and attack the program with adversarial tests
  - [[Engineering Philosophy MOC]] — the entry-point map
  - [[Output Format]] — the four-section template (correctness model → invariants → design → test strategy)
  - [[Principles Index]] — the eight Knuth-style principles
- **Domain Knowledge/** — the math, algorithms, and physics the system computes, distilled from the primary references in `reference/`
  - [[Domain Knowledge MOC]] — the entry-point map
  - [[Three-Layer Stack]] — how the three reference documents compose
  - Per-reference notes: [[SDF Primitives Catalog]], [[Lattice Generation Pipeline]], [[Unit Cell Topologies]], [[Printing Artifacts and Compensation]], [[Pressure Drop Correlation]]
  - Implementation-reasoning notes: [[Meshing Complexity Analysis]], [[Isosurface Extraction Methods]], [[SDF vs Mesh Booleans]], [[Mesh Decimation]], [[GPU Compute Pipeline]], [[AM Platforms]]

## How the four sections relate

| Section | Answers |
|---|---|
| **Architecture** | What does the system look like in the large? |
| **Rust Practices** | What does good code look like at the line level? |
| **Engineering Philosophy** | How do we *reason and test* toward correctness? |
| **Domain Knowledge** | What is the system computing, and why? |

Domain Knowledge is the **specification source** — when a question of the form "what should this function actually compute?" arises, the answer should be traceable to one of those notes and from there to a reference document in `reference/`.

A worked example of how the four interlock: the [[GPU Boundary]] architectural rule is realized at the line level by [[Pure Core Effectful Edges]]; both are *defended* by the adversarial testing posture of [[Tests Should Make Programs Fail]] and the differential oracles of [[Sharp Oracles]] (e.g. comparing GPU output to a [[CPU Reference Path]]); and what the GPU is *actually computing* — exact vs. bound SDFs, tortuosity-based pressure drop, seven-stage lattice pipeline — is defined by the [[Domain Knowledge MOC|Domain Knowledge notes]] against the primary references.
