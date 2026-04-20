---
title: Patterns Index
tags: [index, patterns]
---

# Patterns Index

Cross-cutting design patterns the workspace follows. Each one is a rule with a rationale.

## Architectural rules

- [[GPU Boundary]] — wgpu types stay inside [[gpu crate]]
- [[Headless First]] — windowing is opt-in
- [[CPU Reference Path]] — every GPU algorithm has a CPU oracle
- [[Binary Strategy]] — one binary vs. two

## Cargo / build

- [[Feature Gating]] — features represent real capability slices
- [[Shared Dependencies]] — what belongs in `[workspace.dependencies]`
- [[Workspace Lints and Profiles]] — root-only configuration

## Platform

- [[winit Lifecycle]] — modern event-loop expectations

## Related

- [[Workspace Architecture MOC]]
- [[Core Principles]]
