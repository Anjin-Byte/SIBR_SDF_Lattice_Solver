---
title: Three-Layer Stack
tags: [domain, overview]
---

# Three-Layer Stack

The three primary references decompose into three layers that compose cleanly — each lower layer provides what the next layer consumes.

```
      ┌───────────────────────────────────────────────┐
      │  Layer 3 — Physics consumer                   │
      │  Inayat et al. 2016 — pressure-drop model     │
      │  Inputs: open porosity ε_o, window diameter d_w│
      └───────────────────────────────────────────────┘
                            ▲
                            │  ε_o, d_w (both computable from the geometry)
      ┌───────────────────────────────────────────────┐
      │  Layer 2 — Generation pipeline                │
      │  Woodward & Fromen 2021 — seven stages:       │
      │  full vox → boundary vox → population →       │
      │  boundary population → inclusion →            │
      │  intersection → meshing                       │
      └───────────────────────────────────────────────┘
                            ▲
                            │  SDF evaluations, Boolean ops, repetitions
      ┌───────────────────────────────────────────────┐
      │  Layer 1 — Math primitives                    │
      │  Inigo Quilez — exact SDFs for ~30 shapes,    │
      │  Boolean + smooth combinators, operators,     │
      │  positioning, symmetry, repetition,           │
      │  (deformations — which break exactness)       │
      └───────────────────────────────────────────────┘
```

## Why this decomposition matters

### The layers have different change rates

- **Layer 1** is the slowest-changing: these are mathematical primitives with closed-form derivations. Once ported to WGSL, they should rarely change.
- **Layer 2** is moderate: the pipeline structure is stable, but stage implementations, acceleration strategies, and cell-topology support evolve with the project.
- **Layer 3** is application-facing: pressure-drop prediction is one consumer. Future consumers (mass transport, mechanical analysis) would sit at the same level.

This mirrors the architectural decomposition in [[Workspace Architecture MOC]] — Layer 1 is [[sdf crate]] territory (with [[gpu crate]] hosting the matching WGSL), Layer 2 is [[lattice-gen crate]] and [[gpu crate]] working together, Layer 3 is a module inside [[lattice-gen crate]] (or eventually a sibling crate) consuming computed geometric properties.

### Each layer boundary is a natural oracle

- **Layer 1 ↔ Layer 2 boundary:** SDF correctness can be validated by CPU reference implementation (see [[CPU Reference Path]]). Any Woodward-pipeline stage that uses a primitive depends on the primitive being correct first.
- **Layer 2 ↔ Layer 3 boundary:** given a generated lattice, [[Pressure Drop Correlation]] predicts ΔP/L from ε_o and d_w. A **printed-and-measured** ΔP provides a sharp end-to-end oracle for the full stack — exactly the kind of [[Sharp Oracles|sharp oracle]] that [[Tests Should Make Programs Fail]] wants.

### The layers have different test strategies

| Layer | Test style | Oracle |
|---|---|---|
| 1 | Unit + property-based | Closed-form math (distance to a sphere *is* `|p|-r`); reference GLSL compiles and produces known values |
| 2 | Integration, differential | CPU pipeline reproduction; comparison against [[Unit Cell Topologies|known unit-cell geometries]] |
| 3 | End-to-end, external | Measured ΔP on printed parts (Inayat validates <10% MRPE on foams; our lattices are more regular and should do at least as well) |

## What composes and what doesn't

**Composes cleanly:**
- SDFs → lattice body (repetition + intersection at primitive boundary)
- Lattice body → geometric properties (ε_o, d_w, S_v-geo)
- Geometric properties → pressure-drop prediction

**Does NOT compose through the SDF layer cleanly:**
- Non-Lipschitz deformations ([[SDF Primitives Catalog|opTwist, opBend]]) break the distance property and can only be combined up to bounds. Using these inside a sphere-tracing loop requires a step-size discount; using them inside Boolean operations requires care about which operand is exact and which is a bound. Surface this invariant loudly in the code; see [[Reason About Correctness]].

## Related

- [[Domain Knowledge MOC]]
- [[SDF Primitives Catalog]]
- [[Lattice Generation Pipeline]]
- [[Pressure Drop Correlation]]
- [[CPU Reference Path]]
