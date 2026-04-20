---
title: Domain Knowledge MOC
tags: [moc, domain]
---

# Domain Knowledge — Map of Content

The mathematical, algorithmic, and physical references that anchor what this project actually does. Distilled from primary sources stored in `reference/` at the repo root.

## The thesis (one sentence)

> The SDF Lattice Solver sits in a **three-layer stack**: SDF math primitives at the bottom, a lattice-generation pipeline in the middle, and a downstream physics consumer (pressure-drop prediction) at the top — and each layer has a canonical reference that pins down its required behavior.

See [[Three-Layer Stack]] for the full picture.

## The three primary references

| Reference | File | Role | Vault note |
|---|---|---|---|
| Inigo Quilez, *Distance Functions* (HTML export) | `reference/iq_distfunctions_exact_from_uploaded_html.md` | The kernel library — exact GLSL SDFs for ~30 primitives, combinators, operators, deformations | [[SDF Primitives Catalog]] |
| Woodward & Fromen, *Scalable, process-oriented beam lattices*, Additive Manufacturing 48 (2021) 102386 | `reference/Woodward&Fromen-LatticeGenerationCompensation.pdf` | The seven-stage open-lattice generation pipeline plus a feedforward dimensional-compensation strategy for AM artifacts | [[Lattice Generation Pipeline]], [[Unit Cell Topologies]], [[Printing Artifacts and Compensation]] |
| Inayat et al., *A new pressure drop correlation for open-cell foams*, Chem. Eng. J. 287 (2016) 704–719 | `reference/pressureDropCorrelation.pdf` | A theoretical (no fitted constants) ΔP/L correlation needing only open porosity ε_o and window diameter d_w | [[Pressure Drop Correlation]] |

## Notes in this section

### Reference-derived (from the three primary documents)

- [[Three-Layer Stack]] — how the three references compose
- [[SDF Primitives Catalog]] — the GLSL kernel library to port to WGSL
- [[Lattice Generation Pipeline]] — the seven-stage pipeline (full vox → … → meshing)
- [[Unit Cell Topologies]] — Cubic, Kelvin, BCCxy, RD, WP, sheared/taut variants
- [[Printing Artifacts and Compensation]] — overcure, radial trends, two-iteration feedforward
- [[Pressure Drop Correlation]] — Eq 22 + Eq 35; ε_o and d_w as the only inputs

### Implementation-reasoning (derived during design discussion)

- [[Meshing Complexity Analysis]] — the four candidate meshing methods, time/memory complexity, worked example
- [[Isosurface Extraction Methods]] — MC33 vs DC vs Manifold DC; watertight vs 2-manifold guarantees
- [[SDF vs Mesh Booleans]] — the two operations that share a name and almost nothing else
- [[Mesh Decimation]] — QEM post-process to reduce triangle count to printer resolution
- [[GPU Compute Pipeline]] — the two-pass structure and the three performance pitfalls
- [[AM Platforms]] — per-printer calibration numbers (Formlabs primary, Carbon reference)

## How this section relates to the others

| Section | Answers |
|---|---|
| [[Workspace Architecture MOC]] | What does the system look like in the large? |
| [[Rust Practices MOC]] | What does good code look like at the line level? |
| [[Engineering Philosophy MOC]] | How do we reason and test toward correctness? |
| **Domain Knowledge** (this section) | **What is the system computing, and why?** |

The domain notes are the **specification source**. When a question of the form "what should this function actually compute?" arises, the answer should be traceable to one of these notes — and from there to the primary reference in `reference/`.
