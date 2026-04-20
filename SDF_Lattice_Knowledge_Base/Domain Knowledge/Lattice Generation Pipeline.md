---
title: Lattice Generation Pipeline
tags: [domain, lattice, pipeline, woodward]
source: "Woodward & Fromen 2021, Additive Manufacturing 48, 102386"
source-file: reference/Woodward&Fromen-LatticeGenerationCompensation.pdf
---

# Lattice Generation Pipeline

The seven-stage process for generating a fully **conformal**, **uniform**, **open** lattice from an arbitrary primitive boundary, due to Woodward & Fromen (2021). This is the algorithmic blueprint the [[lattice-gen crate]] implements.

## The pipeline (Woodward Fig. 2)

| Stage | Inputs | Output |
|---|---|---|
| **A. Full voxelization** | Unit cell dimensions, orientation; primitive shape | Voxel grid filling the primitive's bounding region |
| **B. Boundary voxelization** | (same) | The subset of voxels that intersect the primitive's surface |
| **C. Population** | Unit cell geometry (skeleton + connectivity) | Lattice beams placed at every interior voxel |
| **D. Boundary population** | (same) | Lattice beams placed at every boundary voxel, **truncated at the primitive surface** |
| **E. Inclusion** | The full populated lattice; the primitive | Internal-strut subset trimmed to lie inside the primitive |
| **F. Intersection** | The boundary-populated lattice; surfaces describing skin connectivity | The conformal **open lattice skin** — a connected surface formed by joining truncated boundary struts and CAD-edge wireframe elements |
| **G. Meshing** | Strut radius (inflation); Boolean operations | Triangular mesh: `Internal Struts ∪ Wireframe ∪ Boundary Struts ∪ Solid Regions` |

The final output is a triangular mesh suitable for AM (3D printing). The procedure is appropriate for combination strut-surface lattices, solid skins, and solid-lattice hybrids.

## Key properties of the result

The pipeline's output, when correctly implemented, satisfies these invariants:

1. **Conformality.** The lattice's outer boundary follows the primitive surface, regardless of curvature.
2. **Uniformity.** Cell-to-cell, the lattice geometry does not adapt to local features (unlike adaptive meshing). Cell size and orientation are constant inside the primitive.
3. **Openness.** No closed pores. Every void volume is connected to the exterior through some path of windows. This is what makes the lattice usable for fluid applications and resin removal in AM.
4. **Self-supporting (under appropriate cell choice).** No floating beams, no unsupported overhangs at the boundary. Avoids the need for post-processing support removal.
5. **Uniform connectivity.** Every interior strut connects to a node; every truncated boundary strut connects to the skin.

These are the [[Reason About Correctness|invariants]] the implementation must maintain. They are also the properties [[Tests Should Make Programs Fail|adversarial tests]] should attack: degenerate primitives (very thin, very curved, near-zero-volume), unit cells whose periodicity badly aligns with the primitive (forcing many truncated beams), large l_c relative to primitive size, very small l_c (resolution stress).

## Why this pipeline (vs alternatives)

The paper contrasts its approach against three families:

- **Adaptive meshing** (e.g. tetrahedral conformal cell sizing) — produces conformal lattices but breaks uniformity, demands computational resources, and creates per-region variability.
- **Solid skin** — produces fewer paths for resin/material flow, restricts post-processing.
- **Hanging beams** — leaves unsupported truncated members, requiring per-region manual or programmatic supports.

Woodward's approach is *uniform* + *conformal* + *open* + *self-supporting* simultaneously. The price is that it does not adapt to local features — a non-issue for the SDF Lattice Solver's intended applications.

## Software the paper used (and why we are doing it differently)

The paper's implementation chain:
- Rhinoceros 3D 6 + Grasshopper (CAD)
- Crystallon plugin (lattice population)
- Dendro plugin (marching-cubes meshing)
- BullAnt / GeometryGym (Weaire-Phelan generation)
- R for analysis

This is brittle, license-bound, and not headless. **Our project replaces this with a native Rust + wgpu implementation** that:
- Runs headlessly (per [[Headless First]])
- Uses SDFs throughout (per [[SDF Primitives Catalog]]) instead of mesh Booleans
- Is testable per [[Three Levels of Tests]] without the GUI
- Produces deterministic, reproducible output

## Pipeline as state-transition types

Each stage is a natural [[State Transition Types|type-state transition]]. A first-pass shape:

```rust
pub struct PrimitiveSpec { /* ... */ }

pub struct VoxelGrid<Stage> { /* ... */ }
pub struct LatticeSpec<Stage> { /* ... */ }
pub struct Mesh { /* ... */ }

// Marker types for stages
pub struct FullVox;
pub struct BoundaryVox;
pub struct Populated;
pub struct BoundaryPopulated;
pub struct Trimmed;       // after inclusion
pub struct Skinned;       // after intersection

impl PrimitiveSpec {
    pub fn voxelize(&self, cell: CellGeometry) -> VoxelGrid<FullVox>;
}
impl VoxelGrid<FullVox> {
    pub fn extract_boundary(self) -> VoxelGrid<BoundaryVox>;
    pub fn populate(self, cell: CellTopology) -> LatticeSpec<Populated>;
}
// ... and so on through the seven stages
```

This makes every illegal stage transition (e.g. "skin a not-yet-trimmed lattice") a compile error — exactly the [[Push Correctness Left|left-shift]] [[Make Invalid States Unrepresentable]] calls for.

## Parameter conventions (from §3.4 of the paper)

- **L_c** — cell length (mm). Constant across the lattice.
- **r** — strut radius (mm). Constant in a uniform lattice; varies with radial position in functionally-graded lattices ([[Printing Artifacts and Compensation]]).
- **r\*** — dimensionless strut radius, `r* = r / L_c`. Useful for cross-scale comparison.
- **Porosity ρ** — `volume(unit_cell_void) / volume(unit_cell_box)`. Reported per design.

Tested ranges in the paper:
- L_c: 0.5 to 3.5 mm
- Strut diameter: 0.11 to 1.05 mm
- Cells per inch: ~5 (rule of thumb)
- Critical printable strut radius (UMA 90 Black on Carbon M1, ~0.065 mm)

## Maxwell stability classification

All five tested unit cells (Cubic, Kelvin, BCCxy, RD, WP) are **bending-dominated** by Maxwell's stability criterion. This is application-relevant:
- Bending-dominated → lower stiffness, higher cushioning, more deflection under load
- Stretch-dominated → higher stiffness, less cushioning

A future application might call for stretch-dominated cells. The pipeline does not constrain this; the [[Unit Cell Topologies|cell catalog]] should be extensible.

## Project implications

1. **The pipeline is the `lattice-gen` job description.** A `LatticeJob` in [[lattice-gen crate]] describes inputs (primitive, cell, parameters); the engine in [[gpu crate]] executes the seven stages on the GPU; the result is a mesh.
2. **Each stage is testable in isolation.** Per [[Three Levels of Tests]], unit-test each stage's output against a CPU reference; integration-test the full pipeline; doc-test small examples.
3. **The output mesh has measurable properties.** Compute open porosity ε_o and window diameter d_w directly from the geometry — feed those into [[Pressure Drop Correlation]] for downstream prediction.
4. **Adversarial tests** should attack the boundary stages (D, E, F) hardest — they handle the geometric awkwardness of arbitrary primitives.

## Related

- [[Three-Layer Stack]]
- [[Unit Cell Topologies]]
- [[Printing Artifacts and Compensation]]
- [[Pressure Drop Correlation]]
- [[State Transition Types]]
- [[lattice-gen crate]]
- [[sdf crate]]
