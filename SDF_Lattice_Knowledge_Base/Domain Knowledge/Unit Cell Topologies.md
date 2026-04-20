---
title: Unit Cell Topologies
tags: [domain, lattice, cells, woodward]
source: "Woodward & Fromen 2021, §3.2, Fig. 3"
---

# Unit Cell Topologies

The unit cells the solver must support, with their geometric definitions, printability characteristics, and application-relevant properties.

A **unit cell** is the repeat pattern that tiles the interior of the primitive. It has:
- A **skeleton** — the set of nodes and edges (struts) that define connectivity
- A **connectivity description** — which internal struts connect to which boundary elements (needed for skin generation)
- A **bounding box** — usually a cube of side `L_c`, but can be non-cubic (prismatic) in principle

## The standard cell set (Woodward Fig. 3)

### Cubic
- **Topology:** simple cubic lattice — struts along the 12 edges of a cube
- **Connectivity:** every node connects to 6 neighbors (6-valence)
- **Printability:** high (horizontal struts are supported by vertical ones; boundary is well-behaved)
- **Properties:** lowest cell count in the family; highest porosity achievable for given r*; bending-dominated
- **Use as test case:** the paper's primary benchmark (§4.2, §4.3, §4.4). **First cell to implement** — it has the simplest topology and the richest body of reference data.

### Kelvin cell (truncated octahedron)
- **Topology:** 14-faced polyhedron (8 hexagons + 6 squares); struts along the 36 edges
- **Connectivity:** 4-valence nodes; tessellates space by translation alone
- **Printability:** moderate — some struts are near-horizontal
- **Properties:** the classical reference for open-cell foams (see [[Pressure Drop Correlation]]; the Inayat correlation is derived from a tetrakaidecahedron/Kelvin model). Bending-dominated.
- **Relevance:** bridges the synthetic-lattice world of Woodward with the foam literature of Inayat.

### BCCxy (body-centered cubic with XY connections)
- **Also called:** vertex octahedron
- **Topology:** octahedral struts (body diagonals) plus edge struts in the XY plane
- **Connectivity:** high, with struts crossing the cell body
- **Printability:** moderate; XY connections provide lateral support
- **Properties:** denser than cubic at equivalent r*; bending-dominated

### Rhombic Dodecahedron (RD)
- **Topology:** 12-faced polyhedron (all rhombic faces); tessellates by translation alone
- **Connectivity:** 4-valence at vertices, 3-valence at face-centers
- **Printability:** higher porosity failures than BCCxy/Cubic at small L_c in the paper's tests
- **Properties:** bending-dominated; space-filling without gaps

### Weaire-Phelan (WP)
- **Topology:** two cell shapes (irregular dodecahedron + 14-faced Goldberg polyhedron) tessellating together
- **Connectivity:** complex; generated via volume-filling structure
- **Printability:** fewest failures in the paper's tests among the standard five
- **Properties:** lowest surface area for a given volume partition (Kelvin-beating foam geometry); bending-dominated
- **Implementation note:** WP has two distinct cell shapes in the repeat unit — the population stage needs a two-element cell template, not a single one.

## Curved-beam demonstration cells

These show the pipeline's robustness to curvilinear elements. Not a first-wave implementation priority.

### Sheared Cubic
- Cubic topology but with curved struts (cubic edges replaced with smooth curves)
- Demonstrates that the pipeline handles non-straight beams

### Taut Cubic
- Cubic variant with additional support structures in the wireframe description
- Designed to avoid hanging beams along the print axis
- Demonstrates that connectivity-description edits can produce self-supporting variants of otherwise-unsupported topologies

## Maxwell stability

All cells above are **bending-dominated**. This means under load, deflection is dominated by bending of struts rather than axial stretching. Properties:
- Lower stiffness than stretch-dominated equivalents
- Higher cushioning / energy absorption
- More deflection before failure

The paper notes (§4.2) that stretch-dominated structures may be preferable when rigidity is the priority, but all five standard tested cells are bending-dominated.

## Parameter space

For each cell, the design space is parametrized by:
- **L_c** — cell length (uniform; non-uniform prismatic cells exist but are out of first-wave scope)
- **r** — strut radius
- **r\* = r / L_c** — dimensionless strut radius
- Derived: porosity ρ, open porosity ε_o, window diameter d_w, geometric specific surface area S_v-geo

## Suggested type shape in Rust

```rust
pub enum CellTopology {
    Cubic,
    Kelvin,
    BccXy,
    RhombicDodecahedron,
    WeairePhelan,
    // Curved-beam demonstrations (later phase):
    ShearedCubic,
    TautCubic,
}

pub struct UnitCell {
    pub topology: CellTopology,
    pub length: Millimeters,   // L_c
    pub strut_radius: Millimeters, // r
    // derived, computed on demand:
    // dimensionless_radius = strut_radius / length
}
```

This is a simple enum for now; if per-cell data grows large (shape factors, default radius ranges, connectivity templates), factor into a trait or per-variant struct.

## Project implications

1. **Start with Cubic only.** The paper benchmarks against it extensively; it is the simplest test case for the full pipeline.
2. **Kelvin next.** Bridges to the foam literature ([[Pressure Drop Correlation]]) and is the most-cited reference cell.
3. **BCCxy, RD, WP in order.** WP last because its two-cell repeat complicates the population stage.
4. **Curved-beam cells in a later phase.** They test robustness but are not necessary for first-wave applications.
5. **Adversarial test candidates:** for each cell, exercise degenerate dimensions (L_c near strut diameter, strut diameter at printer resolution limit), near-critical r* values (the paper's printability Voronoi diagrams, Fig. 7, are sources of useful failure cases).

## Related

- [[Lattice Generation Pipeline]]
- [[Printing Artifacts and Compensation]]
- [[Pressure Drop Correlation]] — the foam correlations use Kelvin-family (tetrakaidecahedron) geometry
- [[State Transition Types]]
