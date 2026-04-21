//! SDF-based lattice generation.
//!
//! # Overview
//!
//! This crate composes primitives from the [`sdf`] crate into **lattice body
//! SDFs** — the input specifications `lattice-gen` accepts, and the output
//! SDFs that downstream tooling (GPU meshing, CPU reference meshing) consumes.
//!
//! The entry point is [`lattice_body`]: given a validated [`LatticeJob`]
//! (primitive shape + unit cell + strut radius), it returns an `ExactSdf`
//! representing the tiled-and-trimmed lattice body.
//!
//! # Current topology support
//!
//! - **Cubic**, **Kelvin** (truncated octahedron), and **`BCCxy`** (vertex
//!   octahedron) unit cells. Remaining topologies from the Woodward paper
//!   (Rhombic Dodecahedron, Weaire-Phelan, curved-beam demo cells) arrive
//!   in later phases.
//! - **Cube and Cylinder primitive shapes.** Arbitrary user-supplied SDF
//!   primitives arrive later.
//! - **Uniform strut radius.** Functionally graded radius arrives later.
//!
//! # Property-query coverage
//!
//! [`LatticeJob::open_porosity`], [`LatticeJob::window_diameter`], and
//! [`LatticeJob::specific_surface_area`] currently have closed-form
//! implementations only for the cubic topology. Calling them on a Kelvin
//! or `BCCxy` job panics via `todo!`. Deferred to a Phase-2 follow-up; the
//! CLI only logs these for cubic jobs.
//!
//! # Phase 1c+1d additions
//!
//! - [`LatticeJob::open_porosity`] and related geometric property queries.
//! - CPU reference meshing via [`mesh()`](crate::mesh) — marching cubes on
//!   a dense voxel grid. Slow but correct; intended as the oracle for the
//!   later GPU meshing path.
//!
//! # Example
//!
//! ```
//! use glam::{Vec3, vec3};
//! use lattice_gen::{lattice_body, LatticeJob, PrimitiveShape, StrutSpec, UnitCell};
//! use sdf::Sdf;
//!
//! let job = LatticeJob::new(
//!     PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
//!     UnitCell::cubic(2.0).unwrap(),
//!     StrutSpec::uniform(0.2).unwrap(),
//! ).unwrap();
//!
//! let body = lattice_body(&job);
//! // At a query point outside the primitive, the SDF is positive.
//! assert!(body.eval(vec3(100.0, 0.0, 0.0)) > 0.0);
//! ```

#![forbid(unsafe_code)]

pub mod cell;
pub mod error;
pub mod job;
pub mod mesh;
pub mod primitive;
pub mod progress;
pub mod properties;
pub mod strut;

pub use cell::UnitCell;
pub use error::LatticeError;
pub use job::{LatticeJob, lattice_body};
pub use mesh::{GridSpec, Mesh};
pub use primitive::PrimitiveShape;
pub use progress::Progress;
pub use strut::StrutSpec;
