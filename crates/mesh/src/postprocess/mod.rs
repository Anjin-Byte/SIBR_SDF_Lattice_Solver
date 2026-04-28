//! Mesh post-processing — operations that consume a [`Mesh`](super::Mesh)
//! and return a transformed [`Mesh`](super::Mesh) of the same shape.
//!
//! Where extraction ([`super::marching_cubes`]) and welding ([`super::weld`])
//! are *structural* stages of the mesh pipeline, post-processing is the
//! *quality* stage: visual polish (Taubin smoothing), surface refinement
//! (Butterfly subdivision), and (eventually) topology-preserving reduction
//! (QEM decimation, currently in `xtask`).
//!
//! # Pipeline position
//!
//! ```text
//! sdf -> marching_cubes -> weld -> postprocess::* -> export
//! ```
//!
//! Each post-processing module is independent — they consume and produce
//! the same [`Mesh`](super::Mesh) type, can be chained in any order, and
//! none introduce new transitive dependencies. They are not Cargo-feature
//! gated for that reason; per [Feature Gating](../../../../SDF_Lattice_Knowledge_Base/Architecture/Patterns/Feature%20Gating.md)
//! the rule is "real capability slices," not internal toggles. A caller
//! that does not want smoothing simply does not call [`smooth::taubin`].
//!
//! # Module index
//!
//! - [`smooth`] — Taubin λ/μ low-pass smoothing (Taubin 1995). Reduces
//!   high-frequency mesh noise without volume shrinkage. Topology-preserving.
//! - [`subdivide`] — Butterfly interpolating subdivision (Dyn-Levin-Gregory
//!   1990). Quadruples triangle count per iteration; original vertex
//!   positions are preserved exactly.

pub mod smooth;
pub mod subdivide;
