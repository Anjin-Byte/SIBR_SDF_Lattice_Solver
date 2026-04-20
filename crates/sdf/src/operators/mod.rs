//! Operators — wrappers that transform the query point before delegating to
//! an inner SDF.
//!
//! Every operator in this module is **precision-preserving**: an operator over
//! an [`crate::ExactSdf`] is itself `ExactSdf`, and over a [`crate::BoundSdf`]
//! is `BoundSdf`. The transforms are chosen specifically to preserve the
//! 1-Lipschitz property of the inner SDF:
//!
//! - [`Translate`] — isometry (rigid translation).
//! - [`UniformScale`] — Lipschitz scaling: `inner(p/s) * s` preserves distance
//!   correctness because the scaling compensates for the contracted query frame.
//! - [`SymX`], [`SymY`], [`SymZ`] — reflections across axes (isometries).
//! - [`Repeat`] — folds the query point into a fundamental period domain;
//!   representing an infinite tiling of the inner SDF.
//! - [`LimitedRepeat`] — same, but bounded to a finite extent of cells.
//!
//! # Lattice relevance
//!
//! [`Repeat`] and [`LimitedRepeat`] are the operators that make lattice
//! generation computationally feasible. Rather than evaluating `N_struts`
//! primitive SDFs per query, the query point is folded into a single unit cell,
//! reducing per-query cost to `O(k)` (struts per cell) regardless of part size.
//! See the Domain Knowledge note "GPU Compute Pipeline" for the full
//! complexity analysis.

pub mod limited_repeat;
pub mod repeat;
pub mod sym;
pub mod translate;
pub mod uniform_scale;

pub use limited_repeat::LimitedRepeat;
pub use repeat::Repeat;
pub use sym::{SymX, SymY, SymZ};
pub use translate::Translate;
pub use uniform_scale::UniformScale;
