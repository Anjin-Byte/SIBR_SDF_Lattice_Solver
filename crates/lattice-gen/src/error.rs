//! Typed errors returned by `lattice-gen` constructors.
//!
//! Errors from the underlying [`sdf`] crate are wrapped in [`LatticeError::Sdf`]
//! so callers need only match one error type, while still being able to
//! examine the exact underlying SDF construction failure.

use thiserror::Error;

/// Construction failure for lattice jobs and their sub-types.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum LatticeError {
    /// Underlying SDF construction failed. See the wrapped error for details.
    #[error("sdf construction error: {0}")]
    Sdf(#[from] sdf::BuildError),

    /// Strut radius is too large relative to the unit cell length.
    ///
    /// Each topology imposes its own tightest `r_max(L)` — the radius beyond
    /// which non-adjacent struts in the same cell would overlap (or, for
    /// Cubic, struts from adjacent cells would touch across every shared
    /// face). The caller-visible rule: `radius < max_radius`.
    ///
    /// Derivations per topology (see `SDF_Lattice_Knowledge_Base` and the
    /// `job.rs` `r_max_for` helper):
    /// - **Cubic**: `max_radius = L / 2` (inter-cell tiling constraint).
    /// - **Kelvin**: `max_radius = L / (4·√2) ≈ 0.177·L` (binding case is
    ///   parallel square-face edges on opposite sides of a TO square face).
    /// - **`BCCxy`**: `max_radius = L / (2·√6) ≈ 0.204·L` (binding case is the
    ///   skew body-diagonal vs top-face edge sharing no endpoint).
    #[error(
        "strut too thick for {topology}: radius {radius} must be < {max_radius} \
         (cell length = {cell_length})"
    )]
    StrutTooThick {
        /// The rejected strut radius.
        radius: f32,
        /// The cell length it was compared against.
        cell_length: f32,
        /// The topology whose threshold was violated.
        topology: &'static str,
        /// The tightest acceptable radius (exclusive upper bound) for this
        /// topology at the given cell length.
        max_radius: f32,
    },
}
