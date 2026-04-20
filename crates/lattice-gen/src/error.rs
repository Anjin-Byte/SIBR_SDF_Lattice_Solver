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
    /// To guarantee that each strut's non-zero region fits inside a single
    /// periodic cell — the precondition for [`sdf::Repeat`] to yield exact
    /// distances — the radius must be strictly less than half the cell length.
    /// At the limit, struts from adjacent cells would touch across every
    /// shared face, eliminating any open pore structure.
    #[error(
        "strut too thick: radius {radius} must be < half the cell length {cell_length} / 2 = {}",
        cell_length / 2.0
    )]
    StrutTooThick {
        /// The rejected strut radius.
        radius: f32,
        /// The cell length it was compared against.
        cell_length: f32,
    },
}
