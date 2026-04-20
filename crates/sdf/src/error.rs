//! Typed errors returned by fallible constructors.
//!
//! All primitive and operator constructors that can fail return
//! `Result<Self, BuildError>`. Each variant records the field name and the
//! offending value so the caller can produce a precise diagnostic without
//! string parsing.

use thiserror::Error;

/// Construction failure for primitives and operators.
#[derive(Debug, Clone, Error, PartialEq)]
pub enum BuildError {
    /// A value that must be strictly positive was zero or negative.
    #[error("{field}: value {value} must be strictly positive")]
    NonPositive {
        /// The constructor parameter that was invalid.
        field: &'static str,
        /// The offending value.
        value: f32,
    },

    /// A value that must be non-negative was negative.
    #[error("{field}: value {value} must be non-negative")]
    Negative {
        /// The constructor parameter that was invalid.
        field: &'static str,
        /// The offending value.
        value: f32,
    },

    /// A value that must be finite was `NaN` or `±∞`.
    #[error("{field}: value {value} is not finite")]
    NonFinite {
        /// The constructor parameter that was invalid.
        field: &'static str,
        /// The offending value.
        value: f32,
    },

    /// The supplied parameters describe a degenerate shape that has no
    /// well-defined SDF (e.g., coincident capsule endpoints, zero-volume box).
    #[error("degenerate geometry: {reason}")]
    Degenerate {
        /// Short explanation of why the geometry is degenerate.
        reason: &'static str,
    },
}
