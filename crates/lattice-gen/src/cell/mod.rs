//! Unit cell topologies.
//!
//! A **unit cell** is the repeating pattern tiled across the primitive's
//! volume. Phase 1a supports only [`UnitCell::Cubic`]. Additional topologies
//! from the [Domain Knowledge note](../../SDF_Lattice_Knowledge_Base/Domain%20Knowledge/Unit%20Cell%20Topologies.md)
//! — Kelvin, `BccXy`, Rhombic Dodecahedron, Weaire-Phelan — arrive in later phases.

pub mod cubic;

use crate::error::LatticeError;

/// A validated unit cell specification.
///
/// Constructed via topology-specific helpers (e.g. [`UnitCell::cubic`]) that
/// validate parameters before admitting the value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum UnitCell {
    /// Simple cubic topology: three edge-aligned struts per cell.
    ///
    /// The `length` field is the cell's edge length (mm). Every Cubic cell's
    /// bounding box is `[0, length]^3` (origin-aligned) for the purposes of
    /// internal tiling.
    Cubic {
        /// Edge length of the cubic cell in mm. Strictly positive and finite.
        length: f32,
    },
}

impl UnitCell {
    /// Constructs a Cubic unit cell with the given edge length.
    ///
    /// # Errors
    ///
    /// - [`LatticeError::Sdf`] wrapping [`sdf::BuildError::NonFinite`] if
    ///   `length` is `NaN` or `±∞`.
    /// - [`LatticeError::Sdf`] wrapping [`sdf::BuildError::NonPositive`] if
    ///   `length <= 0`.
    pub fn cubic(length: f32) -> Result<Self, LatticeError> {
        if !length.is_finite() {
            return Err(sdf::BuildError::NonFinite {
                field: "cell.length",
                value: length,
            }
            .into());
        }
        if length <= 0.0 {
            return Err(sdf::BuildError::NonPositive {
                field: "cell.length",
                value: length,
            }
            .into());
        }
        Ok(Self::Cubic { length })
    }

    /// Returns the characteristic cell length — the edge length for Cubic,
    /// equivalent dimension for other topologies when added.
    pub fn length(self) -> f32 {
        match self {
            Self::Cubic { length } => length,
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    missing_docs
)]
mod tests {
    use super::*;

    #[test]
    fn cubic_accepts_positive_length() {
        let c = UnitCell::cubic(2.4).unwrap();
        assert_eq!(c.length(), 2.4);
    }

    #[test]
    fn cubic_rejects_zero() {
        assert!(UnitCell::cubic(0.0).is_err());
    }

    #[test]
    fn cubic_rejects_negative() {
        assert!(UnitCell::cubic(-1.0).is_err());
    }

    #[test]
    fn cubic_rejects_nan_and_inf() {
        assert!(UnitCell::cubic(f32::NAN).is_err());
        assert!(UnitCell::cubic(f32::INFINITY).is_err());
    }
}
