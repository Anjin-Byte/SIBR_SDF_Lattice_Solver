//! Unit cell topologies.
//!
//! A **unit cell** is the repeating pattern tiled across the primitive's
//! volume. The current supported topologies are:
//!
//! - [`UnitCell::Cubic`] — 3 axis-centered struts (simplest test case).
//! - [`UnitCell::Kelvin`] — 24 capsules modeling the truncated octahedron
//!   (tetrakaidecahedron) edge graph; bridges to the Inayat foam
//!   pressure-drop correlation.
//! - [`UnitCell::BccXy`] — 8 body-diagonal capsules + 4 top-face edges.
//!
//! Remaining topologies from the [Domain Knowledge note](../../SDF_Lattice_Knowledge_Base/Domain%20Knowledge/Unit%20Cell%20Topologies.md)
//! — Rhombic Dodecahedron, Weaire-Phelan, and curved-beam demo cells —
//! arrive in later phases.

pub mod bccxy;
pub mod cubic;
pub mod kelvin;

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
    /// bounding box is `[-length/2, length/2]^3` (origin-centered) for the
    /// purposes of internal tiling.
    Cubic {
        /// Edge length of the cubic cell in mm. Strictly positive and finite.
        length: f32,
    },
    /// Kelvin cell (truncated octahedron).
    ///
    /// Inscribed in the origin-centered cube `[-length/2, length/2]^3`; 24
    /// vertices at permutations of `(0, ±length/4, ±length/2)` connected by
    /// 36 edges (24 home-assigned to this cell under positive-half
    /// convention, 12 supplied by periodic neighbors). See [`kelvin`].
    Kelvin {
        /// Edge length of the bounding cube in mm. Strictly positive and finite.
        length: f32,
    },
    /// `BCCxy` (vertex octahedron) topology.
    ///
    /// 8 body-diagonal struts from each corner of `[-length/2, length/2]^3`
    /// to the cube center, plus 4 edges on the `z = +length/2` face. Bottom
    /// face comes from the `+z` neighbor under periodic tiling. See [`bccxy`].
    BccXy {
        /// Edge length of the bounding cube in mm. Strictly positive and finite.
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
        validate_length(length, "cell.length")?;
        Ok(Self::Cubic { length })
    }

    /// Constructs a Kelvin (truncated octahedron) unit cell.
    ///
    /// # Errors
    ///
    /// Same as [`UnitCell::cubic`]: non-finite or non-positive `length`.
    pub fn kelvin(length: f32) -> Result<Self, LatticeError> {
        validate_length(length, "cell.length")?;
        Ok(Self::Kelvin { length })
    }

    /// Constructs a `BCCxy` (vertex octahedron) unit cell.
    ///
    /// # Errors
    ///
    /// Same as [`UnitCell::cubic`]: non-finite or non-positive `length`.
    pub fn bccxy(length: f32) -> Result<Self, LatticeError> {
        validate_length(length, "cell.length")?;
        Ok(Self::BccXy { length })
    }

    /// Returns the characteristic cell length — the edge length of the
    /// bounding cube for all current topologies.
    pub fn length(self) -> f32 {
        match self {
            Self::Cubic { length } | Self::Kelvin { length } | Self::BccXy { length } => length,
        }
    }
}

/// Shared length-field validation for every [`UnitCell`] constructor.
fn validate_length(length: f32, field: &'static str) -> Result<(), LatticeError> {
    if !length.is_finite() {
        return Err(sdf::BuildError::NonFinite {
            field,
            value: length,
        }
        .into());
    }
    if length <= 0.0 {
        return Err(sdf::BuildError::NonPositive {
            field,
            value: length,
        }
        .into());
    }
    Ok(())
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

    // ----- Cubic -----

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

    // ----- Kelvin -----

    #[test]
    fn kelvin_accepts_positive_length() {
        let c = UnitCell::kelvin(2.4).unwrap();
        assert_eq!(c.length(), 2.4);
    }

    #[test]
    fn kelvin_rejects_zero() {
        assert!(UnitCell::kelvin(0.0).is_err());
    }

    #[test]
    fn kelvin_rejects_negative() {
        assert!(UnitCell::kelvin(-1.0).is_err());
    }

    #[test]
    fn kelvin_rejects_nan_and_inf() {
        assert!(UnitCell::kelvin(f32::NAN).is_err());
        assert!(UnitCell::kelvin(f32::INFINITY).is_err());
    }

    // ----- BccXy -----

    #[test]
    fn bccxy_accepts_positive_length() {
        let c = UnitCell::bccxy(2.4).unwrap();
        assert_eq!(c.length(), 2.4);
    }

    #[test]
    fn bccxy_rejects_zero() {
        assert!(UnitCell::bccxy(0.0).is_err());
    }

    #[test]
    fn bccxy_rejects_negative() {
        assert!(UnitCell::bccxy(-1.0).is_err());
    }

    #[test]
    fn bccxy_rejects_nan_and_inf() {
        assert!(UnitCell::bccxy(f32::NAN).is_err());
        assert!(UnitCell::bccxy(f32::INFINITY).is_err());
    }
}
