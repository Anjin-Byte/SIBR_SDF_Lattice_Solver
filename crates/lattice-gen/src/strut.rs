//! Strut specification — radius and (later) grading.
//!
//! Phase 1a supports uniform strut radius only. Functional grading (radius
//! as a function of position — see the Domain Knowledge note "Printing
//! Artifacts and Compensation") arrives in a later phase.

use crate::error::LatticeError;

/// The radius (and, in future phases, grading) of lattice struts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrutSpec {
    radius: f32,
}

impl StrutSpec {
    /// Constructs a uniform-radius strut specification.
    ///
    /// # Errors
    ///
    /// - Propagates [`LatticeError::Sdf`] from [`sdf::BuildError::NonFinite`]
    ///   if `radius` is `NaN` or `±∞`.
    /// - Propagates [`LatticeError::Sdf`] from [`sdf::BuildError::NonPositive`]
    ///   if `radius <= 0`.
    pub fn uniform(radius: f32) -> Result<Self, LatticeError> {
        if !radius.is_finite() {
            return Err(sdf::BuildError::NonFinite {
                field: "strut.radius",
                value: radius,
            }
            .into());
        }
        if radius <= 0.0 {
            return Err(sdf::BuildError::NonPositive {
                field: "strut.radius",
                value: radius,
            }
            .into());
        }
        Ok(Self { radius })
    }

    /// Returns the uniform strut radius.
    pub fn radius(self) -> f32 {
        self.radius
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
    fn uniform_accepts_positive_radius() {
        let s = StrutSpec::uniform(0.1).unwrap();
        assert_eq!(s.radius(), 0.1);
    }

    #[test]
    fn uniform_rejects_zero() {
        assert!(StrutSpec::uniform(0.0).is_err());
    }

    #[test]
    fn uniform_rejects_negative() {
        assert!(StrutSpec::uniform(-0.1).is_err());
    }

    #[test]
    fn uniform_rejects_nan() {
        assert!(StrutSpec::uniform(f32::NAN).is_err());
    }

    #[test]
    fn uniform_rejects_infinity() {
        assert!(StrutSpec::uniform(f32::INFINITY).is_err());
    }
}
