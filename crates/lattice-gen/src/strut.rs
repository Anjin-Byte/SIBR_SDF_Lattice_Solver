//! Strut specification — radius and (later) grading.
//!
//! Phase 1a supports uniform strut radius only. Functional grading (radius
//! as a function of position — see the Domain Knowledge note "Printing
//! Artifacts and Compensation") arrives in a later phase.

use crate::error::LatticeError;

/// The radius (and, in future phases, grading) of lattice struts, plus
/// the smoothness of strut-to-strut joints.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StrutSpec {
    radius: f32,
    joint_smoothness: f32,
}

impl StrutSpec {
    /// Constructs a uniform-radius strut specification with sharp
    /// (hard-`min`) joints — i.e., `joint_smoothness = 0`.
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
        Ok(Self {
            radius,
            joint_smoothness: 0.0,
        })
    }

    /// Sets the joint smoothness (smooth-`min` radius) for strut-to-strut
    /// junctions. `0.0` means hard-`min` (creased joints, current
    /// default behavior); positive values fillet the joints with a
    /// blend radius `~k`. Cleaner-looking junctions and avoids the
    /// sub-voxel "noise islands" Marching Cubes can produce near
    /// gradient discontinuities.
    ///
    /// Practical guidance: pick `k` between `0.1 * strut_radius` and
    /// `0.5 * strut_radius`. Smaller `k` → minimally-visible smoothing,
    /// already enough to suppress noise islands at typical extraction
    /// resolutions. Larger `k` → visibly rounded joints, which can be
    /// desirable as a stress-concentration mitigation in printed parts.
    ///
    /// # Errors
    ///
    /// - [`LatticeError::Sdf`] from [`sdf::BuildError::NonFinite`] if
    ///   `k` is `NaN` or `±∞`.
    /// - [`LatticeError::Sdf`] from [`sdf::BuildError::NonPositive`] if
    ///   `k < 0`. (Zero is allowed — it means hard `min`.)
    pub fn with_joint_smoothness(mut self, k: f32) -> Result<Self, LatticeError> {
        if !k.is_finite() {
            return Err(sdf::BuildError::NonFinite {
                field: "strut.joint_smoothness",
                value: k,
            }
            .into());
        }
        if k < 0.0 {
            return Err(sdf::BuildError::NonPositive {
                field: "strut.joint_smoothness",
                value: k,
            }
            .into());
        }
        self.joint_smoothness = k;
        Ok(self)
    }

    /// Returns the uniform strut radius.
    pub fn radius(self) -> f32 {
        self.radius
    }

    /// Returns the joint-smoothing radius. `0.0` means hard `min`.
    pub fn joint_smoothness(self) -> f32 {
        self.joint_smoothness
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
