//! `LimitedRepeat` — bounded-extent domain repetition.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

use super::repeat::validate_period;

/// Tiles the inner SDF on a grid with per-axis period `period`, but only
/// within a bounded range of `±extents` cells from the origin.
///
/// Evaluates `inner.eval(p - period * clamp(round(p / period), -extents, extents))`.
///
/// Points whose nearest cell index is within `[-extents, extents]` behave
/// identically to [`crate::Repeat`]. Points further out are folded to the
/// *boundary* cell (the one at `±extents`) — effectively, the lattice ends
/// there and the outer region is a copy of the boundary-cell's exterior.
///
/// # Lattice relevance
///
/// This is the operator that bounds a lattice to a specific region without
/// needing a separate intersection step. Combined with
/// [`crate::Intersection`] against the primitive boundary, it produces a
/// trimmed conformal lattice body — the core of the Woodward seven-stage
/// pipeline. See the Domain Knowledge "Lattice Generation Pipeline" note.
///
/// # Caller precondition
///
/// Same as [`crate::Repeat`]: the inner SDF's non-zero region must fit inside
/// a single period-cell for the returned distance to be exact outside the
/// boundary region.
///
/// # Invariants
///
/// - Every component of `period` is strictly positive and finite.
/// - Every component of `extents` is non-negative and finite.
///
/// Extents are allowed to be zero — meaning "no repetition along this axis;
/// only the origin cell exists." Fractional extents are accepted but
/// interpreted via rounding in the clamp.
///
/// # Example
///
/// ```
/// use sdf::{Capsule, LimitedRepeat, Sdf};
/// use glam::{Vec3, vec3};
/// // A strut, tiled in a 3×3×3 lattice.
/// let strut = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
/// let bounded = LimitedRepeat::new(Vec3::splat(2.0), Vec3::splat(1.0), strut).unwrap();
/// // Inside the bounded region: tiled strut behavior.
/// assert!(bounded.eval(vec3(2.5, 0.0, 0.0)) < 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LimitedRepeat<T> {
    period: Vec3,
    extents: Vec3,
    inner: T,
}

impl<T> LimitedRepeat<T> {
    /// Constructs a `LimitedRepeat` with the given period, clamp extents,
    /// and inner SDF.
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if any period or extents component is
    ///   `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if any period component is `<= 0`.
    /// - [`BuildError::Negative`] if any extents component is `< 0`.
    pub fn new(period: Vec3, extents: Vec3, inner: T) -> Result<Self, BuildError> {
        validate_period(period)?;
        if !extents.is_finite() {
            let (field, value) = non_finite_extent(extents);
            return Err(BuildError::NonFinite { field, value });
        }
        if extents.x < 0.0 {
            return Err(BuildError::Negative {
                field: "extents.x",
                value: extents.x,
            });
        }
        if extents.y < 0.0 {
            return Err(BuildError::Negative {
                field: "extents.y",
                value: extents.y,
            });
        }
        if extents.z < 0.0 {
            return Err(BuildError::Negative {
                field: "extents.z",
                value: extents.z,
            });
        }
        Ok(Self {
            period,
            extents,
            inner,
        })
    }

    /// Returns the per-axis period.
    pub fn period(&self) -> Vec3 {
        self.period
    }

    /// Returns the per-axis clamp extents (in units of cells from origin).
    pub fn extents(&self) -> Vec3 {
        self.extents
    }

    /// Returns a reference to the inner SDF.
    pub fn inner(&self) -> &T {
        &self.inner
    }
}

fn non_finite_extent(v: Vec3) -> (&'static str, f32) {
    if !v.x.is_finite() {
        ("extents.x", v.x)
    } else if !v.y.is_finite() {
        ("extents.y", v.y)
    } else {
        ("extents.z", v.z)
    }
}

impl<T: Sdf> Sdf for LimitedRepeat<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let s = self.period;
        let l = self.extents;
        let cell = (p / s).round().clamp(-l, l);
        let q = p - s * cell;
        self.inner.eval(q)
    }
}

impl<T: BoundSdf> BoundSdf for LimitedRepeat<T> {}
impl<T: ExactSdf> ExactSdf for LimitedRepeat<T> {}

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
    use crate::{Repeat, Sphere};
    use glam::vec3;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: clamp not applied, wrong cell index at edge.
    // --------------------------------------------------------------

    #[test]
    fn within_extents_matches_unbounded_repeat() {
        // Inside the bounded region, LimitedRepeat should equal Repeat.
        let s = Sphere::new(0.3).unwrap();
        let period = Vec3::splat(2.0);
        let extents = Vec3::splat(2.0);
        let r = Repeat::new(period, s).unwrap();
        let lr = LimitedRepeat::new(period, extents, s).unwrap();
        for p in [
            Vec3::ZERO,
            vec3(2.0, 0.0, 0.0),
            vec3(4.0, 0.0, 0.0), // at the extent boundary
            vec3(-4.0, 0.0, 0.0),
            vec3(2.0, -4.0, 2.0),
        ] {
            assert_eq!(lr.eval(p), r.eval(p));
        }
    }

    #[test]
    fn outside_extents_uses_boundary_cell() {
        let s = Sphere::new(0.3).unwrap();
        let period = Vec3::splat(2.0);
        let extents = Vec3::splat(1.0);
        let lr = LimitedRepeat::new(period, extents, s).unwrap();
        // At [10, 0, 0]: clamp(round(10/2), -1, 1) = clamp(5, -1, 1) = 1
        // → query at 10 - 2 * 1 = 8 → evaluated via inner.eval([8, 0, 0]).
        assert_eq!(lr.eval(vec3(10.0, 0.0, 0.0)), s.eval(vec3(8.0, 0.0, 0.0)));
    }

    #[test]
    fn zero_extents_yields_origin_cell_only() {
        let s = Sphere::new(0.3).unwrap();
        let period = Vec3::splat(2.0);
        let extents = Vec3::ZERO;
        let lr = LimitedRepeat::new(period, extents, s).unwrap();
        // Every query folds to cell 0 → equivalent to the un-repeated inner.
        for p in [
            Vec3::ZERO,
            vec3(2.0, 0.0, 0.0),
            vec3(5.0, 0.0, 0.0),
            vec3(-5.0, 0.0, 0.0),
        ] {
            assert_eq!(lr.eval(p), s.eval(p));
        }
    }

    #[test]
    fn accessors_return_construction_values() {
        let s = Sphere::new(0.3).unwrap();
        let period = vec3(1.0, 2.0, 3.0);
        let extents = vec3(1.0, 2.0, 3.0);
        let lr = LimitedRepeat::new(period, extents, s).unwrap();
        assert_eq!(lr.period(), period);
        assert_eq!(lr.extents(), extents);
        assert_eq!(*lr.inner(), s);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_zero_period() {
        let s = Sphere::new(0.3).unwrap();
        assert!(LimitedRepeat::new(Vec3::ZERO, Vec3::splat(1.0), s).is_err());
    }

    #[test]
    fn new_rejects_negative_extents_component() {
        let s = Sphere::new(0.3).unwrap();
        assert_eq!(
            LimitedRepeat::new(Vec3::splat(2.0), vec3(-1.0, 1.0, 1.0), s),
            Err(BuildError::Negative {
                field: "extents.x",
                value: -1.0,
            })
        );
    }

    #[test]
    fn new_rejects_nan_extents_component() {
        let s = Sphere::new(0.3).unwrap();
        assert!(matches!(
            LimitedRepeat::new(Vec3::splat(2.0), vec3(f32::NAN, 1.0, 1.0), s),
            Err(BuildError::NonFinite {
                field: "extents.x",
                ..
            })
        ));
    }

    #[test]
    fn new_accepts_zero_extents() {
        let s = Sphere::new(0.3).unwrap();
        assert!(LimitedRepeat::new(Vec3::splat(2.0), Vec3::ZERO, s).is_ok());
    }

    #[test]
    fn new_accepts_mixed_extents_including_zero() {
        // extents = [0, 2, 0] → 1D tiling along y only.
        let s = Sphere::new(0.3).unwrap();
        let lr = LimitedRepeat::new(Vec3::splat(2.0), vec3(0.0, 2.0, 0.0), s).unwrap();
        // X, Z always fold to cell 0; Y tiles in [-2, 2].
        assert_eq!(lr.eval(vec3(2.0, 2.0, 2.0)), s.eval(vec3(2.0, 0.0, 2.0)));
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // --------------------------------------------------------------

    #[test]
    fn queries_far_outside_remain_finite() {
        let s = Sphere::new(0.3).unwrap();
        let lr = LimitedRepeat::new(Vec3::splat(2.0), Vec3::splat(1.0), s).unwrap();
        // Very far query: result should still be finite and positive.
        let d = lr.eval(vec3(100.0, 100.0, 100.0));
        assert!(d.is_finite());
        assert!(d > 0.0);
    }

    #[test]
    fn cells_at_extent_boundary_match_internal_cells() {
        // Points at the furthest "valid" cell should produce the same
        // evaluation as equivalent points in the fundamental cell.
        let s = Sphere::new(0.3).unwrap();
        let period = Vec3::splat(2.0);
        let extents = Vec3::splat(3.0);
        let lr = LimitedRepeat::new(period, extents, s).unwrap();
        // At +3 periods from origin: cell = 3. clamp(3, -3, 3) = 3. q = 6 - 6 = 0.
        assert_eq!(lr.eval(vec3(6.0, 0.0, 0.0)), s.eval(Vec3::ZERO));
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Clamp was not applied — `LimitedRepeat` behaved as infinite `Repeat`."
    /// Detection: at a query far outside the extents, the folded point would
    /// be in the fundamental cell for `Repeat` (always inside), but for
    /// `LimitedRepeat` it should be in the boundary cell (outside).
    #[test]
    fn regression_clamp_enforced_outside_extents() {
        let s = Sphere::new(0.1).unwrap();
        let lr = LimitedRepeat::new(Vec3::splat(2.0), Vec3::splat(1.0), s).unwrap();
        // A query at +10 on x-axis — 5 cells out, clamped to 1.
        // Folded q = 10 - 2 * 1 = 8 → outside sphere at origin by 8 - 0.1.
        let d = lr.eval(vec3(10.0, 0.0, 0.0));
        // If clamp was missing (equivalent to Repeat), q would be 0 → d = -0.1.
        assert!((d - s.eval(vec3(8.0, 0.0, 0.0))).abs() < 1e-5);
        assert!(d > 0.0, "expected outside, got {d}");
    }

    /// Regression: "Extents applied to period instead of cell index."
    #[test]
    fn regression_extents_are_cell_units_not_world_units() {
        let s = Sphere::new(0.1).unwrap();
        // extents = 1 cell → bounded region is [-period * 1, +period * 1] = [-2, 2].
        let lr = LimitedRepeat::new(Vec3::splat(2.0), Vec3::splat(1.0), s).unwrap();
        // Query at [1, 0, 0]: within extents → same as inner at [1, 0, 0].
        assert_eq!(lr.eval(vec3(1.0, 0.0, 0.0)), s.eval(vec3(1.0, 0.0, 0.0)));
        // Query at [3, 0, 0]: cell index would be round(1.5) = 2, clamped to 1.
        // q = 3 - 2 * 1 = 1. Equivalent to inner at [1, 0, 0].
        assert_eq!(lr.eval(vec3(3.0, 0.0, 0.0)), s.eval(vec3(1.0, 0.0, 0.0)));
    }

    /// Regression: "Negative extents silently accepted."
    #[test]
    fn regression_negative_extents_rejected() {
        let s = Sphere::new(0.1).unwrap();
        assert!(LimitedRepeat::new(Vec3::splat(2.0), vec3(0.0, -1.0, 0.0), s).is_err());
    }
}
