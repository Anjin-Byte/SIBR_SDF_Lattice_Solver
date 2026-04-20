//! Repeat — infinite domain repetition. The central lattice operator.

use glam::Vec3;

use crate::error::BuildError;
use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Tiles the inner SDF infinitely on a grid with per-axis period `period`.
///
/// Evaluates `inner.eval(p - period * round(p / period))`. The subtraction
/// folds the query point into the "fundamental cell" centered at the origin
/// with extents `±period / 2`.
///
/// # Lattice relevance
///
/// This is the operator that makes lattice generation feasible. Without it,
/// a lattice with N cells requires N primitive evaluations per query; with
/// `Repeat`, per-query cost is reduced to the constant number of primitives
/// in a single cell. See the Domain Knowledge note "GPU Compute Pipeline"
/// for the full analysis.
///
/// # Caller precondition
///
/// `Repeat` returns the correct SDF value **iff the inner SDF's non-zero
/// region fits inside a single period-cell**. If the inner primitive extends
/// beyond half the period in any axis, the folded point lands in the wrong
/// cell and the returned value becomes a bound, not exact.
///
/// For lattice struts (capsule with length ≤ period, radius small), this is
/// true by construction. For larger primitives, wrap the result in
/// [`crate::Intersection`] against a suitable boundary.
///
/// # Invariant
///
/// Every component of `period` is strictly positive and finite.
/// Established by [`Repeat::new`].
///
/// # Example
///
/// ```
/// use sdf::{Capsule, Repeat, Sdf};
/// use glam::{Vec3, vec3};
/// // A capsule strut, repeated every 2 units along each axis.
/// let strut = Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap();
/// let tiled = Repeat::new(Vec3::splat(2.0), strut).unwrap();
/// // Query near cell [0,0,0]: inside the strut.
/// assert!(tiled.eval(vec3(0.5, 0.0, 0.0)) < 0.0);
/// // Query near cell [2,0,0]: also inside a strut (the tiled copy).
/// assert!(tiled.eval(vec3(2.5, 0.0, 0.0)) < 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Repeat<T> {
    period: Vec3,
    inner: T,
}

impl<T> Repeat<T> {
    /// Constructs a `Repeat` with the given per-axis period and inner SDF.
    ///
    /// # Errors
    ///
    /// - [`BuildError::NonFinite`] if any period component is `NaN` or `±∞`.
    /// - [`BuildError::NonPositive`] if any period component is `<= 0`.
    pub fn new(period: Vec3, inner: T) -> Result<Self, BuildError> {
        validate_period(period)?;
        Ok(Self { period, inner })
    }

    /// Returns the per-axis period.
    pub fn period(&self) -> Vec3 {
        self.period
    }

    /// Returns a reference to the inner SDF.
    pub fn inner(&self) -> &T {
        &self.inner
    }
}

/// Validates that every component of `period` is strictly positive and finite.
/// Used by both [`Repeat`] and [`crate::LimitedRepeat`].
pub(super) fn validate_period(period: Vec3) -> Result<(), BuildError> {
    if !period.is_finite() {
        let (field, value) = non_finite_component(period);
        return Err(BuildError::NonFinite { field, value });
    }
    if period.x <= 0.0 {
        return Err(BuildError::NonPositive {
            field: "period.x",
            value: period.x,
        });
    }
    if period.y <= 0.0 {
        return Err(BuildError::NonPositive {
            field: "period.y",
            value: period.y,
        });
    }
    if period.z <= 0.0 {
        return Err(BuildError::NonPositive {
            field: "period.z",
            value: period.z,
        });
    }
    Ok(())
}

fn non_finite_component(v: Vec3) -> (&'static str, f32) {
    if !v.x.is_finite() {
        ("period.x", v.x)
    } else if !v.y.is_finite() {
        ("period.y", v.y)
    } else {
        ("period.z", v.z)
    }
}

impl<T: Sdf> Sdf for Repeat<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let s = self.period;
        let q = p - s * (p / s).round();
        self.inner.eval(q)
    }
}

impl<T: BoundSdf> BoundSdf for Repeat<T> {}
impl<T: ExactSdf> ExactSdf for Repeat<T> {}

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
    use crate::Sphere;
    use glam::vec3;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: fold wrong, wrong cell lookup, per-axis coupling bugs.
    // --------------------------------------------------------------

    #[test]
    fn origin_cell_evaluates_as_inner() {
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(Vec3::splat(2.0), s).unwrap();
        // At origin and nearby points within the origin cell, result == inner.
        for p in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(0.1, 0.1, 0.1)] {
            assert_eq!(r.eval(p), s.eval(p));
        }
    }

    #[test]
    fn adjacent_cell_matches_origin_cell_behavior() {
        // A sphere repeated every 2 units: querying at [2, 0, 0] should give
        // the same result as querying at [0, 0, 0].
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(Vec3::splat(2.0), s).unwrap();
        assert_eq!(r.eval(vec3(2.0, 0.0, 0.0)), s.eval(Vec3::ZERO));
        assert_eq!(r.eval(vec3(-2.0, 0.0, 0.0)), s.eval(Vec3::ZERO));
        assert_eq!(r.eval(vec3(4.0, 4.0, 4.0)), s.eval(Vec3::ZERO));
    }

    #[test]
    fn mid_cell_boundary_points_behave_consistently() {
        // At [1, 0, 0] (exact boundary of origin cell if period=2):
        // round(1 / 2) = round(0.5) = 0 (banker's rounding) or 1 (nearest).
        // Either way, q.x is |1 - 2*0| = 1 or |1 - 2*1| = 1 (symmetric).
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(Vec3::splat(2.0), s).unwrap();
        // Whatever the tie-breaking rule, the magnitude should equal
        // inner.eval(vec3(1, 0, 0)).
        assert!((r.eval(vec3(1.0, 0.0, 0.0)) - s.eval(vec3(1.0, 0.0, 0.0))).abs() < 1e-5);
    }

    #[test]
    fn per_axis_periods_can_differ() {
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(vec3(1.0, 2.0, 4.0), s).unwrap();
        // Each axis folds on its own period.
        assert_eq!(r.eval(vec3(1.0, 0.0, 0.0)), r.eval(Vec3::ZERO));
        assert_eq!(r.eval(vec3(0.0, 2.0, 0.0)), r.eval(Vec3::ZERO));
        assert_eq!(r.eval(vec3(0.0, 0.0, 4.0)), r.eval(Vec3::ZERO));
    }

    #[test]
    fn accessors_return_construction_values() {
        let s = Sphere::new(0.3).unwrap();
        let period = vec3(1.0, 2.0, 3.0);
        let r = Repeat::new(period, s).unwrap();
        assert_eq!(r.period(), period);
        assert_eq!(*r.inner(), s);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_zero_period_x() {
        let s = Sphere::new(0.3).unwrap();
        assert_eq!(
            Repeat::new(vec3(0.0, 1.0, 1.0), s),
            Err(BuildError::NonPositive {
                field: "period.x",
                value: 0.0,
            })
        );
    }

    #[test]
    fn new_rejects_negative_period_component() {
        let s = Sphere::new(0.3).unwrap();
        assert!(Repeat::new(vec3(1.0, -1.0, 1.0), s).is_err());
    }

    #[test]
    fn new_rejects_nan_period_component() {
        let s = Sphere::new(0.3).unwrap();
        assert!(matches!(
            Repeat::new(vec3(f32::NAN, 1.0, 1.0), s),
            Err(BuildError::NonFinite {
                field: "period.x",
                ..
            })
        ));
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: periodic Lipschitz violation, cell-boundary discontinuity.
    // --------------------------------------------------------------

    #[test]
    fn periodic_lipschitz_inside_single_cell() {
        // Within a single cell, Repeat is isometric to its inner, so Lipschitz
        // holds exactly as for the inner.
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(Vec3::splat(2.0), s).unwrap();
        // Sample pairs within the central cell [-1, 1].
        let points = [
            vec3(0.0, 0.0, 0.0),
            vec3(0.5, 0.0, 0.0),
            vec3(0.9, 0.9, 0.9),
            vec3(-0.8, 0.3, -0.5),
        ];
        for &p in &points {
            for &q in &points {
                let df = (r.eval(p) - r.eval(q)).abs();
                let dpq = (p - q).length();
                assert!(df <= dpq + 1e-5);
            }
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Repeat formula used p / s - round instead of p - s * round — wrong units."
    /// Detection: inner sphere at origin, repeated every 2 units. Query at
    /// [2, 0, 0] should produce the same result as at [0, 0, 0]. If the
    /// formula was wrong, the result would differ drastically.
    #[test]
    fn regression_fold_formula_correct() {
        let s = Sphere::new(0.3).unwrap();
        let r = Repeat::new(Vec3::splat(2.0), s).unwrap();
        assert_eq!(r.eval(vec3(2.0, 0.0, 0.0)), s.eval(Vec3::ZERO));
    }

    /// Regression: "Period components coupled across axes — folding in x
    /// used y's period or similar mixup."
    /// Detection: use distinct per-axis periods; moving along each axis by
    /// its own period should yield the same result.
    #[test]
    fn regression_per_axis_periods_independent() {
        let s = Sphere::new(0.1).unwrap();
        let r = Repeat::new(vec3(1.0, 3.0, 5.0), s).unwrap();
        // Moving one period in x:
        assert_eq!(r.eval(vec3(1.0, 0.0, 0.0)), s.eval(Vec3::ZERO));
        // Moving one period in y:
        assert_eq!(r.eval(vec3(0.0, 3.0, 0.0)), s.eval(Vec3::ZERO));
        // Moving one period in z:
        assert_eq!(r.eval(vec3(0.0, 0.0, 5.0)), s.eval(Vec3::ZERO));
    }

    /// Regression: "Zero period accepted, producing divide-by-zero at eval."
    #[test]
    fn regression_zero_period_rejected() {
        let s = Sphere::new(0.3).unwrap();
        assert!(Repeat::new(Vec3::ZERO, s).is_err());
    }
}
