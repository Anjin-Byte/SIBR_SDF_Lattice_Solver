//! Boolean intersection — `max(a, b)`. Precision-preserving.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// The Boolean intersection of two SDFs.
///
/// Evaluates `max(a.eval(p), b.eval(p))` — the further surface constraint
/// at any query point. A point is inside the intersection iff it is inside
/// *both* operands.
///
/// # Precision tracking
///
/// - `Intersection<A: ExactSdf, B: ExactSdf>: ExactSdf`.
/// - `Intersection<A: BoundSdf, B: BoundSdf>: BoundSdf`.
///
/// # Example
///
/// ```
/// use sdf::{Intersection, Sphere, Sdf};
/// use glam::vec3;
/// let a = Sphere::new(1.0).unwrap();
/// let b = Sphere::new(2.0).unwrap();
/// // Intersection of two concentric spheres is the smaller.
/// let i = Intersection { a, b };
/// // Outside the smaller sphere is outside the intersection.
/// assert!(i.eval(vec3(1.5, 0.0, 0.0)) > 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Intersection<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: Sdf, B: Sdf> Sdf for Intersection<A, B> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        self.a.eval(p).max(self.b.eval(p))
    }
}

impl<A: BoundSdf, B: BoundSdf> BoundSdf for Intersection<A, B> {}
impl<A: ExactSdf, B: ExactSdf> ExactSdf for Intersection<A, B> {}

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

    #[test]
    fn eval_returns_max_of_inputs() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let i = Intersection { a, b };
        // At [1.5, 0, 0]: outside a (+0.5), inside b (-0.5). max = +0.5.
        assert!((i.eval(vec3(1.5, 0.0, 0.0)) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn eval_concentric_spheres_yields_smaller() {
        let small = Sphere::new(0.5).unwrap();
        let large = Sphere::new(2.0).unwrap();
        let i = Intersection { a: small, b: large };
        // Intersection of concentric spheres is the smaller one.
        // At center: both have min sign; max of two negatives is the less negative
        // one — which for concentric spheres is the smaller radius.
        assert_eq!(i.eval(vec3(0.0, 0.0, 0.0)), -0.5);
    }

    #[test]
    fn eval_idempotent_when_operands_equal() {
        let s = Sphere::new(1.0).unwrap();
        let i = Intersection { a: s, b: s };
        assert_eq!(i.eval(vec3(0.5, 0.0, 0.0)), s.eval(vec3(0.5, 0.0, 0.0)));
    }

    #[test]
    fn precision_tracking_compiles_for_exact_exact() {
        fn requires_exact<T: ExactSdf>(_: &T) {}
        let i = Intersection {
            a: Sphere::new(1.0).unwrap(),
            b: Sphere::new(2.0).unwrap(),
        };
        requires_exact(&i);
    }

    /// Regression: "Intersection used `min` instead of `max`."
    /// At a point outside a small sphere but inside a large one, intersection
    /// must report outside (positive), not inside.
    #[test]
    fn regression_max_not_min() {
        let small = Sphere::new(0.5).unwrap();
        let large = Sphere::new(2.0).unwrap();
        let i = Intersection { a: small, b: large };
        // At [1.0, 0, 0]: outside small (+0.5), inside large (-1.0).
        // max(0.5, -1.0) = 0.5 → outside the intersection.
        // min would give -1.0 → inside, which would be wrong.
        assert_eq!(i.eval(vec3(1.0, 0.0, 0.0)), 0.5);
    }

    /// Regression: "Deep intersection nesting produced wrong result."
    #[test]
    fn regression_deep_nesting_equals_single_for_identical_operands() {
        let s = Sphere::new(1.0).unwrap();
        let i2 = Intersection { a: s, b: s };
        let i4 = Intersection { a: i2, b: i2 };
        assert_eq!(i4.eval(vec3(0.5, 0.0, 0.0)), s.eval(vec3(0.5, 0.0, 0.0)));
    }
}
