//! Boolean union — `min(a, b)`. Precision-preserving.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// The Boolean union of two SDFs.
///
/// Evaluates `min(a.eval(p), b.eval(p))` — the closer of the two surfaces
/// at any query point. This is the correct distance to the union iff both
/// inputs satisfy the `ExactSdf` or `BoundSdf` contract; for arbitrary `Sdf`
/// inputs (e.g., after a bound-producing deformation), the result is still
/// a bound but without stronger guarantees.
///
/// # Precision tracking
///
/// - `Union<A: ExactSdf, B: ExactSdf>: ExactSdf` — safe inside further Boolean ops.
/// - `Union<A: BoundSdf, B: BoundSdf>: BoundSdf` — safe for sphere tracing.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, Union, Sdf};
/// use glam::vec3;
/// let a = Sphere::new(1.0).unwrap();
/// let b = Sphere::new(1.0).unwrap();
/// let u = Union { a, b };  // two coincident spheres
/// // At the shared center, inside both.
/// assert!(u.eval(vec3(0.0, 0.0, 0.0)) < 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Union<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: Sdf, B: Sdf> Sdf for Union<A, B> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        self.a.eval(p).min(self.b.eval(p))
    }
}

impl<A: BoundSdf, B: BoundSdf> BoundSdf for Union<A, B> {}
impl<A: ExactSdf, B: ExactSdf> ExactSdf for Union<A, B> {}

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
    // Target failure class: min vs max confusion, wrong combinator semantics.
    // --------------------------------------------------------------

    #[test]
    fn eval_returns_min_of_inputs() {
        // Two disjoint spheres on the x-axis, at [0] and [3].
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(1.0).unwrap();
        // Note: we construct Union over two Sphere values directly; one is at
        // the origin (implicit from Sphere's center being origin), the other
        // we'd normally translate. Since we have no Translate yet in Session 2,
        // we use a point that reveals min behavior with two centered spheres:
        // at the origin, both return -1.0, and min(-1, -1) = -1.
        let u = Union { a, b };
        assert_eq!(u.eval(vec3(0.0, 0.0, 0.0)), -1.0);
        // Outside both (at distance 2): min(1, 1) = 1.
        assert_eq!(u.eval(vec3(2.0, 0.0, 0.0)), 1.0);
    }

    /// Distinguishes union from intersection by choosing operands whose
    /// interior/exterior classifications disagree at the query point.
    /// A point inside only A: union says inside, intersection says outside.
    #[test]
    fn eval_distinguishes_union_from_intersection_via_mixed_radii() {
        let small = Sphere::new(0.5).unwrap();
        let large = Sphere::new(2.0).unwrap();
        // At [1.0, 0, 0]: outside small (d=0.5), inside large (d=-1.0).
        // Union (min): -1.0 — inside.
        let u = Union { a: small, b: large };
        assert_eq!(u.eval(vec3(1.0, 0.0, 0.0)), -1.0);
    }

    #[test]
    fn eval_idempotent_when_operands_equal() {
        let s = Sphere::new(1.0).unwrap();
        let u = Union { a: s, b: s };
        // Since both operands are identical, union evaluates same as either.
        assert_eq!(u.eval(vec3(0.5, 0.0, 0.0)), s.eval(vec3(0.5, 0.0, 0.0)));
        assert_eq!(u.eval(vec3(2.0, 0.0, 0.0)), s.eval(vec3(2.0, 0.0, 0.0)));
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: identity-input handling, extreme operand combinations.
    // --------------------------------------------------------------

    #[test]
    fn precision_tracking_compiles_for_exact_exact() {
        fn requires_exact<T: ExactSdf>(_: &T) {}
        let u = Union {
            a: Sphere::new(1.0).unwrap(),
            b: Sphere::new(2.0).unwrap(),
        };
        requires_exact(&u); // Must compile
    }

    #[test]
    fn precision_tracking_compiles_for_bound_bound() {
        fn requires_bound<T: BoundSdf>(_: &T) {}
        let u = Union {
            a: Sphere::new(1.0).unwrap(),
            b: Sphere::new(2.0).unwrap(),
        };
        requires_bound(&u); // Must compile (Sphere is Exact → also Bound)
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: deep nesting stack overflow, mixed-precision bugs.
    // --------------------------------------------------------------

    #[test]
    fn deep_nested_union_evaluates_correctly() {
        // Nested unions of identical spheres — result must equal the sphere.
        let s = Sphere::new(1.0).unwrap();
        let u2 = Union { a: s, b: s };
        let u4 = Union { a: u2, b: u2 };
        let u8 = Union { a: u4, b: u4 };
        assert_eq!(u8.eval(vec3(0.5, 0.0, 0.0)), s.eval(vec3(0.5, 0.0, 0.0)));
    }

    // --------------------------------------------------------------
    // d. Regression tests (seeded on plausible bugs)
    // --------------------------------------------------------------

    /// Regression: "Union used `max` instead of `min`."
    /// For two coincident spheres at the center, `min` gives -r (inside);
    /// `max` would also give -r in this specific case, so we use disjoint
    /// regions: large sphere minus a small sphere, query a point that's
    /// inside large and outside small.
    #[test]
    fn regression_min_not_max() {
        let small = Sphere::new(0.5).unwrap();
        let large = Sphere::new(2.0).unwrap();
        let u = Union { a: small, b: large };
        // At [1, 0, 0]: outside small (+0.5), inside large (-1.0).
        // min(0.5, -1.0) = -1.0 → inside the union.
        // max(0.5, -1.0) = 0.5 → outside, which would be wrong for union.
        assert_eq!(u.eval(vec3(1.0, 0.0, 0.0)), -1.0);
    }

    /// Regression: "Union did not preserve exactness — the `ExactSdf` impl was
    /// accidentally dropped during a refactor."
    /// This test asserts at compile time that exactness survives.
    #[test]
    fn regression_exactness_preserved() {
        fn requires_exact<T: ExactSdf>(_: &T) {}
        let u = Union {
            a: Sphere::new(1.0).unwrap(),
            b: Sphere::new(2.0).unwrap(),
        };
        requires_exact(&u);
    }
}
