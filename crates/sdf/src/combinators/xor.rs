//! Boolean symmetric difference (xor) — `max(min(a,b), -max(a,b))`.
//! Precision-preserving.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// The Boolean symmetric difference of two SDFs.
///
/// Evaluates `max(min(a, b), -max(a, b))`. The result represents the region
/// that is **inside exactly one of `a` or `b`, but not both**.
///
/// # Precision tracking
///
/// - `Xor<A: ExactSdf, B: ExactSdf>: ExactSdf`.
/// - `Xor<A: BoundSdf, B: BoundSdf>: BoundSdf`.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, Xor, Sdf};
/// use glam::vec3;
/// let a = Sphere::new(0.5).unwrap();
/// let b = Sphere::new(1.0).unwrap();
/// // Both concentric.  The "xor" region is the annular shell (inside b, outside a).
/// let shell = Xor { a, b };
/// // At the center: inside both → outside xor.
/// assert!(shell.eval(vec3(0.0, 0.0, 0.0)) > 0.0);
/// // In the shell: inside b only → inside xor.
/// assert!(shell.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Xor<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
}

impl<A: Sdf, B: Sdf> Sdf for Xor<A, B> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let da = self.a.eval(p);
        let db = self.b.eval(p);
        da.min(db).max(-da.max(db))
    }
}

impl<A: BoundSdf, B: BoundSdf> BoundSdf for Xor<A, B> {}
impl<A: ExactSdf, B: ExactSdf> ExactSdf for Xor<A, B> {}

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
    fn xor_of_concentric_spheres_is_shell() {
        let small = Sphere::new(0.5).unwrap();
        let large = Sphere::new(1.0).unwrap();
        let shell = Xor { a: small, b: large };
        // Center: inside both (both operands' SDFs negative) → xor outside.
        assert!(shell.eval(vec3(0.0, 0.0, 0.0)) > 0.0);
        // Shell interior: inside large only → xor inside.
        assert!(shell.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
        // Far outside: outside both → xor outside.
        assert!(shell.eval(vec3(2.0, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn xor_is_symmetric() {
        // xor(a, b) should equal xor(b, a) at every query point.
        let a = Sphere::new(0.5).unwrap();
        let b = Sphere::new(1.0).unwrap();
        let xab = Xor { a, b };
        let xba = Xor { a: b, b: a };
        for t in [0.0, 0.5, 0.75, 1.0, 1.5] {
            let p = vec3(t, 0.0, 0.0);
            assert_eq!(xab.eval(p), xba.eval(p));
        }
    }

    #[test]
    fn precision_tracking_compiles_for_exact_exact() {
        fn requires_exact<T: ExactSdf>(_: &T) {}
        let x = Xor {
            a: Sphere::new(1.0).unwrap(),
            b: Sphere::new(2.0).unwrap(),
        };
        requires_exact(&x);
    }

    /// Regression: "Xor produced holes — a point inside both operands was
    /// classified as inside xor rather than outside."
    /// The center of concentric spheres is the canonical test for this.
    #[test]
    fn regression_point_inside_both_is_outside_xor() {
        let a = Sphere::new(0.5).unwrap();
        let b = Sphere::new(1.0).unwrap();
        let x = Xor { a, b };
        assert!(x.eval(Vec3::ZERO) > 0.0);
    }

    /// Regression: "Xor formula lost the negation, producing `max(min, max)` =
    /// `max` (intersection) instead of `max(min, -max)`."
    #[test]
    fn regression_negation_on_max_preserved() {
        let a = Sphere::new(0.5).unwrap();
        let b = Sphere::new(1.0).unwrap();
        let x = Xor { a, b };
        // In the shell region, both da and db are negative (inside b only is
        // impossible here since they're concentric with b larger containing a).
        // Wait — at [0.75, 0, 0]: da = +0.25 (outside a), db = -0.25 (inside b).
        // min(da, db) = -0.25; max(da, db) = 0.25; -max = -0.25;
        // max(min, -max) = max(-0.25, -0.25) = -0.25 → inside xor. Correct.
        assert!(x.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
    }
}
