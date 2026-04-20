//! Boolean subtraction — `max(-remove, from)`. Carves `remove` out of `from`.
//! Precision-preserving.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Boolean subtraction: `from` minus `remove`.
///
/// Evaluates `max(-remove.eval(p), from.eval(p))`. The result represents the
/// region that is **inside `from` and outside `remove`**.
///
/// The field names `remove` and `from` are deliberate — the operation is
/// directional and asymmetric. `Subtraction { remove: a, from: b }` reads as
/// "remove `a` from `b`", which carves `a` out of `b` and keeps the `b` shell.
/// This avoids the `a`/`b` ambiguity of the IQ-catalog `opSubtraction(a, b)`.
///
/// # Precision tracking
///
/// - `Subtraction<Remove: ExactSdf, From: ExactSdf>: ExactSdf`.
/// - `Subtraction<Remove: BoundSdf, From: BoundSdf>: BoundSdf`.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, Subtraction, Sdf};
/// use glam::vec3;
/// let outer = Sphere::new(1.0).unwrap();
/// let inner = Sphere::new(0.5).unwrap();
/// // Shell: inside outer, outside inner.
/// let shell = Subtraction { remove: inner, from: outer };
/// // At [0.75, 0, 0]: inside the shell.
/// assert!(shell.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
/// // At the center: removed by inner.
/// assert!(shell.eval(vec3(0.0, 0.0, 0.0)) > 0.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Subtraction<Remove, From> {
    /// The SDF whose interior is removed (carved out).
    pub remove: Remove,
    /// The SDF from which `remove` is carved.
    pub from: From,
}

impl<Remove: Sdf, From: Sdf> Sdf for Subtraction<Remove, From> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        (-self.remove.eval(p)).max(self.from.eval(p))
    }
}

impl<Remove: BoundSdf, From: BoundSdf> BoundSdf for Subtraction<Remove, From> {}
impl<Remove: ExactSdf, From: ExactSdf> ExactSdf for Subtraction<Remove, From> {}

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
    // --------------------------------------------------------------

    #[test]
    fn shell_inside_means_in_from_and_not_in_remove() {
        let outer = Sphere::new(1.0).unwrap();
        let inner = Sphere::new(0.5).unwrap();
        let shell = Subtraction {
            remove: inner,
            from: outer,
        };
        // Center: inside outer, inside inner → removed → outside shell.
        assert!(shell.eval(vec3(0.0, 0.0, 0.0)) > 0.0);
        // Between radii: inside outer, outside inner → inside shell.
        assert!(shell.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
        // Past outer: outside outer → outside shell.
        assert!(shell.eval(vec3(2.0, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn subtraction_of_zero_overlap_equals_from() {
        // When `remove` is outside the query point, result == from.eval(p).
        let outer = Sphere::new(1.0).unwrap();
        let tiny = Sphere::new(0.1).unwrap(); // concentric, smaller
        let s = Subtraction {
            remove: tiny,
            from: outer,
        };
        // At [0.8, 0, 0]: inside outer (-0.2), outside tiny (+0.7).
        // -remove.eval(p) = -0.7; max(-0.7, -0.2) = -0.2 = from.eval(p). ✓
        assert_eq!(s.eval(vec3(0.8, 0.0, 0.0)), outer.eval(vec3(0.8, 0.0, 0.0)));
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // --------------------------------------------------------------

    #[test]
    fn precision_tracking_compiles_for_exact_exact() {
        fn requires_exact<T: ExactSdf>(_: &T) {}
        let s = Subtraction {
            remove: Sphere::new(0.5).unwrap(),
            from: Sphere::new(1.0).unwrap(),
        };
        requires_exact(&s);
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Subtraction argument order confused — `remove` and `from`
    /// swapped, producing a shape that is the inner region *minus* the outer."
    ///
    /// The field-name API (`remove: ..., from: ...`) exists precisely to make
    /// this mistake syntactically visible at the call site. This test confirms
    /// the arithmetic matches the field semantics.
    #[test]
    fn regression_remove_from_not_swapped() {
        let outer = Sphere::new(1.0).unwrap();
        let inner = Sphere::new(0.5).unwrap();
        let correct = Subtraction {
            remove: inner,
            from: outer,
        };
        let swapped = Subtraction {
            remove: outer,
            from: inner,
        };
        // At [0.75, 0, 0] (shell point):
        // correct: inside outer, outside inner → inside shell (negative).
        // swapped: "inside inner" [but point is outside inner] ∩ "outside outer" [but it IS inside outer]
        //          → inside the inner-minus-outer (which is empty for concentric),
        //          should be OUTSIDE (positive).
        assert!(correct.eval(vec3(0.75, 0.0, 0.0)) < 0.0);
        assert!(swapped.eval(vec3(0.75, 0.0, 0.0)) > 0.0);
    }

    /// Regression: "Negation sign dropped — `max(-remove, from)` became
    /// `max(remove, from)`, which is intersection not subtraction."
    #[test]
    fn regression_negation_applied_to_remove() {
        let outer = Sphere::new(1.0).unwrap();
        let inner = Sphere::new(0.5).unwrap();
        let s = Subtraction {
            remove: inner,
            from: outer,
        };
        // At origin: outer = -1, inner = -0.5. Correct subtraction: max(-(-0.5), -1) = max(0.5, -1) = 0.5 (outside).
        // Bug: max(-0.5, -1) = -0.5 (inside), wrong.
        assert!(s.eval(vec3(0.0, 0.0, 0.0)) > 0.0);
    }
}
