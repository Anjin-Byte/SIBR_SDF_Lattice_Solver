//! Translate — rigid translation of the inner SDF's frame.

use glam::Vec3;

use crate::traits::{BoundSdf, ExactSdf, Sdf};

/// Shifts the inner SDF by `offset`. Evaluates `inner.eval(p - offset)`.
///
/// Translation is a rigid isometry, so distance is preserved exactly.
/// Precision is inherited from the inner SDF: `Translate<ExactSdf>` is exact,
/// `Translate<BoundSdf>` is bound.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, Translate, Sdf};
/// use glam::vec3;
/// let s = Sphere::new(1.0).unwrap();
/// let shifted = Translate {
///     offset: vec3(3.0, 0.0, 0.0),
///     inner: s,
/// };
/// // Center of the translated sphere is now at [3, 0, 0].
/// assert_eq!(shifted.eval(vec3(3.0, 0.0, 0.0)), -1.0);
/// // Original sphere had center at origin; translated sphere does not.
/// assert_eq!(shifted.eval(vec3(0.0, 0.0, 0.0)), 2.0);
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Translate<T> {
    /// The translation offset.
    pub offset: Vec3,
    /// The inner SDF being translated.
    pub inner: T,
}

impl<T: Sdf> Sdf for Translate<T> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        self.inner.eval(p - self.offset)
    }
}

impl<T: BoundSdf> BoundSdf for Translate<T> {}
impl<T: ExactSdf> ExactSdf for Translate<T> {}

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
    use proptest::prelude::*;

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: wrong sign on offset, inner not delegated.
    // --------------------------------------------------------------

    #[test]
    fn zero_offset_equals_inner_sdf() {
        let s = Sphere::new(1.0).unwrap();
        let t = Translate {
            offset: Vec3::ZERO,
            inner: s,
        };
        for p in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(2.0, 1.0, 0.0)] {
            assert_eq!(t.eval(p), s.eval(p));
        }
    }

    #[test]
    fn nonzero_offset_shifts_center() {
        let s = Sphere::new(1.0).unwrap();
        let offset = vec3(3.0, 0.0, 0.0);
        let t = Translate { offset, inner: s };
        // Querying at the new center returns -radius.
        assert_eq!(t.eval(offset), -1.0);
        // Querying at the original center returns distance to the shifted sphere.
        assert_eq!(t.eval(Vec3::ZERO), offset.length() - 1.0);
    }

    #[test]
    fn evaluate_shifted_plus_offset_equals_inner_plus_zero() {
        // For any p, Translate{offset, inner}.eval(offset + q) == inner.eval(q).
        let s = Sphere::new(1.0).unwrap();
        let offset = vec3(2.0, -1.0, 0.5);
        let t = Translate { offset, inner: s };
        for q in [Vec3::ZERO, vec3(0.5, 0.0, 0.0), vec3(2.0, 1.0, -0.5)] {
            assert_eq!(t.eval(offset + q), s.eval(q));
        }
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: identity behavior.
    // --------------------------------------------------------------

    #[test]
    fn zero_offset_is_identity() {
        let s = Sphere::new(2.5).unwrap();
        let t = Translate {
            offset: Vec3::ZERO,
            inner: s,
        };
        assert_eq!(t.eval(vec3(1.0, 1.0, 1.0)), s.eval(vec3(1.0, 1.0, 1.0)));
    }

    #[test]
    fn double_translate_composes_additively() {
        let s = Sphere::new(1.0).unwrap();
        let t1 = Translate {
            offset: vec3(1.0, 0.0, 0.0),
            inner: s,
        };
        let t2 = Translate {
            offset: vec3(2.0, 0.0, 0.0),
            inner: t1,
        };
        // Net translation should be [3, 0, 0].
        assert_eq!(t2.eval(vec3(3.0, 0.0, 0.0)), -1.0);
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // Target failure class: Lipschitz violation under translation (should never happen).
    // --------------------------------------------------------------

    proptest! {
        /// Translation preserves Lipschitz exactly (it's an isometry).
        #[test]
        fn prop_translation_preserves_lipschitz(
            r in 1e-3f32..1e3f32,
            ofs in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
            pc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
            qc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
        ) {
            let t = Translate {
                offset: vec3(ofs.0, ofs.1, ofs.2),
                inner: Sphere::new(r).unwrap(),
            };
            let p = vec3(pc.0, pc.1, pc.2);
            let q = vec3(qc.0, qc.1, qc.2);
            let df = (t.eval(p) - t.eval(q)).abs();
            let dpq = (p - q).length();
            prop_assert!(df <= dpq + 1e-4);
        }

        /// Translating by `d` and querying at `p + d` must equal querying
        /// the untranslated inner at `p`, up to float32 rounding in the
        /// round-trip `(p + offset) - offset` that is not present in `p` alone.
        #[test]
        fn prop_delegation_equivalent(
            r in 1e-3f32..1e3f32,
            ofs in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
            pc in (-100.0f32..100.0, -100.0f32..100.0, -100.0f32..100.0),
        ) {
            let s = Sphere::new(r).unwrap();
            let offset = vec3(ofs.0, ofs.1, ofs.2);
            let t = Translate { offset, inner: s };
            let p = vec3(pc.0, pc.1, pc.2);
            // Slack accounts for f32 rounding in `p + offset - offset`.
            // With values up to ~100, the worst-case rounding is a few ULPs
            // of magnitude ~100, i.e. ~1e-5.
            let slack = 1e-4 * (p.length() + offset.length()).max(1.0);
            prop_assert!((t.eval(p + offset) - s.eval(p)).abs() <= slack);
        }
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "Sign of offset flipped — implementation added instead of subtracted."
    /// `Translate { offset, inner }.eval(offset)` should return the inner SDF
    /// at the origin (i.e., at the inner's center), NOT at `2 * offset`.
    #[test]
    fn regression_offset_sign_correct() {
        let s = Sphere::new(1.0).unwrap();
        let offset = vec3(3.0, 0.0, 0.0);
        let t = Translate { offset, inner: s };
        // Correct: t.eval(offset) = s.eval(offset - offset) = s.eval(ZERO) = -1.
        // Wrong sign would give s.eval(2*offset) ≈ 5.
        assert_eq!(t.eval(offset), -1.0);
    }

    /// Regression: "Translate applied twice produced wrong composite offset."
    #[test]
    fn regression_translate_composition() {
        let s = Sphere::new(1.0).unwrap();
        let t1 = Translate {
            offset: vec3(1.0, 0.0, 0.0),
            inner: s,
        };
        let t2 = Translate {
            offset: vec3(0.0, 2.0, 0.0),
            inner: t1,
        };
        // Net center at [1, 2, 0].
        assert_eq!(t2.eval(vec3(1.0, 2.0, 0.0)), -1.0);
    }
}
