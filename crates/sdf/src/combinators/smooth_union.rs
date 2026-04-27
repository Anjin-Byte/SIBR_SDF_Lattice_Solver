//! Smooth Boolean union — polynomial smooth-min variant.
//!
//! Where [`super::Union`] uses hard `min(a, b)` and produces gradient
//! discontinuities (creases) at the locus where the two operand SDFs are
//! equal, [`SmoothUnion`] interpolates smoothly across that locus,
//! producing a fillet of approximate radius `k` at every joint between
//! the two surfaces. Cleaner-looking joints; more importantly for the
//! lattice-generator use case, **the resulting SDF is C1-continuous** so
//! the trilinear interpolant used by Marching Cubes can faithfully
//! approximate it without producing the sub-voxel "noise islands" that
//! creased SDFs sometimes generate near junction points.

use glam::Vec3;

use crate::traits::Sdf;

/// The smooth Boolean union of two SDFs.
///
/// Evaluates the polynomial smooth-min (Inigo Quilez, 2013):
///
/// ```text
/// h = clamp(0.5 + 0.5 * (b - a) / k, 0, 1)
/// smin(a, b, k) = mix(b, a, h) - k * h * (1 - h)
/// ```
///
/// where `mix(b, a, h) = b * (1 - h) + a * h`.
///
/// Properties:
/// - **At points outside the joint region** (`|a - b| > k`), the formula
///   degenerates to `min(a, b)` exactly — same as [`super::Union`].
/// - **Inside the joint region** (`|a - b| < k`), the result is up to
///   `k/4` less than `min(a, b)` — i.e., more negative than hard min.
/// - **C1-continuous** in the operand values (no gradient discontinuity).
/// - **`k = 0` falls back to bit-exact `min(a, b)`** so a nested chain of
///   `SmoothUnion`s with `k = 0` is bitwise-equivalent to a chain of
///   plain `Union`s. This is the load-bearing property that lets the
///   lattice-generator switch between hard and smooth combinators
///   without changing types.
///
/// # Precision tracking
///
/// Implements [`Sdf`] only — **not [`super::BoundSdf`] or
/// [`super::ExactSdf`]**. Smooth-min applied to two SDFs produces a
/// function whose zero-set is a *new* surface (the filleted shape), not
/// `min(a, b)`'s surface. We have no Lipschitz proof for the smoothed
/// function with respect to that new surface for arbitrary operand
/// pairs, so we conservatively expose only the `Sdf` trait. Use this
/// combinator at the final stage before `eval`-only consumers
/// (Marching Cubes sampling, rendering); never as input to a Boolean
/// op that requires `ExactSdf` or `BoundSdf`.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, SmoothUnion, Sdf};
/// use glam::vec3;
/// let a = Sphere::new(1.0).unwrap();
/// let b = Sphere::new(1.0).unwrap();
/// // Two identical unit spheres at the origin, smoothed with k=0.1.
/// let u = SmoothUnion { a, b, k: 0.1 };
/// // At the joint locus (where a.eval == b.eval, here the whole space),
/// // smooth-min returns less than min(a, b) by k/4 = 0.025.
/// let v = u.eval(vec3(0.5, 0.0, 0.0));
/// let hard = a.eval(vec3(0.5, 0.0, 0.0)); // -0.5
/// assert!(v < hard, "smooth value {v} should be < hard min {hard}");
/// assert!((v - (hard - 0.025)).abs() < 1e-5, "expected hard-0.025, got {v}");
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothUnion<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Smoothing radius. Distances within `k` of the joint locus are
    /// blended into a fillet. `k <= 0` is bit-exact equivalent to
    /// hard [`super::Union`].
    pub k: f32,
}

impl<A: Sdf, B: Sdf> Sdf for SmoothUnion<A, B> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let a = self.a.eval(p);
        let b = self.b.eval(p);
        if self.k <= 0.0 {
            // Bit-exact hard-min fallback. Critical for "k = 0 means
            // identical to Union" — used when smoothing is disabled in
            // the lattice generator.
            return a.min(b);
        }
        let h = (0.5 + 0.5 * (b - a) / self.k).clamp(0.0, 1.0);
        let mix = b * (1.0 - h) + a * h;
        mix - self.k * h * (1.0 - h)
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
    use crate::Sphere;
    use glam::vec3;

    // Polynomial smin uses ~10 f32 ops (sub, div, add, clamp, two muls,
    // sub-from-1, mul, FMA-style combine). Each op is <= 0.5 ULP; cumulative
    // worst case ~5 ULPs ≈ 6e-7 at magnitude 1. We use 1e-6 to leave headroom
    // for FMA reassociation under -O2.
    const SMOOTH_F32_TOL: f32 = 1.0e-6;

    // At points far from the joint locus, h saturates to 0 or 1, so the
    // k*h*(1-h) residual is zero exactly. Only mix vs. a/b selection drift
    // remains: 1 ULP at the operand's magnitude. At magnitude ~9.5 that is
    // ~1e-6.
    const FAR_LOCUS_TOL: f32 = 1.0e-6;

    /// k = 0 must produce the same f32 bits as hard `Union`. This is
    /// the load-bearing property that lets the lattice generator swap
    /// `Union` ↔ `SmoothUnion(k=0)` without any behavioral drift.
    #[test]
    fn k_zero_is_bit_identical_to_hard_min() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothUnion { a, b, k: 0.0 };
        for q in [
            vec3(0.0, 0.0, 0.0),
            vec3(0.5, 0.0, 0.0),
            vec3(1.5, 0.0, 0.0),
            vec3(2.5, 0.0, 0.0),
            vec3(0.3, 0.7, 0.5),
        ] {
            let s = smooth.eval(q);
            let hard = a.eval(q).min(b.eval(q));
            assert_eq!(s.to_bits(), hard.to_bits(), "k=0 drift at {q:?}");
        }
    }

    /// k < 0 also takes the hard-min fallback branch. (The branch
    /// guards `<= 0` so negative inputs are tolerated.)
    #[test]
    fn negative_k_falls_back_to_hard_min() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothUnion { a, b, k: -0.5 };
        let s = smooth.eval(vec3(1.0, 0.0, 0.0));
        let hard = a.eval(vec3(1.0, 0.0, 0.0)).min(b.eval(vec3(1.0, 0.0, 0.0)));
        assert_eq!(s.to_bits(), hard.to_bits());
    }

    /// At the joint locus where `a == b`, smooth-min returns
    /// exactly `a - k/4` (the depth of the maximum smoothing pull).
    #[test]
    fn at_equal_locus_smooth_pulls_by_k_over_four() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(1.0).unwrap();
        let smooth = SmoothUnion { a, b, k: 0.4 };
        let q = vec3(0.5, 0.0, 0.0);
        let s = smooth.eval(q);
        let hard = a.eval(q).min(b.eval(q));
        // At a == b, h = 0.5, so smooth = a - k * 0.25 = a - 0.1.
        assert!(
            (s - (hard - 0.1)).abs() < SMOOTH_F32_TOL,
            "expected hard - k/4 = {} - 0.1 = {}, got {s}",
            hard,
            hard - 0.1
        );
    }

    /// Far outside the smoothing radius (|a - b| ≫ k), the formula
    /// recovers hard `min` to within ~1 ULP.
    #[test]
    fn far_from_joint_locus_recovers_hard_min() {
        // Two very-disparate-radius spheres so |a - b| ≫ k at most points.
        let a = Sphere::new(0.5).unwrap();
        let b = Sphere::new(10.0).unwrap();
        let smooth = SmoothUnion { a, b, k: 0.1 };
        // At a point well outside both surfaces, the difference between
        // a.eval and b.eval is ~9.5 ≫ k=0.1.
        let q = vec3(20.0, 0.0, 0.0);
        let s = smooth.eval(q);
        let hard = a.eval(q).min(b.eval(q));
        assert!(
            (s - hard).abs() < FAR_LOCUS_TOL,
            "expected near-hard {hard}, got {s} (delta {})",
            (s - hard).abs()
        );
    }

    /// Symmetry: smooth-min should be symmetric under operand swap to
    /// within f32 precision. (The polynomial form is mathematically
    /// symmetric; small ULP drift from operand-order in arithmetic is
    /// acceptable.)
    #[test]
    fn smooth_min_is_approximately_symmetric() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let q = vec3(1.5, 0.0, 0.0);
        let ab = SmoothUnion { a, b, k: 0.2 }.eval(q);
        let ba = SmoothUnion { a: b, b: a, k: 0.2 }.eval(q);
        assert!(
            (ab - ba).abs() < SMOOTH_F32_TOL,
            "smooth-min symmetry broken: ab={ab}, ba={ba}"
        );
    }

    // ------------------------------------------------------------------
    // Adversarial tests — pin behavior at the awkward edges of the input
    // space so a future refactor can't silently change it.
    // ------------------------------------------------------------------

    /// NaN operand: pinned behavior is "propagates NaN through the
    /// result." Documenting it as a test stops a future refactor from
    /// silently changing the behavior to a fallback or filter.
    #[test]
    fn nan_operand_propagates() {
        struct NanSdf;
        impl Sdf for NanSdf {
            fn eval(&self, _: Vec3) -> f32 {
                f32::NAN
            }
        }
        let smooth = SmoothUnion {
            a: Sphere::new(1.0).unwrap(),
            b: NanSdf,
            k: 0.1,
        };
        assert!(smooth.eval(vec3(0.5, 0.0, 0.0)).is_nan());
    }

    /// +Inf operand: smin(a, +inf, k) yields NaN. The clamp pulls h to
    /// 1, but `mix = b*(1-h) + a*h` evaluates `∞ * 0`, which IEEE-754
    /// defines as NaN, contaminating the result. Pinning this behavior
    /// stops a future "guard against infinity" patch from silently
    /// shipping without a deliberate decision about what should happen.
    /// Callers that may produce ±Inf must handle this above the
    /// combinator (e.g., clamp the operand SDFs to a finite range).
    #[test]
    fn positive_infinity_operand_produces_nan() {
        struct InfSdf;
        impl Sdf for InfSdf {
            fn eval(&self, _: Vec3) -> f32 {
                f32::INFINITY
            }
        }
        let a = Sphere::new(1.0).unwrap();
        let smooth = SmoothUnion {
            a,
            b: InfSdf,
            k: 0.1,
        };
        assert!(smooth.eval(vec3(0.5, 0.0, 0.0)).is_nan());
    }

    /// Subnormal k (still > 0): the (b - a)/k term overflows but clamp
    /// pulls h to a saturating value; the formula must not produce NaN.
    /// At this k, smin should be indistinguishable from hard min.
    #[test]
    fn subnormal_k_does_not_nan() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothUnion { a, b, k: 1.0e-30 };
        let q = vec3(0.5, 0.0, 0.0);
        let v = smooth.eval(q);
        assert!(
            v.is_finite(),
            "smin produced non-finite for subnormal k: {v}"
        );
        let hard = a.eval(q).min(b.eval(q));
        assert!((v - hard).abs() < SMOOTH_F32_TOL);
    }

    /// k ≫ |a - b|: deep-fillet regime. h sits near 0.5; result tracks
    /// the k/4 pull at the equal-locus dual. Assert finiteness and the
    /// hard upper bound (smooth-min must never exceed hard min).
    #[test]
    fn k_much_greater_than_operand_difference() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(1.0001).unwrap();
        let smooth = SmoothUnion { a, b, k: 100.0 };
        let q = vec3(0.5, 0.0, 0.0);
        let v = smooth.eval(q);
        let hard = a.eval(q).min(b.eval(q));
        assert!(v.is_finite());
        assert!(
            v <= hard,
            "smooth-min must never exceed hard min: {v} > {hard}"
        );
    }

    // ------------------------------------------------------------------
    // Reference oracle — gold value computed in higher precision from
    // the IQ formula. This is a genuine differential check (different
    // precision, same algebra), not a self-test.
    // ------------------------------------------------------------------

    /// Gold value: `smin(0.3, 0.5, 0.4)` from the IQ polynomial form.
    /// h = clamp(0.5 + 0.5*(0.5-0.3)/0.4, 0, 1) = 0.75
    /// mix = 0.5*0.25 + 0.3*0.75 = 0.35
    /// result = 0.35 - 0.4*0.75*0.25 = 0.275
    #[test]
    fn matches_iq_polynomial_smin_gold_value() {
        struct Const(f32);
        impl Sdf for Const {
            fn eval(&self, _: Vec3) -> f32 {
                self.0
            }
        }
        let smooth = SmoothUnion {
            a: Const(0.3),
            b: Const(0.5),
            k: 0.4,
        };
        let v = smooth.eval(Vec3::ZERO);
        let expected: f32 = 0.275;
        assert!(
            (v - expected).abs() < SMOOTH_F32_TOL,
            "expected {expected}, got {v}"
        );
    }
}
