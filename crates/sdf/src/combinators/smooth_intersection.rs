//! Smooth Boolean intersection — polynomial smooth-max variant.
//!
//! Where [`super::Intersection`] uses hard `max(a, b)` and produces gradient
//! discontinuities (creases) at the locus where the two operand SDFs are
//! equal, [`SmoothIntersection`] interpolates smoothly across that locus,
//! producing an inward fillet of approximate radius `k` at every edge
//! between the two surfaces.
//!
//! # Why this matters for Marching Cubes
//!
//! At a hard `max` crease, the SDF gradient flips direction abruptly. The
//! trilinear interpolant the MC sampler uses can faithfully approximate
//! C1 functions but **cannot** faithfully approximate creased ones — at
//! the discontinuity, it samples one branch on each side of the corner
//! and linearly blends across the voxel, producing sub-voxel sign flips
//! the underlying SDF doesn't have. Those spurious flips manifest as
//! tiny isolated "noise islands" in the extracted mesh. Replacing `max`
//! with `smax(_, _, k)` for `k` ≥ a few ULPs widens the crease into a
//! C1-smooth band, eliminating the islands.

use glam::Vec3;

use crate::traits::Sdf;

/// The smooth Boolean intersection of two SDFs.
///
/// Evaluates the polynomial smooth-max (Inigo Quilez), the dual of the
/// polynomial smooth-min in [`super::SmoothUnion`]:
///
/// ```text
/// h = clamp(0.5 + 0.5 * (a - b) / k, 0, 1)
/// smax(a, b, k) = mix(b, a, h) + k * h * (1 - h)
/// ```
///
/// where `mix(b, a, h) = b * (1 - h) + a * h`. Note the `+` (vs `−` for
/// smooth-min) and the `(a - b)` (vs `(b - a)` for smooth-min).
///
/// Properties:
/// - **At points outside the joint region** (`|a - b| > k`), the formula
///   degenerates to `max(a, b)` exactly — same as [`super::Intersection`].
/// - **Inside the joint region** (`|a - b| < k`), the result is up to
///   `k/4` *greater* than `max(a, b)` — i.e., more positive, meaning the
///   body is slightly thinned at the intersection edge.
/// - **C1-continuous** in the operand values (no gradient discontinuity).
/// - **`k = 0` falls back to bit-exact `max(a, b)`** so a chain of
///   `SmoothIntersection`s with `k = 0` is bitwise-equivalent to a chain
///   of plain `Intersection`s. Load-bearing for default-behavior preservation.
///
/// # Precision tracking
///
/// Implements [`Sdf`] only — **not [`super::BoundSdf`] or
/// [`super::ExactSdf`]**. Smooth-max applied to two SDFs produces a
/// function whose zero-set is a *new* surface (the inwardly-filleted
/// shape), not `max(a, b)`'s surface. We have no Lipschitz proof for
/// the smoothed function with respect to that new surface for
/// arbitrary operand pairs, so we conservatively expose only the `Sdf`
/// trait. Use this combinator at the final stage before `eval`-only
/// consumers (Marching Cubes sampling, rendering); never as input to a
/// Boolean op that requires `ExactSdf` or `BoundSdf`.
///
/// # Example
///
/// ```
/// use sdf::{Sphere, SmoothIntersection, Sdf};
/// use glam::vec3;
/// let a = Sphere::new(1.0).unwrap();
/// let b = Sphere::new(1.0).unwrap();
/// // Two identical unit spheres; smooth-intersect them with k=0.1.
/// let i = SmoothIntersection { a, b, k: 0.1 };
/// // At the joint locus (a == b everywhere here), smax = a + k/4.
/// let v = i.eval(vec3(0.5, 0.0, 0.0));
/// let hard = a.eval(vec3(0.5, 0.0, 0.0)); // -0.5
/// assert!(v > hard, "smooth value {v} should be > hard max {hard}");
/// assert!((v - (hard + 0.025)).abs() < 1e-5, "expected hard+0.025, got {v}");
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SmoothIntersection<A, B> {
    /// First operand.
    pub a: A,
    /// Second operand.
    pub b: B,
    /// Smoothing radius. Distances within `k` of the joint locus are
    /// blended into a fillet. `k <= 0` is bit-exact equivalent to
    /// hard [`super::Intersection`].
    pub k: f32,
}

impl<A: Sdf, B: Sdf> Sdf for SmoothIntersection<A, B> {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        let a = self.a.eval(p);
        let b = self.b.eval(p);
        if self.k <= 0.0 {
            // Bit-exact hard-max fallback. Critical for "k = 0 means
            // identical to Intersection" — used when smoothing is
            // disabled in the lattice generator (the default).
            return a.max(b);
        }
        let h = (0.5 + 0.5 * (a - b) / self.k).clamp(0.0, 1.0);
        let mix = b * (1.0 - h) + a * h;
        mix + self.k * h * (1.0 - h)
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

    // Polynomial smax uses ~10 f32 ops (sub, div, add, clamp, two muls,
    // sub-from-1, mul, FMA-style combine). Each op is <= 0.5 ULP; cumulative
    // worst case ~5 ULPs ≈ 6e-7 at magnitude 1. We use 1e-6 to leave headroom
    // for FMA reassociation under -O2.
    const SMOOTH_F32_TOL: f32 = 1.0e-6;

    // At points far from the joint locus, h saturates to 0 or 1, so the
    // k*h*(1-h) residual is zero exactly. Only mix vs. a/b selection drift
    // remains: 1 ULP at the operand's magnitude. At magnitude ~9.5 that is
    // ~1e-6.
    const FAR_LOCUS_TOL: f32 = 1.0e-6;

    /// k = 0 must produce the same f32 bits as hard `Intersection`. This
    /// is the load-bearing property that lets the lattice generator swap
    /// `Intersection` ↔ `SmoothIntersection(k=0)` without any behavioral
    /// drift.
    #[test]
    fn k_zero_is_bit_identical_to_hard_max() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothIntersection { a, b, k: 0.0 };
        for q in [
            vec3(0.0, 0.0, 0.0),
            vec3(0.5, 0.0, 0.0),
            vec3(1.5, 0.0, 0.0),
            vec3(2.5, 0.0, 0.0),
            vec3(0.3, 0.7, 0.5),
        ] {
            let s = smooth.eval(q);
            let hard = a.eval(q).max(b.eval(q));
            assert_eq!(s.to_bits(), hard.to_bits(), "k=0 drift at {q:?}");
        }
    }

    /// k < 0 also takes the hard-max fallback branch.
    #[test]
    fn negative_k_falls_back_to_hard_max() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothIntersection { a, b, k: -0.5 };
        let s = smooth.eval(vec3(1.0, 0.0, 0.0));
        let hard = a.eval(vec3(1.0, 0.0, 0.0)).max(b.eval(vec3(1.0, 0.0, 0.0)));
        assert_eq!(s.to_bits(), hard.to_bits());
    }

    /// At the joint locus where `a == b`, smooth-max returns
    /// exactly `a + k/4` (the depth of the maximum smoothing push). The
    /// dual of `smin = a - k/4` for `SmoothUnion`.
    #[test]
    fn at_equal_locus_smooth_pushes_by_k_over_four() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(1.0).unwrap();
        let smooth = SmoothIntersection { a, b, k: 0.4 };
        let q = vec3(0.5, 0.0, 0.0);
        let s = smooth.eval(q);
        let hard = a.eval(q).max(b.eval(q));
        // At a == b, h = 0.5, so smax = a + k * 0.25 = a + 0.1.
        assert!(
            (s - (hard + 0.1)).abs() < SMOOTH_F32_TOL,
            "expected hard + k/4 = {} + 0.1 = {}, got {s}",
            hard,
            hard + 0.1
        );
    }

    /// Far outside the smoothing radius (|a - b| ≫ k), the formula
    /// recovers hard `max` to within ~1 ULP.
    #[test]
    fn far_from_joint_locus_recovers_hard_max() {
        // Two very-disparate-radius spheres so |a - b| ≫ k at most points.
        let a = Sphere::new(0.5).unwrap();
        let b = Sphere::new(10.0).unwrap();
        let smooth = SmoothIntersection { a, b, k: 0.1 };
        // At a point well outside both surfaces, the difference between
        // a.eval and b.eval is ~9.5 ≫ k=0.1.
        let q = vec3(20.0, 0.0, 0.0);
        let s = smooth.eval(q);
        let hard = a.eval(q).max(b.eval(q));
        assert!(
            (s - hard).abs() < FAR_LOCUS_TOL,
            "expected near-hard {hard}, got {s} (delta {})",
            (s - hard).abs()
        );
    }

    /// Symmetry: smooth-max should be symmetric under operand swap to
    /// within f32 precision.
    #[test]
    fn smooth_max_is_approximately_symmetric() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let q = vec3(1.5, 0.0, 0.0);
        let ab = SmoothIntersection { a, b, k: 0.2 }.eval(q);
        let ba = SmoothIntersection { a: b, b: a, k: 0.2 }.eval(q);
        assert!(
            (ab - ba).abs() < SMOOTH_F32_TOL,
            "smooth-max symmetry broken: ab={ab}, ba={ba}"
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
        let smooth = SmoothIntersection {
            a: Sphere::new(1.0).unwrap(),
            b: NanSdf,
            k: 0.1,
        };
        assert!(smooth.eval(vec3(0.5, 0.0, 0.0)).is_nan());
    }

    /// -Inf operand: smax(a, -inf, k) yields NaN. The clamp pulls h to
    /// 1, but `mix = b*(1-h) + a*h` evaluates `(-∞) * 0`, which
    /// IEEE-754 defines as NaN, contaminating the result. Pinning this
    /// behavior stops a future "guard against infinity" patch from
    /// silently shipping without a deliberate decision about what
    /// should happen. Callers that may produce ±Inf must handle this
    /// above the combinator (e.g., clamp the operand SDFs to a finite
    /// range).
    #[test]
    fn negative_infinity_operand_produces_nan() {
        struct NegInfSdf;
        impl Sdf for NegInfSdf {
            fn eval(&self, _: Vec3) -> f32 {
                f32::NEG_INFINITY
            }
        }
        let a = Sphere::new(1.0).unwrap();
        let smooth = SmoothIntersection {
            a,
            b: NegInfSdf,
            k: 0.1,
        };
        assert!(smooth.eval(vec3(0.5, 0.0, 0.0)).is_nan());
    }

    /// Subnormal k (still > 0): the (a - b)/k term overflows but clamp
    /// pulls h to a saturating value; the formula must not produce NaN.
    /// At this k, smax should be indistinguishable from hard max.
    #[test]
    fn subnormal_k_does_not_nan() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(2.0).unwrap();
        let smooth = SmoothIntersection { a, b, k: 1.0e-30 };
        let q = vec3(0.5, 0.0, 0.0);
        let v = smooth.eval(q);
        assert!(
            v.is_finite(),
            "smax produced non-finite for subnormal k: {v}"
        );
        let hard = a.eval(q).max(b.eval(q));
        assert!((v - hard).abs() < SMOOTH_F32_TOL);
    }

    /// k ≫ |a - b|: deep-fillet regime. h sits near 0.5; result tracks
    /// the k/4 push at the equal-locus dual. Assert finiteness and the
    /// hard lower bound (smooth-max must never drop below hard max).
    #[test]
    fn k_much_greater_than_operand_difference() {
        let a = Sphere::new(1.0).unwrap();
        let b = Sphere::new(1.0001).unwrap();
        let smooth = SmoothIntersection { a, b, k: 100.0 };
        let q = vec3(0.5, 0.0, 0.0);
        let v = smooth.eval(q);
        let hard = a.eval(q).max(b.eval(q));
        assert!(v.is_finite());
        assert!(
            v >= hard,
            "smooth-max must never drop below hard max: {v} < {hard}"
        );
    }

    // ------------------------------------------------------------------
    // Reference oracle — gold value computed in higher precision from
    // the IQ formula. This is a genuine differential check (different
    // precision, same algebra), not a self-test.
    // ------------------------------------------------------------------

    /// Gold value: `smax(0.3, 0.5, 0.4)` from the IQ polynomial form.
    /// h = clamp(0.5 + 0.5*(0.3-0.5)/0.4, 0, 1) = 0.25
    /// mix = 0.5*0.75 + 0.3*0.25 = 0.45
    /// result = 0.45 + 0.4*0.25*0.75 = 0.525
    #[test]
    fn matches_iq_polynomial_smax_gold_value() {
        struct Const(f32);
        impl Sdf for Const {
            fn eval(&self, _: Vec3) -> f32 {
                self.0
            }
        }
        let smooth = SmoothIntersection {
            a: Const(0.3),
            b: Const(0.5),
            k: 0.4,
        };
        let v = smooth.eval(Vec3::ZERO);
        let expected: f32 = 0.525;
        assert!(
            (v - expected).abs() < SMOOTH_F32_TOL,
            "expected {expected}, got {v}"
        );
    }
}
