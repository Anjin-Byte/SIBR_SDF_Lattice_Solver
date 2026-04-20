//! Trait hierarchy for signed-distance functions, with precision tracked at
//! the type level.
//!
//! # Hierarchy
//!
//! ```text
//!         Sdf            (anything that maps Vec3 -> f32)
//!           |
//!           v
//!        BoundSdf        (eval(p) <= distance(p, surface), 1-Lipschitz)
//!           |
//!           v
//!        ExactSdf        (eval(p) == distance(p, surface), to f32 precision)
//! ```
//!
//! [`ExactSdf`] is a supertrait refinement of [`BoundSdf`]: an exact distance
//! is a trivial bound, so any type that impls `ExactSdf` must also impl
//! `BoundSdf`. This lets combinators accept `BoundSdf` inputs and produce
//! `BoundSdf` outputs while a separate `ExactSdf` impl adds the tighter
//! guarantee when both inputs are exact.
//!
//! # Policy: never impl both as a lie
//!
//! The trait system cannot prevent a type from implementing [`Sdf`] but
//! **claiming** to be [`BoundSdf`] when it is not. Violating this is a
//! programmer error, not a type-system error. Each primitive's doc comment
//! states which guarantees its `eval` satisfies.

use glam::Vec3;

/// A scalar function from a 3D query point to `f32`.
///
/// By itself, `Sdf` does not guarantee distance-function properties.
/// Use [`BoundSdf`] or [`ExactSdf`] to assert stronger contracts.
pub trait Sdf {
    /// Evaluates the function at `p`.
    ///
    /// # Sign convention
    ///
    /// - `eval(p) < 0.0`: `p` is inside the surface.
    /// - `eval(p) == 0.0`: `p` is on the surface (within numerical tolerance).
    /// - `eval(p) > 0.0`: `p` is outside the surface.
    fn eval(&self, p: Vec3) -> f32;
}

/// Marker trait: `eval(p) <= distance(p, surface)` and the function is
/// 1-Lipschitz, i.e. `|eval(p) - eval(q)| <= |p - q|` for all `p, q`.
///
/// Values implementing this trait are usable for sphere tracing, though the
/// trace may take smaller steps than an [`ExactSdf`] would allow.
///
/// # Implementer's contract
///
/// Types that impl `BoundSdf` must guarantee their `eval`:
///
/// 1. Never returns `NaN` for finite input.
/// 2. Returns a value `<=` the true Euclidean signed distance at `p`.
/// 3. Is 1-Lipschitz with respect to the Euclidean norm on `Vec3`.
///
/// These guarantees are not enforced by the type system; they are a trust
/// boundary. Violating them produces silently wrong sphere tracing.
pub trait BoundSdf: Sdf {}

/// Marker trait: `eval(p)` is the true Euclidean distance to the surface,
/// to within `f32` rounding.
///
/// Values implementing this trait are safe to compose via Boolean combinators
/// (union, intersection, subtraction, xor) without silently losing exactness.
///
/// Every `ExactSdf` is a [`BoundSdf`] (an exact distance is a trivial bound).
///
/// # Implementer's contract
///
/// In addition to the [`BoundSdf`] contract, types that impl `ExactSdf` must
/// guarantee their `eval` returns exactly the Euclidean signed distance, with
/// error not exceeding a few `f32` ULPs. Specifically:
///
/// 1. For any point `p` on the represented surface, `|eval(p)| < ε` for some
///    small `ε` documented per primitive.
/// 2. For any `p`, `eval(p)` has the same sign as the true signed distance.
pub trait ExactSdf: BoundSdf {}
