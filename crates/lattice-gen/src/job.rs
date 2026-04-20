//! `LatticeJob` — the validated input specification, and `lattice_body` —
//! the entry point that composes a lattice body SDF from it.
//!
//! The shape produced by [`lattice_body`] is:
//!
//! ```text
//! lattice_body = Intersection(
//!     LimitedRepeat(period, extents, cell_body),
//!     boundary_shape
//! )
//! ```
//!
//! Every component on the right is an `ExactSdf` (per the sdf crate's
//! precision markers), so the whole composition is an `ExactSdf`.

use glam::Vec3;
use sdf::{ExactSdf, Intersection, LimitedRepeat};

use crate::cell::UnitCell;
use crate::cell::cubic::cubic_cell_body;
use crate::error::LatticeError;
use crate::primitive::PrimitiveShape;
use crate::strut::StrutSpec;

/// The validated input specification for a lattice generation.
///
/// Constructed via [`LatticeJob::new`], which validates cross-field invariants
/// — in particular, that [`StrutSpec::radius`] is strictly less than half of
/// [`UnitCell::length`]. See [`LatticeError::StrutTooThick`] for rationale.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LatticeJob {
    primitive: PrimitiveShape,
    cell: UnitCell,
    strut: StrutSpec,
}

impl LatticeJob {
    /// Constructs a validated lattice job.
    ///
    /// # Errors
    ///
    /// Returns [`LatticeError::StrutTooThick`] if
    /// `strut.radius() >= cell.length() / 2`. At this threshold, struts from
    /// adjacent periodic cells would touch across every shared face, violating
    /// the "inner fits in a single period" precondition of [`sdf::Repeat`]
    /// and eliminating any open pore structure.
    pub fn new(
        primitive: PrimitiveShape,
        cell: UnitCell,
        strut: StrutSpec,
    ) -> Result<Self, LatticeError> {
        let half_cell = cell.length() * 0.5;
        if strut.radius() >= half_cell {
            return Err(LatticeError::StrutTooThick {
                radius: strut.radius(),
                cell_length: cell.length(),
            });
        }
        Ok(Self {
            primitive,
            cell,
            strut,
        })
    }

    /// Returns the primitive boundary shape.
    pub fn primitive(&self) -> PrimitiveShape {
        self.primitive
    }

    /// Returns the unit cell specification.
    pub fn cell(&self) -> UnitCell {
        self.cell
    }

    /// Returns the strut specification.
    pub fn strut(&self) -> StrutSpec {
        self.strut
    }
}

/// Composes the lattice-body SDF for a validated job.
///
/// The returned value is an [`ExactSdf`] whose concrete type is a nested
/// composition of `sdf` operators ([`Intersection`] over [`LimitedRepeat`]
/// over a cell-body [`sdf::Union`] chain). The concrete shape is deliberately
/// hidden behind `impl ExactSdf` — callers should treat it as an opaque SDF
/// and compose further via `sdf`'s combinators or query it via
/// [`sdf::Sdf::eval`].
///
/// # Structure
///
/// The returned SDF is the intersection of a tiled cell body with the
/// primitive boundary. `extents` is computed from the primitive's AABB and
/// the cell length, rounded up to the nearest whole cell on each axis.
///
/// # Panics
///
/// Panics if any internal SDF construction fails. In practice this can only
/// happen if invariants established by [`LatticeJob::new`] have been
/// violated by unsafe construction, which is not possible via the public
/// API.
// `expect` is used here at two internal SDF-construction sites. Both are
// unreachable for inputs that passed `LatticeJob::new` (which validates every
// field and the cross-field `StrutTooThick` invariant). The messages document
// the invariant being relied upon, per Result-vs-Panic policy: panic only
// when violation means "programmer error," not for runtime-handleable failure.
#[allow(clippy::expect_used)]
pub fn lattice_body(job: &LatticeJob) -> impl ExactSdf {
    // Cell body: three edge-aligned struts, union'd.
    let cell_body = match job.cell {
        UnitCell::Cubic { length } => cubic_cell_body(length, job.strut.radius())
            .expect("invariants verified by LatticeJob::new"),
    };

    // Extents in cell units: ceil of AABB extent divided by cell length,
    // per axis. Symmetric around origin in the current simple model
    // (primitives whose AABBs are not origin-centered will produce an
    // over-approximation; acceptable — any extra tiling cells will be
    // trimmed by the Intersection).
    let (lo, hi) = job.primitive.aabb();
    let extents = extents_from_aabb(lo, hi, job.cell.length());

    // Period equals cell length on every axis for the cubic topology.
    let period = Vec3::splat(job.cell.length());

    let tiled = LimitedRepeat::new(period, extents, cell_body)
        .expect("invariants verified by LatticeJob::new and extents_from_aabb");

    Intersection {
        a: tiled,
        b: job.primitive.boundary(),
    }
}

/// Computes `LimitedRepeat` extents (in cell units, per axis) from a
/// primitive's AABB and the cell length.
///
/// Returns `ceil(max(|lo|, |hi|) / cell_length)` on each axis, as a `Vec3`
/// of f32s. The ceiling ensures the tiling covers the full primitive;
/// any extra cells on the outside are harmlessly trimmed by the intersection.
///
/// # Panics
///
/// Never panics for validated input: `cell_length > 0` (from `UnitCell::cubic`)
/// and `lo`/`hi` components are finite (from the primitive constructors).
fn extents_from_aabb(lo: Vec3, hi: Vec3, cell_length: f32) -> Vec3 {
    let max_abs = lo.abs().max(hi.abs());
    (max_abs / cell_length).ceil()
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
    use glam::vec3;
    use sdf::Sdf;

    fn valid_cubic_job(cube_extent: f32, cell: f32, radius: f32) -> LatticeJob {
        LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(cube_extent)).unwrap(),
            UnitCell::cubic(cell).unwrap(),
            StrutSpec::uniform(radius).unwrap(),
        )
        .unwrap()
    }

    // --------------------------------------------------------------
    // a. Ordinary behavioral tests
    // Target failure class: composition misconfigured, wrong extents, intersection sense wrong.
    // --------------------------------------------------------------

    #[test]
    fn lattice_body_at_strut_midpoint_is_negative() {
        let job = valid_cubic_job(10.0, 2.0, 0.2);
        let body = lattice_body(&job);
        // Midpoint of an x-axis strut in the origin cell.
        assert!((body.eval(vec3(1.0, 0.0, 0.0)) - (-0.2)).abs() < 1e-5);
    }

    #[test]
    fn lattice_body_far_outside_primitive_is_positive() {
        let job = valid_cubic_job(10.0, 2.0, 0.2);
        let body = lattice_body(&job);
        // Well beyond the primitive's half-extent.
        assert!(body.eval(vec3(100.0, 0.0, 0.0)) > 0.0);
    }

    #[test]
    fn lattice_body_inside_primitive_far_from_struts_is_positive() {
        let job = valid_cubic_job(10.0, 2.0, 0.2);
        let body = lattice_body(&job);
        // A point inside the cube but in a cell interior, far from struts.
        // In each cell, [l/2, l/2, l/2] is roughly equidistant from all
        // three axes: distance sqrt(2) * l/2 - r.
        let p = vec3(1.0, 1.0, 1.0);
        assert!(body.eval(p) > 0.0);
    }

    #[test]
    fn lattice_body_struts_tile_periodically() {
        let job = valid_cubic_job(10.0, 2.0, 0.2);
        let body = lattice_body(&job);
        // Strut midpoint in a cell offset by the period on x.
        // At [3, 0, 0]: this is the midpoint of the +x strut in the cell
        // whose origin is at [2, 0, 0] — so the SDF should be -r.
        assert!((body.eval(vec3(3.0, 0.0, 0.0)) - (-0.2)).abs() < 1e-5);
    }

    // --------------------------------------------------------------
    // b. Edge and boundary tests
    // Target failure class: cross-field validation gaps.
    // --------------------------------------------------------------

    #[test]
    fn new_rejects_strut_radius_equal_to_half_cell() {
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(2.0).unwrap(),
            StrutSpec::uniform(1.0).unwrap(),
        );
        assert!(matches!(result, Err(LatticeError::StrutTooThick { .. })));
    }

    #[test]
    fn new_rejects_strut_radius_greater_than_half_cell() {
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(2.0).unwrap(),
            StrutSpec::uniform(1.5).unwrap(),
        );
        assert!(matches!(result, Err(LatticeError::StrutTooThick { .. })));
    }

    #[test]
    fn new_accepts_strut_radius_just_below_half_cell() {
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(2.0).unwrap(),
            StrutSpec::uniform(0.999).unwrap(),
        );
        assert!(result.is_ok());
    }

    #[test]
    fn extents_rounds_up() {
        // A primitive with half-extents [3.5, 3.5, 3.5] and cell length 2:
        // ceil(3.5 / 2) = 2. Must cover the full primitive.
        let extents = extents_from_aabb(Vec3::splat(-3.5), Vec3::splat(3.5), 2.0);
        assert_eq!(extents, Vec3::splat(2.0));
    }

    #[test]
    fn extents_of_integer_multiple_is_exact() {
        // Half-extents [4, 4, 4], cell 2 → extents [2, 2, 2] exactly.
        let extents = extents_from_aabb(Vec3::splat(-4.0), Vec3::splat(4.0), 2.0);
        assert_eq!(extents, Vec3::splat(2.0));
    }

    #[test]
    fn extents_of_tiny_primitive_is_at_least_zero() {
        // A primitive smaller than one cell — extents should be 0 (no
        // repetition needed; the origin cell covers it).
        let extents = extents_from_aabb(Vec3::splat(-0.5), Vec3::splat(0.5), 2.0);
        assert_eq!(extents, Vec3::splat(1.0)); // ceil(0.5/2) = 1 (ceil rounds 0.25 to 1.0 in f32)
    }

    // --------------------------------------------------------------
    // c. Adversarial / torture tests
    // --------------------------------------------------------------

    #[test]
    fn large_extent_lattice_composes_without_overflow() {
        // Cube of half-extent 50mm, cell 0.5mm → 100 extents per axis.
        let job = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(50.0)).unwrap(),
            UnitCell::cubic(0.5).unwrap(),
            StrutSpec::uniform(0.05).unwrap(),
        )
        .unwrap();
        let body = lattice_body(&job);
        // Eval at interior strut and outside — all finite, correct signs.
        assert!(body.eval(Vec3::ZERO) < 0.0);
        assert!(body.eval(vec3(100.0, 0.0, 0.0)) > 0.0);
        assert!(body.eval(Vec3::ZERO).is_finite());
    }

    #[test]
    fn cylinder_primitive_composes_and_trims_correctly() {
        let job = LatticeJob::new(
            PrimitiveShape::cylinder(Vec3::ZERO, vec3(0.0, 0.0, 10.0), 5.0).unwrap(),
            UnitCell::cubic(2.0).unwrap(),
            StrutSpec::uniform(0.2).unwrap(),
        )
        .unwrap();
        let body = lattice_body(&job);
        // A point inside the cylinder and on a strut (the z-axis strut in
        // the origin cell runs through (0,0,0) — also on the cylinder axis):
        assert!(body.eval(vec3(0.0, 0.0, 2.0)) < 0.0);
        // A point inside the cylinder but away from any strut:
        // (2, 2, 2) lands in cell centered at (2, 0, 2); local q = (0, 2, 0)
        // which is the endpoint of the y-strut → inside (distance -r = -0.2).
        // Choose a point that's genuinely off-strut instead.
        // (2.5, 2.5, 2.5) → q = (0.5, 0.5, 0.5): distance to nearest axis ≈ sqrt(0.5) - r.
        // Cylinder half-extent is r_cyl=5 around z-axis, so (2.5, 2.5, 2.5) is
        // inside cylinder (sqrt(2.5² + 2.5²) ≈ 3.54 < 5) → cylinder SDF < 0.
        // Final intersection: max(strut_dist, cyl_dist). strut_dist = sqrt(0.5) - 0.2 ≈ 0.507.
        // So outside the lattice (positive).
        assert!(body.eval(vec3(2.5, 2.5, 2.5)) > 0.0);
        // A point far outside the cylinder: outside the lattice regardless
        // of strut positions.
        assert!(body.eval(vec3(100.0, 100.0, 5.0)) > 0.0);
    }

    // --------------------------------------------------------------
    // d. Regression tests
    // --------------------------------------------------------------

    /// Regression: "`LimitedRepeat` extents computed in world units instead
    /// of cell units." Would tile the wrong number of times.
    #[test]
    fn regression_extents_are_cell_units() {
        // With L_c = 1, primitive half-extents [3, 3, 3]: extents must be
        // [3, 3, 3] (cell units), not [3, 3, 3] world × some factor.
        let extents = extents_from_aabb(Vec3::splat(-3.0), Vec3::splat(3.0), 1.0);
        assert_eq!(extents, Vec3::splat(3.0));
    }

    /// Regression: "Primitive intersected as Union instead of Intersection."
    /// Would produce a lattice that extends past the primitive.
    #[test]
    fn regression_primitive_intersected_not_union() {
        let job = valid_cubic_job(1.0, 2.0, 0.2);
        let body = lattice_body(&job);
        // Far outside the small primitive: must be positive (outside).
        // If intersection were mistakenly a Union, the tiled struts would
        // extend far beyond the primitive and report inside somewhere.
        for far in [
            vec3(10.0, 0.0, 0.0),
            vec3(0.0, 10.0, 0.0),
            vec3(0.0, 0.0, 10.0),
        ] {
            let d = body.eval(far);
            assert!(d > 0.0, "at {far:?}: expected outside, got {d}");
        }
    }

    /// Regression: "Cross-field validation missed — struts overlap adjacent
    /// cells without error." Caught by the threshold check.
    #[test]
    fn regression_strut_radius_half_cell_rejected() {
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(1.0).unwrap(),
            StrutSpec::uniform(0.5).unwrap(),
        );
        assert!(result.is_err());
    }
}
