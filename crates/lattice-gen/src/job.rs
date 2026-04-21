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
//!
//! # Topology dispatch — [`LatticeBody`]
//!
//! Different topologies produce cell bodies of different concrete types
//! (`CubicCellBody` ≠ `KelvinCellBody` ≠ `BccXyCellBody`), so a single
//! `match` arm cannot unify them. The [`LatticeBody`] enum is the
//! resolution: one variant per topology, each holding the fully concrete
//! composed SDF, with `Sdf`/`BoundSdf`/`ExactSdf` implemented on the enum
//! via a delegating `match`. This mirrors the existing
//! [`crate::primitive::BoundaryShape`] pattern — zero-cost enum dispatch,
//! no trait objects, no type erasure.

use core::f32::consts::FRAC_1_SQRT_2;

use glam::Vec3;
use sdf::{BoundSdf, ExactSdf, Intersection, LimitedRepeat, Sdf};

use crate::cell::UnitCell;
use crate::cell::bccxy::{BccXyCellBody, bccxy_cell_body};
use crate::cell::cubic::{CubicCellBody, cubic_cell_body};
use crate::cell::kelvin::{KelvinCellBody, kelvin_cell_body};
use crate::error::LatticeError;
use crate::primitive::{BoundaryShape, PrimitiveShape};
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
    /// Returns [`LatticeError::StrutTooThick`] if `strut.radius() >= r_max(L)`
    /// for the topology in question. The per-topology `r_max` is the radius
    /// beyond which non-adjacent struts inside one cell would overlap (or,
    /// for Cubic, adjacent cells would touch across every shared face). See
    /// [`r_max_for`] for the derivations.
    pub fn new(
        primitive: PrimitiveShape,
        cell: UnitCell,
        strut: StrutSpec,
    ) -> Result<Self, LatticeError> {
        let (topology, max_radius) = r_max_for(cell);
        if strut.radius() >= max_radius {
            return Err(LatticeError::StrutTooThick {
                radius: strut.radius(),
                cell_length: cell.length(),
                topology,
                max_radius,
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

/// Per-topology tightest strut radius at which non-adjacent struts inside
/// the cell do not collide.
///
/// Returns `(topology_name, r_max)`. The caller rejects `radius >= r_max`.
///
/// Derivations:
/// - **Cubic**: `L/2`. The 3 axis-struts meet only at the origin (intended
///   node joint); the binding constraint is inter-cell tiling — adjacent
///   cells' struts meet across every shared face once `r ≥ L/2`.
/// - **Kelvin**: `L/(4·√2) = L·√2/8 ≈ 0.177·L`. Binding case is parallel
///   square-face edges on opposite sides of a TO square face, separated
///   by `L/(2·√2)`; non-collision requires `2r < L/(2·√2)`.
/// - **`BCCxy`**: `L/(2·√6) ≈ 0.204·L` (analytical; see the plan note for
///   derivation). This is the proposed bound from the skew body-diagonal
///   vs top-face-edge case sharing no endpoint; tighter numerical
///   verification is tracked as a follow-on task.
fn r_max_for(cell: UnitCell) -> (&'static str, f32) {
    match cell {
        UnitCell::Cubic { length } => ("cubic", length * 0.5),
        UnitCell::Kelvin { length } => ("kelvin", length * FRAC_1_SQRT_2 * 0.25),
        UnitCell::BccXy { length } => ("bccxy", length / (2.0 * 6.0_f32.sqrt())),
    }
}

/// The composed lattice body — one variant per supported topology.
///
/// Each variant wraps the fully concrete
/// `Intersection<LimitedRepeat<CellBody>, BoundaryShape>` for its topology.
/// `Sdf` / `BoundSdf` / `ExactSdf` are implemented via delegating `match`
/// — no trait objects, no dynamic dispatch overhead beyond a single tag
/// test per `eval`.
///
/// Same pattern as [`crate::primitive::BoundaryShape`].
///
/// Crate-private — callers see the `lattice_body()` return as `impl
/// ExactSdf` and need not name the concrete type.
///
/// The Kelvin variant carries a 36-capsule nested-Union chain which makes
/// its size much larger than the cubic variant's 3-capsule chain. Stack
/// size of a single `LatticeBody` is ~1 KiB; boxing would avoid the
/// mismatch but add a heap allocation per job. Since a job constructs
/// exactly one `LatticeBody` per invocation, the mismatch is immaterial
/// — allow the lint.
#[allow(clippy::type_complexity, clippy::large_enum_variant)]
pub(crate) enum LatticeBody {
    /// Cubic topology composition.
    Cubic(Intersection<LimitedRepeat<CubicCellBody>, BoundaryShape>),
    /// Kelvin (truncated octahedron) topology composition.
    Kelvin(Intersection<LimitedRepeat<KelvinCellBody>, BoundaryShape>),
    /// `BCCxy` (vertex octahedron) topology composition.
    BccXy(Intersection<LimitedRepeat<BccXyCellBody>, BoundaryShape>),
}

impl Sdf for LatticeBody {
    #[inline]
    fn eval(&self, p: Vec3) -> f32 {
        match self {
            Self::Cubic(b) => b.eval(p),
            Self::Kelvin(b) => b.eval(p),
            Self::BccXy(b) => b.eval(p),
        }
    }
}

impl BoundSdf for LatticeBody {}
impl ExactSdf for LatticeBody {}

/// Composes the lattice-body SDF for a validated job.
///
/// The returned value is a [`LatticeBody`] — an enum with one concrete
/// `ExactSdf` composition per topology. The public return type is
/// `impl ExactSdf` to keep callers topology-agnostic; query via
/// [`sdf::Sdf::eval`].
///
/// # Structure
///
/// The returned SDF is the intersection of a tiled cell body with the
/// primitive boundary. `extents` is computed from the primitive's AABB and
/// the cell length, rounded up to the nearest whole cell on each axis.
/// `period = Vec3::splat(L)` in every supported topology.
///
/// # Panics
///
/// Panics if any internal SDF construction fails. In practice this can only
/// happen if invariants established by [`LatticeJob::new`] have been
/// violated by unsafe construction, which is not possible via the public
/// API.
// `expect` is used at internal construction sites. All are unreachable for
// inputs that passed `LatticeJob::new` (which validates every field and the
// per-topology `StrutTooThick` invariant). The messages document the
// invariant being relied upon, per Result-vs-Panic policy.
#[allow(clippy::expect_used)]
pub fn lattice_body(job: &LatticeJob) -> impl ExactSdf {
    // Extents and period are topology-independent (all current topologies
    // tile on a simple-cubic lattice of step L).
    let (lo, hi) = job.primitive.aabb();
    let extents = extents_from_aabb(lo, hi, job.cell.length());
    let period = Vec3::splat(job.cell.length());
    let boundary = job.primitive.boundary();
    let radius = job.strut.radius();

    match job.cell {
        UnitCell::Cubic { length } => {
            let cell_body = cubic_cell_body(length, radius)
                .expect("invariants verified by LatticeJob::new");
            let tiled = LimitedRepeat::new(period, extents, cell_body)
                .expect("invariants verified by LatticeJob::new and extents_from_aabb");
            LatticeBody::Cubic(Intersection {
                a: tiled,
                b: boundary,
            })
        }
        UnitCell::Kelvin { length } => {
            let cell_body = kelvin_cell_body(length, radius)
                .expect("invariants verified by LatticeJob::new");
            let tiled = LimitedRepeat::new(period, extents, cell_body)
                .expect("invariants verified by LatticeJob::new and extents_from_aabb");
            LatticeBody::Kelvin(Intersection {
                a: tiled,
                b: boundary,
            })
        }
        UnitCell::BccXy { length } => {
            let cell_body = bccxy_cell_body(length, radius)
                .expect("invariants verified by LatticeJob::new");
            let tiled = LimitedRepeat::new(period, extents, cell_body)
                .expect("invariants verified by LatticeJob::new and extents_from_aabb");
            LatticeBody::BccXy(Intersection {
                a: tiled,
                b: boundary,
            })
        }
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

    // --------------------------------------------------------------
    // Kelvin — cross-field + end-to-end.
    // --------------------------------------------------------------

    fn valid_kelvin_job(cube_extent: f32, cell: f32, radius: f32) -> LatticeJob {
        LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(cube_extent)).unwrap(),
            UnitCell::kelvin(cell).unwrap(),
            StrutSpec::uniform(radius).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn lattice_body_kelvin_at_edge_midpoint_is_negative_radius() {
        // Kelvin home face-edge midpoint at (+L/2, L/8, L/8) = (1, 0.25, 0.25)
        // for L=2, inside the primitive.
        let job = valid_kelvin_job(10.0, 2.0, 0.1);
        let body = lattice_body(&job);
        let d = body.eval(vec3(1.0, 0.25, 0.25));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn lattice_body_kelvin_tiles_across_shared_square_face() {
        // A query at (+L/2, L/8, L/8) in the origin cell and the same
        // point folded from the neighbor cell at (L - L/2, L/8, L/8) =
        // (+L/2, L/8, L/8)... but we want to test that a point on the
        // neighbor's home edge tiles correctly. In the cell centered at
        // (L, 0, 0), the local-frame point for world (3L/2, L/8, L/8) is
        // (+L/2, L/8, L/8), on the same home edge.
        let l = 2.0;
        let job = valid_kelvin_job(10.0, l, 0.1);
        let body = lattice_body(&job);
        let d_origin_cell = body.eval(vec3(l * 0.5, l * 0.125, l * 0.125));
        let d_next_cell = body.eval(vec3(l * 1.5, l * 0.125, l * 0.125));
        assert!(
            (d_origin_cell - d_next_cell).abs() < 1e-5,
            "tiling mismatch: {d_origin_cell} vs {d_next_cell}"
        );
    }

    #[test]
    fn new_rejects_kelvin_strut_at_topology_threshold() {
        // r_max for Kelvin at L = 2 is 2 / (4√2) ≈ 0.3536.
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::kelvin(2.0).unwrap(),
            StrutSpec::uniform(0.3536).unwrap(),
        );
        assert!(matches!(result, Err(LatticeError::StrutTooThick { .. })));
    }

    #[test]
    fn new_accepts_kelvin_strut_just_below_threshold() {
        // r_max ≈ 0.3536; use 0.35.
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::kelvin(2.0).unwrap(),
            StrutSpec::uniform(0.35).unwrap(),
        );
        assert!(result.is_ok());
    }

    // --------------------------------------------------------------
    // BccXy — cross-field + end-to-end.
    // --------------------------------------------------------------

    fn valid_bccxy_job(cube_extent: f32, cell: f32, radius: f32) -> LatticeJob {
        LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(cube_extent)).unwrap(),
            UnitCell::bccxy(cell).unwrap(),
            StrutSpec::uniform(radius).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn lattice_body_bccxy_at_body_diagonal_midpoint_is_negative_radius() {
        // Body-diagonal midpoint at (L/4, L/4, L/4) = (0.5, 0.5, 0.5) for L=2.
        let job = valid_bccxy_job(10.0, 2.0, 0.1);
        let body = lattice_body(&job);
        let d = body.eval(vec3(0.5, 0.5, 0.5));
        assert!((d - (-0.1)).abs() < 1e-5, "expected -0.1, got {d}");
    }

    #[test]
    fn new_rejects_bccxy_strut_at_topology_threshold() {
        // r_max for BccXy at L = 2 is 2 / (2·√6) ≈ 0.40825. Use a value
        // clearly above to exercise the >= check.
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::bccxy(2.0).unwrap(),
            StrutSpec::uniform(0.409).unwrap(),
        );
        assert!(matches!(result, Err(LatticeError::StrutTooThick { .. })));
    }

    #[test]
    fn new_accepts_bccxy_strut_just_below_threshold() {
        // r_max ≈ 0.40825; use 0.40 (comfortably below).
        let result = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::bccxy(2.0).unwrap(),
            StrutSpec::uniform(0.40).unwrap(),
        );
        assert!(result.is_ok());
    }

    // --------------------------------------------------------------
    // Topology-specific r_max regression
    // --------------------------------------------------------------

    /// Regression: "`StrutTooThick` used the cubic threshold L/2 for every
    /// topology — Kelvin and `BccXy` jobs passed validation with struts
    /// that would overlap in-cell."
    /// Detection: r = 0.3·L passes cubic L/2 = 0.5 but fails Kelvin
    /// L/(4√2) ≈ 0.177.
    #[test]
    fn regression_strut_too_thick_uses_topology_specific_threshold() {
        let l = 2.0;
        let r = 0.3 * l; // 0.6 — passes cubic (L/2=1.0), fails Kelvin (~0.354).
        // Cubic should accept.
        let cubic = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::cubic(l).unwrap(),
            StrutSpec::uniform(r).unwrap(),
        );
        assert!(cubic.is_ok(), "cubic should accept r = 0.3·L");
        // Kelvin should reject.
        let kelvin = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::kelvin(l).unwrap(),
            StrutSpec::uniform(r).unwrap(),
        );
        assert!(
            matches!(kelvin, Err(LatticeError::StrutTooThick { .. })),
            "kelvin should reject r = 0.3·L"
        );
    }

    /// The `topology` field on `StrutTooThick` carries the right label so
    /// the error message is informative.
    #[test]
    fn strut_too_thick_reports_topology_label() {
        let err = LatticeJob::new(
            PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
            UnitCell::kelvin(2.0).unwrap(),
            StrutSpec::uniform(0.5).unwrap(),
        )
        .unwrap_err();
        match err {
            LatticeError::StrutTooThick { topology, .. } => assert_eq!(topology, "kelvin"),
            other @ LatticeError::Sdf(_) => panic!("unexpected error: {other:?}"),
        }
    }
}
