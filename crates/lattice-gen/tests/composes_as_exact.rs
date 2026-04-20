//! Compile-time verification that `lattice_body()`'s return type is
//! `ExactSdf` for every supported `PrimitiveShape × UnitCell` combination.
//!
//! If a future refactor drops an `ExactSdf` impl anywhere in the chain —
//! from `sdf` primitives through the cubic cell composition to the final
//! `Intersection` — this file fails to compile. A compile-time regression
//! check, in the spirit of *Make Invalid States Unrepresentable*.

#![allow(clippy::unwrap_used, dead_code)]

use glam::{Vec3, vec3};
use lattice_gen::{LatticeJob, PrimitiveShape, StrutSpec, UnitCell, lattice_body};
use sdf::{ExactSdf, Sdf};

fn requires_exact<T: ExactSdf>(_: &T) {}

#[test]
fn cube_primitive_cubic_cell_composes_as_exact() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
        UnitCell::cubic(2.0).unwrap(),
        StrutSpec::uniform(0.2).unwrap(),
    )
    .unwrap();
    let body = lattice_body(&job);
    requires_exact(&body);
}

#[test]
fn cylinder_primitive_cubic_cell_composes_as_exact() {
    let job = LatticeJob::new(
        PrimitiveShape::cylinder(Vec3::ZERO, vec3(0.0, 0.0, 10.0), 5.0).unwrap(),
        UnitCell::cubic(2.0).unwrap(),
        StrutSpec::uniform(0.2).unwrap(),
    )
    .unwrap();
    let body = lattice_body(&job);
    requires_exact(&body);
}

/// Reasonable runtime check: the end-to-end composition evaluates to finite,
/// correctly-signed values at representative points. Complements the
/// compile-time assertion above with a minimal runtime smoke test.
#[test]
fn end_to_end_smoke_runtime() {
    let job = LatticeJob::new(
        PrimitiveShape::cube(Vec3::splat(10.0)).unwrap(),
        UnitCell::cubic(2.0).unwrap(),
        StrutSpec::uniform(0.2).unwrap(),
    )
    .unwrap();
    let body = lattice_body(&job);

    // Inside a strut at the origin node.
    assert!(body.eval(Vec3::ZERO) < 0.0);
    // Strut midpoint.
    assert!((body.eval(vec3(1.0, 0.0, 0.0)) - (-0.2)).abs() < 1e-5);
    // Far outside the primitive.
    assert!(body.eval(vec3(100.0, 0.0, 0.0)) > 0.0);
    // All finite.
    assert!(body.eval(Vec3::ZERO).is_finite());
}
