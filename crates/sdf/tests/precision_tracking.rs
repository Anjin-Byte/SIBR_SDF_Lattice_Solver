//! Compile-time verification that the precision-tracking type hierarchy
//! composes correctly across the crate's public API.
//!
//! Every function below succeeds iff the trait impls on the combinators,
//! operators, and primitives preserve precision as documented. If a future
//! refactor drops an `ExactSdf` impl or weakens a precision guarantee, this
//! file will fail to compile — a compile-time regression check in the
//! spirit of *Make Invalid States Unrepresentable*.
//!
//! The crown test is [`full_cubic_lattice_composes_as_exact`]: it constructs
//! the canonical lattice-body shape that `lattice-gen` will use, and asserts
//! at compile time that it types as [`sdf::ExactSdf`]. If this fails, the
//! entire lattice-generation approach is broken at the type level.

#![allow(clippy::unwrap_used, dead_code)]

use glam::{Vec3, vec3};
use sdf::{
    AxisBox, BoundSdf, CappedCylinder, Capsule, ExactSdf, Intersection, LimitedRepeat, Repeat,
    Sphere, Subtraction, SymX, SymY, SymZ, Translate, UniformScale, Union, Xor,
};

/// Takes any exact SDF.
fn requires_exact<T: ExactSdf>(_: &T) {}

/// Takes any bound SDF (including exact).
fn requires_bound<T: BoundSdf>(_: &T) {}

// ============================================================
// Primitives
// ============================================================

#[test]
fn primitives_are_exact() {
    requires_exact(&Sphere::new(1.0).unwrap());
    requires_exact(&AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap());
    requires_exact(&Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap());
    requires_exact(&CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap());
}

#[test]
fn primitives_are_bound() {
    // Every ExactSdf is also a BoundSdf (by supertrait).
    requires_bound(&Sphere::new(1.0).unwrap());
    requires_bound(&AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap());
    requires_bound(&Capsule::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap());
    requires_bound(&CappedCylinder::new(Vec3::ZERO, vec3(1.0, 0.0, 0.0), 0.1).unwrap());
}

// ============================================================
// Combinators
// ============================================================

#[test]
fn union_of_exact_is_exact() {
    let u = Union {
        a: Sphere::new(1.0).unwrap(),
        b: AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap(),
    };
    requires_exact(&u);
}

#[test]
fn intersection_of_exact_is_exact() {
    let i = Intersection {
        a: Sphere::new(1.0).unwrap(),
        b: AxisBox::new(vec3(1.0, 1.0, 1.0)).unwrap(),
    };
    requires_exact(&i);
}

#[test]
fn subtraction_of_exact_is_exact() {
    let s = Subtraction {
        remove: Sphere::new(0.5).unwrap(),
        from: Sphere::new(1.0).unwrap(),
    };
    requires_exact(&s);
}

#[test]
fn xor_of_exact_is_exact() {
    let x = Xor {
        a: Sphere::new(0.5).unwrap(),
        b: Sphere::new(1.0).unwrap(),
    };
    requires_exact(&x);
}

// ============================================================
// Operators
// ============================================================

#[test]
fn translate_of_exact_is_exact() {
    let t = Translate {
        offset: vec3(1.0, 0.0, 0.0),
        inner: Sphere::new(1.0).unwrap(),
    };
    requires_exact(&t);
}

#[test]
fn uniform_scale_of_exact_is_exact() {
    let u = UniformScale::new(2.0, Sphere::new(1.0).unwrap()).unwrap();
    requires_exact(&u);
}

#[test]
fn sym_x_of_exact_is_exact() {
    requires_exact(&SymX {
        inner: Sphere::new(1.0).unwrap(),
    });
}

#[test]
fn sym_y_of_exact_is_exact() {
    requires_exact(&SymY {
        inner: Sphere::new(1.0).unwrap(),
    });
}

#[test]
fn sym_z_of_exact_is_exact() {
    requires_exact(&SymZ {
        inner: Sphere::new(1.0).unwrap(),
    });
}

#[test]
fn repeat_of_exact_is_exact() {
    let r = Repeat::new(Vec3::splat(2.0), Sphere::new(0.3).unwrap()).unwrap();
    requires_exact(&r);
}

#[test]
fn limited_repeat_of_exact_is_exact() {
    let lr = LimitedRepeat::new(
        Vec3::splat(2.0),
        Vec3::splat(1.0),
        Sphere::new(0.3).unwrap(),
    )
    .unwrap();
    requires_exact(&lr);
}

// ============================================================
// The crown test — the full cubic lattice body
// ============================================================

/// The canonical lattice-body composition, end-to-end.
///
/// An intersection of a bounding primitive (`AxisBox`) with a limited-repeat
/// tiling of a union of three capsule struts (x/y/z cell edges). This is
/// exactly the shape `lattice-gen` will construct for a Cubic unit cell in
/// its first implementation.
///
/// If this test fails to compile, the lattice-generation approach is
/// structurally broken at the type level and nothing downstream will work.
#[test]
fn full_cubic_lattice_composes_as_exact() {
    let l_c: f32 = 2.4; // cell length in mm
    let r: f32 = 0.264; // strut radius in mm
    let extents = Vec3::splat(5.0); // 11×11×11 cell lattice

    // The three edge-aligned struts in a cubic cell.
    let strut_x = Capsule::new(Vec3::ZERO, vec3(l_c, 0.0, 0.0), r).unwrap();
    let strut_y = Capsule::new(Vec3::ZERO, vec3(0.0, l_c, 0.0), r).unwrap();
    let strut_z = Capsule::new(Vec3::ZERO, vec3(0.0, 0.0, l_c), r).unwrap();

    let cell_body = Union {
        a: Union {
            a: strut_x,
            b: strut_y,
        },
        b: strut_z,
    };

    // Tile the cell body in 3D.
    let tiled = LimitedRepeat::new(Vec3::splat(l_c), extents, cell_body).unwrap();

    // Trim to a bounding primitive. In a real lattice-gen workload, `boundary`
    // would be the user-supplied PrimitiveShape (a sphere, cylinder, or
    // arbitrary convex primitive). For this test, an axis-aligned box suffices.
    let boundary = AxisBox::new(vec3(6.0 * l_c, 6.0 * l_c, 6.0 * l_c)).unwrap();
    let lattice = Intersection {
        a: tiled,
        b: boundary,
    };

    // Must type-check as ExactSdf.
    requires_exact(&lattice);
}

// ============================================================
// Deep composition smoke tests
// ============================================================

#[test]
fn deeply_nested_union_remains_exact() {
    // 16-deep Union nesting — the ExactSdf impl must propagate through every level.
    let s = Sphere::new(1.0).unwrap();
    let u1 = Union { a: s, b: s };
    let u2 = Union { a: u1, b: u1 };
    let u3 = Union { a: u2, b: u2 };
    let u4 = Union { a: u3, b: u3 };
    requires_exact(&u4);
}

#[test]
fn operator_stack_preserves_exactness() {
    // Stack several operators deep over an exact primitive.
    let s = Sphere::new(1.0).unwrap();
    let translated = Translate {
        offset: vec3(1.0, 0.0, 0.0),
        inner: s,
    };
    let scaled = UniformScale::new(2.0, translated).unwrap();
    let mirrored = SymX { inner: scaled };
    let tiled = Repeat::new(Vec3::splat(5.0), mirrored).unwrap();
    let bounded = LimitedRepeat::new(Vec3::splat(10.0), Vec3::splat(2.0), tiled).unwrap();
    requires_exact(&bounded);
}
