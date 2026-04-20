//! Primitive SDFs — exact, closed-form distance functions for simple shapes.
//!
//! Each primitive:
//!
//! - Has a validated, fallible constructor that rejects degenerate inputs.
//! - Implements [`Sdf`](crate::Sdf), [`BoundSdf`](crate::BoundSdf), and
//!   [`ExactSdf`](crate::ExactSdf).
//! - Is documented with the closed-form formula it evaluates.
//! - Is covered by unit, property, and regression tests.
//!
//! See the project's Obsidian vault "Domain Knowledge / SDF Primitives Catalog"
//! for the mathematical derivations these implementations port.

pub mod axis_box;
pub mod capsule;
pub mod cylinder;
pub mod sphere;

pub use axis_box::AxisBox;
pub use capsule::Capsule;
pub use cylinder::CappedCylinder;
pub use sphere::Sphere;
