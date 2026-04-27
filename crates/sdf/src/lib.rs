//! Signed-distance functions with type-level precision tracking.
//!
//! # Overview
//!
//! An **SDF** is a function `f: R^3 -> R` that, for a given query point `p`,
//! returns a scalar whose sign indicates inside/outside relative to some
//! surface and whose magnitude relates to the distance from `p` to that
//! surface. This crate provides:
//!
//! - A trait hierarchy ([`Sdf`], [`BoundSdf`], [`ExactSdf`]) that distinguishes
//!   exact distance functions from conservative bounds **at the type level**.
//! - Primitive SDFs (sphere, box, capsule, cylinder) with validated
//!   constructors.
//! - Boolean combinators (`Union`, `Intersection`, `Subtraction`, `Xor`) that
//!   preserve precision when their inputs do.
//! - Operators (`Translate`, `UniformScale`, `SymX/Y/Z`, `Repeat`,
//!   `LimitedRepeat`) that transform the query frame while preserving precision.
//!
//! # Precision
//!
//! The crate enforces the difference between:
//!
//! - [`ExactSdf`] — `eval(p)` is the true Euclidean distance to the surface,
//!   to within float32 rounding. 1-Lipschitz. Safe inside Boolean combinators.
//! - [`BoundSdf`] — `eval(p) <= distance(p, surface)`. Usable for sphere tracing
//!   with a step discount. **Not** safe inside Boolean combinators (produces
//!   wrong distances silently).
//!
//! Any `ExactSdf` is also a `BoundSdf` (an exact distance is a trivial bound).
//! Combinators and operators that lose exactness — smooth unions, most
//! deformations — produce `BoundSdf` rather than `ExactSdf`. The compiler
//! rejects the composition of a `BoundSdf` into a position that requires
//! `ExactSdf`, so contamination is a compile error, not a silent runtime bug.
//!
//! # Non-goals
//!
//! - No GPU types. Zero `wgpu` dependency. WGSL ports of these primitives
//!   live in the `gpu` crate under `shaders/`; cross-crate differential tests
//!   verify CPU and GPU agree.
//! - No lattice-specific logic. Unit cells and the seven-stage generation
//!   pipeline live in `lattice-gen`.
//! - No `async`, no interior mutability, no globals. `eval` is a pure function.

#![forbid(unsafe_code)]

pub mod combinators;
pub mod error;
pub mod operators;
pub mod primitives;
pub mod traits;

pub use combinators::{Intersection, SmoothIntersection, SmoothUnion, Subtraction, Union, Xor};
pub use error::BuildError;
pub use operators::{LimitedRepeat, Repeat, SymX, SymY, SymZ, Translate, UniformScale};
pub use primitives::{AxisBox, CappedCylinder, Capsule, Sphere};
pub use traits::{BoundSdf, ExactSdf, Sdf};
