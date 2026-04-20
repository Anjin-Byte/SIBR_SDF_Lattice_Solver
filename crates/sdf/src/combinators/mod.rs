//! Boolean combinators — union, intersection, subtraction, and symmetric
//! difference (xor) over SDFs.
//!
//! Each combinator is **precision-preserving**: if both operands are
//! [`crate::ExactSdf`], the combined SDF is also `ExactSdf`; if both are
//! [`crate::BoundSdf`], the combined SDF is `BoundSdf`. The type system
//! enforces this at compile time — composing a `BoundSdf` into a position
//! expecting `ExactSdf` is a compile error, not a silent correctness loss.
//!
//! # The SDF vs mesh distinction
//!
//! These are **SDF Boolean operations** — point-wise scalar operations with
//! cost equivalent to a single ALU instruction per query. They are not the
//! fragile, numerically-unstable "mesh Boolean" operations that operate on
//! triangle meshes. See the project's `SDF vs Mesh Booleans` Domain Knowledge
//! note.

pub mod intersection;
pub mod subtraction;
pub mod union;
pub mod xor;

pub use intersection::Intersection;
pub use subtraction::Subtraction;
pub use union::Union;
pub use xor::Xor;
