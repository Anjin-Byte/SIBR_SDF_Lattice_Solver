//! Boolean combinators — union, intersection, subtraction, and symmetric
//! difference (xor) over SDFs.
//!
//! The hard combinators ([`Union`], [`Intersection`], [`Subtraction`],
//! [`Xor`]) are **precision-preserving**: if both operands are
//! [`crate::ExactSdf`], the combined SDF is also `ExactSdf`; if both are
//! [`crate::BoundSdf`], the combined SDF is `BoundSdf`. The type system
//! enforces this at compile time — composing a `BoundSdf` into a position
//! expecting `ExactSdf` is a compile error, not a silent correctness loss.
//!
//! The smooth variants ([`SmoothUnion`], [`SmoothIntersection`]) are
//! **`Sdf`-only** by design — they do not preserve `BoundSdf` or
//! `ExactSdf`, and the compiler will refuse to use them in a position
//! that requires either marker. Smooth combinators belong at the final
//! stage before `eval`-only consumers (Marching Cubes sampling,
//! rendering); they are not safe inputs to further Boolean ops that
//! depend on bound-correctness. See the `SDF Primitives Catalog` and
//! `Three-Layer Stack` Domain Knowledge notes.
//!
//! # The SDF vs mesh distinction
//!
//! These are **SDF Boolean operations** — point-wise scalar operations with
//! cost equivalent to a single ALU instruction per query. They are not the
//! fragile, numerically-unstable "mesh Boolean" operations that operate on
//! triangle meshes. See the project's `SDF vs Mesh Booleans` Domain Knowledge
//! note.

pub mod intersection;
pub mod smooth_intersection;
pub mod smooth_union;
pub mod subtraction;
pub mod union;
pub mod xor;

pub use intersection::Intersection;
pub use smooth_intersection::SmoothIntersection;
pub use smooth_union::SmoothUnion;
pub use subtraction::Subtraction;
pub use union::Union;
pub use xor::Xor;
