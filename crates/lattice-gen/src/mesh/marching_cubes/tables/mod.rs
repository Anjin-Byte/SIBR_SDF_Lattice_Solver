//! Marching-cubes lookup tables.
//!
//! The cube geometry (corner numbering, edge-to-corner map) is fixed by our
//! chosen convention and shared across extraction methods. Method-specific
//! tables (triangulations, ambiguity deciders) live in submodules.
//!
//! # Corner numbering
//!
//! ```text
//!      4 ---------- 5
//!     /|          /|
//!    / |         / |
//!   7 ---------- 6 |
//!   |  |        |  |
//!   |  0 -------|- 1
//!   | /         | /
//!   |/          |/
//!   3 ---------- 2
//! ```
//!
//! Corner → offset (from cell origin, in unit-cell coordinates):
//!
//! | # | Offset    |
//! |---|-----------|
//! | 0 | `(0,0,0)` |
//! | 1 | `(1,0,0)` |
//! | 2 | `(1,1,0)` |
//! | 3 | `(0,1,0)` |
//! | 4 | `(0,0,1)` |
//! | 5 | `(1,0,1)` |
//! | 6 | `(1,1,1)` |
//! | 7 | `(0,1,1)` |
//!
//! # Edge numbering
//!
//! - 0–3: bottom face (z=0), counter-clockwise viewed from +z
//! - 4–7: top face (z=1), same winding
//! - 8–11: vertical edges (0→4, 1→5, 2→6, 3→7)

pub mod classic;

/// For each corner, its integer offset from the cell origin.
pub const CORNER_OFFSETS: [(u32, u32, u32); 8] = [
    (0, 0, 0), // 0
    (1, 0, 0), // 1
    (1, 1, 0), // 2
    (0, 1, 0), // 3
    (0, 0, 1), // 4
    (1, 0, 1), // 5
    (1, 1, 1), // 6
    (0, 1, 1), // 7
];

/// For each edge, the two corner indices it connects.
/// Edge index → `(corner_a, corner_b)`.
pub const EDGE_CORNERS: [(u8, u8); 12] = [
    (0, 1), // 0: 0→1
    (1, 2), // 1: 1→2
    (2, 3), // 2: 2→3
    (3, 0), // 3: 3→0
    (4, 5), // 4: 4→5
    (5, 6), // 5: 5→6
    (6, 7), // 6: 6→7
    (7, 4), // 7: 7→4
    (0, 4), // 8: 0→4
    (1, 5), // 9: 1→5
    (2, 6), // 10: 2→6
    (3, 7), // 11: 3→7
];
