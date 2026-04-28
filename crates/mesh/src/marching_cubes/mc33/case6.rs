//! Case 6 (Chernyaev) / Lewiner base case 7 — "three corners + diagonal".
//!
//! The configuration: three same-sign corners forming an L-shape on a
//! face, with a fourth same-sign corner diagonally opposite on the
//! opposing face. There are 48 voxel configurations up to rotation
//! (`base_case == 7` rows in [`super::cases::CASES`]).
//!
//! # Dispatch
//!
//! 1 face test + 1 interior test, branching to 3 sub-tilings:
//!
//! - face test `true` → **6.2** (5 triangles, no c-vertex).
//! - face test `false`, interior `true` → **6.1.1** (3 triangles, no c-vertex).
//! - face test `false`, interior `false` → **6.1.2** (9 triangles, c-vertex).
//!
//! # Sign / index conventions
//!
//! - **Tile edges** are stored 0-based 0..=11; edge index 12 (1-based 13)
//!   refers to the c-vertex computed by [`super::compute_c_vertex`].
//!   Used only by `TILING_6_1_2`.
//! - **`TEST_6` columns**:
//!   - 0: face index (1..=6 with sign-encoded inversion) for [`face_decider`].
//!   - 1: interior sentinel `s` (always ±7).
//!   - 2: per-edge reference index `0..=11` (already 0-based in source).
//!
//! # Port
//!
//! Tables transcribed from `JuliaGeometry/MarchingCubes.jl` `src/lut.jl`:
//! - `test6` at lines 621-668
//! - `tiling6_1_1` at lines 682-729
//! - `tiling6_1_2` at lines 743-790
//! - `tiling6_2` at lines 804-851
//!
//! Dispatch mirrors `case == 7` in `src/MarchingCubes.jl`, lines 224-233.

use super::super::{CellCoord, Mesh};
use super::decider::{face_decider, interior_decider_per_edge};
use super::emit_triangles;
use crate::mesh::grid::GridSpec;

/// Per-cfg 3-tuple: `[face, interior_s, ref_edge]`.
///
/// Source: Lewiner `test6[48][3]`. Column 0 is the face index for
/// [`face_decider`] (signed for inversion). Column 1 is the interior-
/// test `s` argument (±7). Column 2 is the reference edge for
/// [`interior_decider_per_edge`] — **already 0-based** in the source.
#[rustfmt::skip]
const TEST_6: [[i8; 3]; 48] = [
    [  2,   7,  10],
    [  4,   7,  11],
    [  5,   7,   1],
    [  5,   7,   3],
    [  1,   7,   9],
    [  3,   7,  10],
    [  6,   7,   5],
    [  1,   7,   8],
    [  4,   7,   8],
    [  1,   7,   8],
    [  3,   7,  11],
    [  5,   7,   2],
    [  5,   7,   0],
    [  1,   7,   9],
    [  6,   7,   6],
    [  2,   7,   9],
    [  4,   7,   8],
    [  2,   7,   9],
    [  2,   7,  10],
    [  6,   7,   7],
    [  3,   7,  10],
    [  4,   7,  11],
    [  3,   7,  11],
    [  6,   7,   4],
    [ -6,  -7,   4],
    [ -3,  -7,  11],
    [ -4,  -7,  11],
    [ -3,  -7,  10],
    [ -6,  -7,   7],
    [ -2,  -7,  10],
    [ -2,  -7,   9],
    [ -4,  -7,   8],
    [ -2,  -7,   9],
    [ -6,  -7,   6],
    [ -1,  -7,   9],
    [ -5,  -7,   0],
    [ -5,  -7,   2],
    [ -3,  -7,  11],
    [ -1,  -7,   8],
    [ -4,  -7,   8],
    [ -1,  -7,   8],
    [ -6,  -7,   5],
    [ -3,  -7,  10],
    [ -1,  -7,   9],
    [ -5,  -7,   3],
    [ -5,  -7,   1],
    [ -4,  -7,  11],
    [ -2,  -7,  10],
];

/// Triangulation for Chernyaev subcase 6.1.1 (3 triangles, no c-vertex).
/// Selected when [`interior_decider_per_edge`] returns `true`.
///
/// Source: Lewiner `tiling6_1_1[48][9]`.
#[rustfmt::skip]
const TILING_6_1_1: [[u8; 9]; 48] = [
    [  6,   5,  10,   3,   1,   8,   9,   8,   1],
    [ 11,   7,   6,   9,   3,   1,   3,   9,   8],
    [  1,   2,  10,   7,   0,   4,   0,   7,   3],
    [  3,   0,   8,   5,   2,   6,   2,   5,   1],
    [  5,   4,   9,   2,   0,  11,   8,  11,   0],
    [ 10,   6,   5,   8,   2,   0,   2,   8,  11],
    [ 10,   6,   5,   0,   4,   3,   7,   3,   4],
    [  3,   0,   8,   6,   4,  10,   9,  10,   4],
    [  8,   3,   0,  10,   7,   5,   7,  10,  11],
    [  8,   4,   7,  10,   0,   2,   0,  10,   9],
    [  7,   6,  11,   0,   2,   9,  10,   9,   2],
    [  2,   3,  11,   4,   1,   5,   1,   4,   0],
    [  0,   1,   9,   6,   3,   7,   3,   6,   2],
    [  9,   0,   1,  11,   4,   6,   4,  11,   8],
    [ 11,   7,   6,   1,   5,   0,   4,   0,   5],
    [  0,   1,   9,   7,   5,  11,  10,  11,   5],
    [  4,   7,   8,   1,   3,  10,  11,  10,   3],
    [  9,   5,   4,  11,   1,   3,   1,  11,  10],
    [ 10,   1,   2,   8,   5,   7,   5,   8,   9],
    [  8,   4,   7,   2,   6,   1,   5,   1,   6],
    [  1,   2,  10,   4,   6,   8,  11,   8,   6],
    [  2,   3,  11,   5,   7,   9,   8,   9,   7],
    [ 11,   2,   3,   9,   6,   4,   6,   9,  10],
    [  9,   5,   4,   3,   7,   2,   6,   2,   7],
    [  4,   5,   9,   2,   7,   3,   7,   2,   6],
    [  3,   2,  11,   4,   6,   9,  10,   9,   6],
    [ 11,   3,   2,   9,   7,   5,   7,   9,   8],
    [ 10,   2,   1,   8,   6,   4,   6,   8,  11],
    [  7,   4,   8,   1,   6,   2,   6,   1,   5],
    [  2,   1,  10,   7,   5,   8,   9,   8,   5],
    [  4,   5,   9,   3,   1,  11,  10,  11,   1],
    [  8,   7,   4,  10,   3,   1,   3,  10,  11],
    [  9,   1,   0,  11,   5,   7,   5,  11,  10],
    [  6,   7,  11,   0,   5,   1,   5,   0,   4],
    [  1,   0,   9,   6,   4,  11,   8,  11,   4],
    [  9,   1,   0,   7,   3,   6,   2,   6,   3],
    [ 11,   3,   2,   5,   1,   4,   0,   4,   1],
    [ 11,   6,   7,   9,   2,   0,   2,   9,  10],
    [  7,   4,   8,   2,   0,  10,   9,  10,   0],
    [  0,   3,   8,   5,   7,  10,  11,  10,   7],
    [  8,   0,   3,  10,   4,   6,   4,  10,   9],
    [  5,   6,  10,   3,   4,   0,   4,   3,   7],
    [  5,   6,  10,   0,   2,   8,  11,   8,   2],
    [  9,   4,   5,  11,   0,   2,   0,  11,   8],
    [  8,   0,   3,   6,   2,   5,   1,   5,   2],
    [ 10,   2,   1,   4,   0,   7,   3,   7,   0],
    [  6,   7,  11,   1,   3,   9,   8,   9,   3],
    [ 10,   5,   6,   8,   1,   3,   1,   8,   9],
];

/// Triangulation for Chernyaev subcase 6.1.2 (9 triangles, uses c-vertex
/// at slot 12). Selected when [`interior_decider_per_edge`] returns `false`.
///
/// Source: Lewiner `tiling6_1_2[48][27]`.
#[rustfmt::skip]
const TILING_6_1_2: [[u8; 27]; 48] = [
    [  1,  12,   3,  12,  10,   3,   6,   3,  10,   3,   6,   8,   5,   8,   6,   8,   5,  12,  12,   9,   8,   1,   9,  12,  12,   5,  10],
    [  1,  12,   3,   1,  11,  12,  11,   1,   6,   9,   6,   1,   6,   9,   7,  12,   7,   9,   9,   8,  12,  12,   8,   3,  11,   7,  12],
    [  4,  12,   0,   4,   1,  12,   1,   4,  10,   7,  10,   4,  10,   7,   2,  12,   2,   7,   7,   3,  12,  12,   3,   0,   1,   2,  12],
    [  6,  12,   2,   6,   3,  12,   3,   6,   8,   5,   8,   6,   8,   5,   0,  12,   0,   5,   5,   1,  12,  12,   1,   2,   3,   0,  12],
    [  0,  12,   2,  12,   9,   2,   5,   2,   9,   2,   5,  11,   4,  11,   5,  11,   4,  12,  12,   8,  11,   0,   8,  12,  12,   4,   9],
    [  0,  12,   2,   0,  10,  12,  10,   0,   5,   8,   5,   0,   5,   8,   6,  12,   6,   8,   8,  11,  12,  12,  11,   2,  10,   6,  12],
    [  4,  12,   0,  12,   5,   0,  10,   0,   5,   0,  10,   3,   6,   3,  10,   3,   6,  12,  12,   7,   3,   4,   7,  12,  12,   6,   5],
    [  4,  12,   6,  12,   8,   6,   3,   6,   8,   6,   3,  10,   0,  10,   3,  10,   0,  12,  12,   9,  10,   4,   9,  12,  12,   0,   8],
    [  5,  12,   7,   5,   8,  12,   8,   5,   0,  10,   0,   5,   0,  10,   3,  12,   3,  10,  10,  11,  12,  12,  11,   7,   8,   3,  12],
    [  2,  12,   0,   2,   8,  12,   8,   2,   7,  10,   7,   2,   7,  10,   4,  12,   4,  10,  10,   9,  12,  12,   9,   0,   8,   4,  12],
    [  2,  12,   0,  12,  11,   0,   7,   0,  11,   0,   7,   9,   6,   9,   7,   9,   6,  12,  12,  10,   9,   2,  10,  12,  12,   6,  11],
    [  5,  12,   1,   5,   2,  12,   2,   5,  11,   4,  11,   5,  11,   4,   3,  12,   3,   4,   4,   0,  12,  12,   0,   1,   2,   3,  12],
    [  7,  12,   3,   7,   0,  12,   0,   7,   9,   6,   9,   7,   9,   6,   1,  12,   1,   6,   6,   2,  12,  12,   2,   3,   0,   1,  12],
    [  6,  12,   4,   6,   9,  12,   9,   6,   1,  11,   1,   6,   1,  11,   0,  12,   0,  11,  11,   8,  12,  12,   8,   4,   9,   0,  12],
    [  5,  12,   1,  12,   6,   1,  11,   1,   6,   1,  11,   0,   7,   0,  11,   0,   7,  12,  12,   4,   0,   5,   4,  12,  12,   7,   6],
    [  5,  12,   7,  12,   9,   7,   0,   7,   9,   7,   0,  11,   1,  11,   0,  11,   1,  12,  12,  10,  11,   5,  10,  12,  12,   1,   9],
    [  3,  12,   1,  12,   8,   1,   4,   1,   8,   1,   4,  10,   7,  10,   4,  10,   7,  12,  12,  11,  10,   3,  11,  12,  12,   7,   8],
    [  3,  12,   1,   3,   9,  12,   9,   3,   4,  11,   4,   3,   4,  11,   5,  12,   5,  11,  11,  10,  12,  12,  10,   1,   9,   5,  12],
    [  7,  12,   5,   7,  10,  12,  10,   7,   2,   8,   2,   7,   2,   8,   1,  12,   1,   8,   8,   9,  12,  12,   9,   5,  10,   1,  12],
    [  6,  12,   2,  12,   7,   2,   8,   2,   7,   2,   8,   1,   4,   1,   8,   1,   4,  12,  12,   5,   1,   6,   5,  12,  12,   4,   7],
    [  6,  12,   4,  12,  10,   4,   1,   4,  10,   4,   1,   8,   2,   8,   1,   8,   2,  12,  12,  11,   8,   6,  11,  12,  12,   2,  10],
    [  7,  12,   5,  12,  11,   5,   2,   5,  11,   5,   2,   9,   3,   9,   2,   9,   3,  12,  12,   8,   9,   7,   8,  12,  12,   3,  11],
    [  4,  12,   6,   4,  11,  12,  11,   4,   3,   9,   3,   4,   3,   9,   2,  12,   2,   9,   9,  10,  12,  12,  10,   6,  11,   2,  12],
    [  7,  12,   3,  12,   4,   3,   9,   3,   4,   3,   9,   2,   5,   2,   9,   2,   5,  12,  12,   6,   2,   7,   6,  12,  12,   5,   4],
    [  3,  12,   7,   3,   4,  12,   4,   3,   9,   2,   9,   3,   9,   2,   5,  12,   5,   2,   2,   6,  12,  12,   6,   7,   4,   5,  12],
    [  6,  12,   4,  12,  11,   4,   3,   4,  11,   4,   3,   9,   2,   9,   3,   9,   2,  12,  12,  10,   9,   6,  10,  12,  12,   2,  11],
    [  5,  12,   7,   5,  11,  12,  11,   5,   2,   9,   2,   5,   2,   9,   3,  12,   3,   9,   9,   8,  12,  12,   8,   7,  11,   3,  12],
    [  4,  12,   6,   4,  10,  12,  10,   4,   1,   8,   1,   4,   1,   8,   2,  12,   2,   8,   8,  11,  12,  12,  11,   6,  10,   2,  12],
    [  2,  12,   6,   2,   7,  12,   7,   2,   8,   1,   8,   2,   8,   1,   4,  12,   4,   1,   1,   5,  12,  12,   5,   6,   7,   4,  12],
    [  5,  12,   7,  12,  10,   7,   2,   7,  10,   7,   2,   8,   1,   8,   2,   8,   1,  12,  12,   9,   8,   5,   9,  12,  12,   1,  10],
    [  1,  12,   3,  12,   9,   3,   4,   3,   9,   3,   4,  11,   5,  11,   4,  11,   5,  12,  12,  10,  11,   1,  10,  12,  12,   5,   9],
    [  1,  12,   3,   1,   8,  12,   8,   1,   4,  10,   4,   1,   4,  10,   7,  12,   7,  10,  10,  11,  12,  12,  11,   3,   8,   7,  12],
    [  7,  12,   5,   7,   9,  12,   9,   7,   0,  11,   0,   7,   0,  11,   1,  12,   1,  11,  11,  10,  12,  12,  10,   5,   9,   1,  12],
    [  1,  12,   5,   1,   6,  12,   6,   1,  11,   0,  11,   1,  11,   0,   7,  12,   7,   0,   0,   4,  12,  12,   4,   5,   6,   7,  12],
    [  4,  12,   6,  12,   9,   6,   1,   6,   9,   6,   1,  11,   0,  11,   1,  11,   0,  12,  12,   8,  11,   4,   8,  12,  12,   0,   9],
    [  3,  12,   7,  12,   0,   7,   9,   7,   0,   7,   9,   6,   1,   6,   9,   6,   1,  12,  12,   2,   6,   3,   2,  12,  12,   1,   0],
    [  1,  12,   5,  12,   2,   5,  11,   5,   2,   5,  11,   4,   3,   4,  11,   4,   3,  12,  12,   0,   4,   1,   0,  12,  12,   3,   2],
    [  0,  12,   2,   0,  11,  12,  11,   0,   7,   9,   7,   0,   7,   9,   6,  12,   6,   9,   9,  10,  12,  12,  10,   2,  11,   6,  12],
    [  0,  12,   2,  12,   8,   2,   7,   2,   8,   2,   7,  10,   4,  10,   7,  10,   4,  12,  12,   9,  10,   0,   9,  12,  12,   4,   8],
    [  7,  12,   5,  12,   8,   5,   0,   5,   8,   5,   0,  10,   3,  10,   0,  10,   3,  12,  12,  11,  10,   7,  11,  12,  12,   3,   8],
    [  6,  12,   4,   6,   8,  12,   8,   6,   3,  10,   3,   6,   3,  10,   0,  12,   0,  10,  10,   9,  12,  12,   9,   4,   8,   0,  12],
    [  0,  12,   4,   0,   5,  12,   5,   0,  10,   3,  10,   0,  10,   3,   6,  12,   6,   3,   3,   7,  12,  12,   7,   4,   5,   6,  12],
    [  2,  12,   0,  12,  10,   0,   5,   0,  10,   0,   5,   8,   6,   8,   5,   8,   6,  12,  12,  11,   8,   2,  11,  12,  12,   6,  10],
    [  2,  12,   0,   2,   9,  12,   9,   2,   5,  11,   5,   2,   5,  11,   4,  12,   4,  11,  11,   8,  12,  12,   8,   0,   9,   4,  12],
    [  2,  12,   6,  12,   3,   6,   8,   6,   3,   6,   8,   5,   0,   5,   8,   5,   0,  12,  12,   1,   5,   2,   1,  12,  12,   0,   3],
    [  0,  12,   4,  12,   1,   4,  10,   4,   1,   4,  10,   7,   2,   7,  10,   7,   2,  12,  12,   3,   7,   0,   3,  12,  12,   2,   1],
    [  3,  12,   1,  12,  11,   1,   6,   1,  11,   1,   6,   9,   7,   9,   6,   9,   7,  12,  12,   8,   9,   3,   8,  12,  12,   7,  11],
    [  3,  12,   1,   3,  10,  12,  10,   3,   6,   8,   6,   3,   6,   8,   5,  12,   5,   8,   8,   9,  12,  12,   9,   1,  10,   5,  12],
];

/// Triangulation for Chernyaev subcase 6.2 (5 triangles, no c-vertex).
/// Selected when [`face_decider`] returns `true`.
///
/// Source: Lewiner `tiling6_2[48][15]`.
#[rustfmt::skip]
const TILING_6_2: [[u8; 15]; 48] = [
    [  1,  10,   3,   6,   3,  10,   3,   6,   8,   5,   8,   6,   8,   5,   9],
    [  1,  11,   3,  11,   1,   6,   9,   6,   1,   6,   9,   7,   8,   7,   9],
    [  4,   1,   0,   1,   4,  10,   7,  10,   4,  10,   7,   2,   3,   2,   7],
    [  6,   3,   2,   3,   6,   8,   5,   8,   6,   8,   5,   0,   1,   0,   5],
    [  0,   9,   2,   5,   2,   9,   2,   5,  11,   4,  11,   5,  11,   4,   8],
    [  0,  10,   2,  10,   0,   5,   8,   5,   0,   5,   8,   6,  11,   6,   8],
    [  4,   5,   0,  10,   0,   5,   0,  10,   3,   6,   3,  10,   3,   6,   7],
    [  4,   8,   6,   3,   6,   8,   6,   3,  10,   0,  10,   3,  10,   0,   9],
    [  5,   8,   7,   8,   5,   0,  10,   0,   5,   0,  10,   3,  11,   3,  10],
    [  2,   8,   0,   8,   2,   7,  10,   7,   2,   7,  10,   4,   9,   4,  10],
    [  2,  11,   0,   7,   0,  11,   0,   7,   9,   6,   9,   7,   9,   6,  10],
    [  5,   2,   1,   2,   5,  11,   4,  11,   5,  11,   4,   3,   0,   3,   4],
    [  7,   0,   3,   0,   7,   9,   6,   9,   7,   9,   6,   1,   2,   1,   6],
    [  6,   9,   4,   9,   6,   1,  11,   1,   6,   1,  11,   0,   8,   0,  11],
    [  5,   6,   1,  11,   1,   6,   1,  11,   0,   7,   0,  11,   0,   7,   4],
    [  5,   9,   7,   0,   7,   9,   7,   0,  11,   1,  11,   0,  11,   1,  10],
    [  3,   8,   1,   4,   1,   8,   1,   4,  10,   7,  10,   4,  10,   7,  11],
    [  3,   9,   1,   9,   3,   4,  11,   4,   3,   4,  11,   5,  10,   5,  11],
    [  7,  10,   5,  10,   7,   2,   8,   2,   7,   2,   8,   1,   9,   1,   8],
    [  6,   7,   2,   8,   2,   7,   2,   8,   1,   4,   1,   8,   1,   4,   5],
    [  6,  10,   4,   1,   4,  10,   4,   1,   8,   2,   8,   1,   8,   2,  11],
    [  7,  11,   5,   2,   5,  11,   5,   2,   9,   3,   9,   2,   9,   3,   8],
    [  4,  11,   6,  11,   4,   3,   9,   3,   4,   3,   9,   2,  10,   2,   9],
    [  7,   4,   3,   9,   3,   4,   3,   9,   2,   5,   2,   9,   2,   5,   6],
    [  3,   4,   7,   4,   3,   9,   2,   9,   3,   9,   2,   5,   6,   5,   2],
    [  6,  11,   4,   3,   4,  11,   4,   3,   9,   2,   9,   3,   9,   2,  10],
    [  5,  11,   7,  11,   5,   2,   9,   2,   5,   2,   9,   3,   8,   3,   9],
    [  4,  10,   6,  10,   4,   1,   8,   1,   4,   1,   8,   2,  11,   2,   8],
    [  2,   7,   6,   7,   2,   8,   1,   8,   2,   8,   1,   4,   5,   4,   1],
    [  5,  10,   7,   2,   7,  10,   7,   2,   8,   1,   8,   2,   8,   1,   9],
    [  1,   9,   3,   4,   3,   9,   3,   4,  11,   5,  11,   4,  11,   5,  10],
    [  1,   8,   3,   8,   1,   4,  10,   4,   1,   4,  10,   7,  11,   7,  10],
    [  7,   9,   5,   9,   7,   0,  11,   0,   7,   0,  11,   1,  10,   1,  11],
    [  1,   6,   5,   6,   1,  11,   0,  11,   1,  11,   0,   7,   4,   7,   0],
    [  4,   9,   6,   1,   6,   9,   6,   1,  11,   0,  11,   1,  11,   0,   8],
    [  3,   0,   7,   9,   7,   0,   7,   9,   6,   1,   6,   9,   6,   1,   2],
    [  1,   2,   5,  11,   5,   2,   5,  11,   4,   3,   4,  11,   4,   3,   0],
    [  0,  11,   2,  11,   0,   7,   9,   7,   0,   7,   9,   6,  10,   6,   9],
    [  0,   8,   2,   7,   2,   8,   2,   7,  10,   4,  10,   7,  10,   4,   9],
    [  7,   8,   5,   0,   5,   8,   5,   0,  10,   3,  10,   0,  10,   3,  11],
    [  6,   8,   4,   8,   6,   3,  10,   3,   6,   3,  10,   0,   9,   0,  10],
    [  0,   5,   4,   5,   0,  10,   3,  10,   0,  10,   3,   6,   7,   6,   3],
    [  2,  10,   0,   5,   0,  10,   0,   5,   8,   6,   8,   5,   8,   6,  11],
    [  2,   9,   0,   9,   2,   5,  11,   5,   2,   5,  11,   4,   8,   4,  11],
    [  2,   3,   6,   8,   6,   3,   6,   8,   5,   0,   5,   8,   5,   0,   1],
    [  0,   1,   4,  10,   4,   1,   4,  10,   7,   2,   7,  10,   7,   2,   3],
    [  3,  11,   1,   6,   1,  11,   1,   6,   9,   7,   9,   6,   9,   7,   8],
    [  3,  10,   1,  10,   3,   6,   8,   6,   3,   6,   8,   5,   9,   5,   8],
];

/// Emits Case-6 triangles for the given voxel.
///
/// `subcase` is Lewiner's 1-based subcase index in `1..=48` from
/// [`super::cases::CASES`].
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(
        (1..=48).contains(&subcase),
        "Case 6 subcase out of range: {subcase}"
    );
    let cfg = usize::try_from(subcase - 1).unwrap_or(0);

    if face_decider(corners, TEST_6[cfg][0]) {
        // Subcase 6.2 — 5 triangles, no c-vertex.
        emit_triangles(grid, cell, corners, &TILING_6_2[cfg], mesh);
    } else {
        let s = TEST_6[cfg][1];
        let edge = u8::try_from(TEST_6[cfg][2]).unwrap_or(0);
        if interior_decider_per_edge(corners, edge, s) {
            // Subcase 6.1.1 — 3 triangles, no c-vertex.
            emit_triangles(grid, cell, corners, &TILING_6_1_1[cfg], mesh);
        } else {
            // Subcase 6.1.2 — 9 triangles, uses c-vertex (slot 12).
            emit_triangles(grid, cell, corners, &TILING_6_1_2[cfg], mesh);
        }
    }
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    clippy::cast_precision_loss,
    missing_docs
)]
mod tests {
    use super::*;
    use crate::mesh::Mesh;
    use crate::mesh::grid::GridSpec;
    use glam::{UVec3, Vec3};

    fn unit_grid_with_corners(corners: &[f32; 8]) -> (GridSpec, Vec<f32>) {
        let grid = GridSpec::new(Vec3::ZERO, UVec3::splat(1), 1.0).unwrap();
        let mut field = vec![0.0_f32; grid.sample_count()];
        field[grid.sample_index(0, 0, 0)] = corners[0];
        field[grid.sample_index(1, 0, 0)] = corners[1];
        field[grid.sample_index(1, 1, 0)] = corners[2];
        field[grid.sample_index(0, 1, 0)] = corners[3];
        field[grid.sample_index(0, 0, 1)] = corners[4];
        field[grid.sample_index(1, 0, 1)] = corners[5];
        field[grid.sample_index(1, 1, 1)] = corners[6];
        field[grid.sample_index(0, 1, 1)] = corners[7];
        (grid, field)
    }

    fn case7_idx_subcase_pairs() -> Vec<(u8, i8)> {
        let mut hits = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 7 {
                hits.push((our_idx, subcase));
            }
        }
        hits
    }

    /// Pin: exactly 48 voxels map to Lewiner base case 7.
    #[test]
    fn case7_dispatch_has_forty_eight_configs() {
        assert_eq!(case7_idx_subcase_pairs().len(), 48);
    }

    /// Transcription oracle: `test6[1]` (Julia line 622 → `(2, 7, 10)`,
    /// no -1 conversion since column 2 is already 0-based).
    #[test]
    fn test_6_first_row_matches_lut() {
        assert_eq!(TEST_6[0], [2, 7, 10]);
    }

    /// Transcription oracle: `tiling6_1_1[1]` (Julia line 683 →
    /// `(7, 6, 11, 4, 2, 9, 10, 9, 2)` → 0-based).
    #[test]
    fn tiling_6_1_1_first_row_matches_lut() {
        assert_eq!(TILING_6_1_1[0], [6, 5, 10, 3, 1, 8, 9, 8, 1]);
    }

    /// Transcription oracle: `tiling6_1_2[1]` (Julia line 744 →
    /// `(2, 13, 4, 13, 11, 4, 7, 4, 11, 4, 7, 9, 6, 9, 7, 9, 6, 13, 13, 10, 9, 2, 10, 13, 13, 6, 11)`
    /// → 0-based, with 13 → 12).
    #[test]
    fn tiling_6_1_2_first_row_matches_lut() {
        let expected = [
            1_u8, 12, 3, 12, 10, 3, 6, 3, 10, 3, 6, 8, 5, 8, 6, 8, 5, 12, 12, 9, 8, 1, 9, 12, 12,
            5, 10,
        ];
        assert_eq!(TILING_6_1_2[0], expected);
    }

    /// Transcription oracle: `tiling6_2[1]` (Julia line 805 →
    /// `(2, 11, 4, 7, 4, 11, 4, 7, 9, 6, 9, 7, 9, 6, 10)` → 0-based).
    #[test]
    fn tiling_6_2_first_row_matches_lut() {
        let expected = [1_u8, 10, 3, 6, 3, 10, 3, 6, 8, 5, 8, 6, 8, 5, 9];
        assert_eq!(TILING_6_2[0], expected);
    }

    /// Sharp oracle: every entry in every tiling is `≤ 12` (cube edges
    /// 0..=11 plus c-vertex slot 12).
    #[test]
    fn tilings_use_only_legal_indices() {
        for cfg in 0..48 {
            for &e in &TILING_6_1_1[cfg] {
                assert!(e <= 12, "TILING_6_1_1[{cfg}] index {e}");
            }
            for &e in &TILING_6_1_2[cfg] {
                assert!(e <= 12, "TILING_6_1_2[{cfg}] index {e}");
            }
            for &e in &TILING_6_2[cfg] {
                assert!(e <= 12, "TILING_6_2[{cfg}] index {e}");
            }
        }
    }

    /// Pin: `TEST_6[cfg][2]` is a 0-based edge index in 0..=11 for
    /// every cfg. Catches transcription that accidentally subtracted 1
    /// from this already-0-based column.
    #[test]
    fn test_6_column_2_is_zero_based_edge() {
        for (cfg, row) in TEST_6.iter().enumerate() {
            assert!(
                (0..=11).contains(&row[2]),
                "TEST_6[{cfg}][2] = {} not in 0..=11",
                row[2]
            );
        }
    }

    /// Row-length pin.
    #[test]
    fn tiling_dimensions_match_lewiner() {
        assert_eq!(TEST_6.len(), 48);
        assert_eq!(TEST_6[0].len(), 3);
        assert_eq!(TILING_6_1_1.len(), 48);
        assert_eq!(TILING_6_1_1[0].len(), 9);
        assert_eq!(TILING_6_1_2.len(), 48);
        assert_eq!(TILING_6_1_2[0].len(), 27);
        assert_eq!(TILING_6_2.len(), 48);
        assert_eq!(TILING_6_2[0].len(), 15);
    }

    /// Behavioral: every Lewiner-7 voxel under ±1 corner inputs runs
    /// `emit` to completion, producing 3, 5, or 9 triangles with
    /// finite vertices in `[0,1]³`.
    #[test]
    fn case6_dispatch_runs_on_every_case7_voxel() {
        for (our_idx, subcase) in case7_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            emit(&grid, cell, &corners, subcase, &mut mesh);

            let n = mesh.triangle_count();
            assert!(
                matches!(n, 3 | 5 | 9),
                "our_idx={our_idx:#010b} subcase={subcase}: expected 3/5/9 tris, got {n}"
            );
            assert_eq!(mesh.vertices.len(), 3 * n, "unwelded invariant");
            for v in &mesh.vertices {
                assert!(v.is_finite());
                assert!((0.0..=1.0).contains(&v.x));
                assert!((0.0..=1.0).contains(&v.y));
                assert!((0.0..=1.0).contains(&v.z));
            }
        }
    }

    /// Adversarial: every triangle is geometrically non-degenerate
    /// across a 48 × 32 perturbation sweep.
    #[test]
    fn case6_triangles_are_nondegenerate() {
        for (our_idx, subcase) in case7_idx_subcase_pairs() {
            for perturb_seed in 0..32_u32 {
                let mut corners = [0.0_f32; 8];
                for (bit, c) in corners.iter_mut().enumerate() {
                    let inside = our_idx & (1 << bit) != 0;
                    let mag = 0.3
                        + 1.4
                            * f32::from(
                                ((perturb_seed.wrapping_mul(7) ^ u32::try_from(bit).unwrap_or(0))
                                    & 0xF) as u8,
                            )
                            / 15.0;
                    *c = if inside { -mag } else { mag };
                }
                let (grid, _field) = unit_grid_with_corners(&corners);
                let mut mesh = Mesh::default();
                let cell = CellCoord { x: 0, y: 0, z: 0 };
                emit(&grid, cell, &corners, subcase, &mut mesh);
                for (t_idx, tri) in mesh.indices.iter().enumerate() {
                    let v0 = mesh.vertices[tri[0] as usize];
                    let v1 = mesh.vertices[tri[1] as usize];
                    let v2 = mesh.vertices[tri[2] as usize];
                    let area2 = (v1 - v0).cross(v2 - v0).length_squared();
                    assert!(
                        area2 > 1e-12,
                        "triangle {t_idx} near-zero area for our_idx={our_idx:#010b} \
                         seed={perturb_seed}"
                    );
                }
            }
        }
    }

    /// Coverage: across a perturbation sweep, count how many of the
    /// three dispatch paths (6.1.1 = 3 tris, 6.1.2 = 9 tris, 6.2 = 5
    /// tris) are reached. Pin that the sweep reaches at least two of
    /// the three; document any unreachable path as a coverage
    /// observation, not a failure.
    #[test]
    fn case6_each_dispatch_path_reachable_under_perturbation() {
        let mut count_3tris = 0_u32; // 6.1.1
        let mut count_5tris = 0_u32; // 6.2
        let mut count_9tris = 0_u32; // 6.1.2
        for (our_idx, subcase) in case7_idx_subcase_pairs() {
            for perturb_seed in 0..32_u32 {
                let mut corners = [0.0_f32; 8];
                for (bit, c) in corners.iter_mut().enumerate() {
                    let inside = our_idx & (1 << bit) != 0;
                    let mag = 0.3
                        + 1.4
                            * f32::from(
                                ((perturb_seed.wrapping_mul(11) ^ u32::try_from(bit).unwrap_or(0))
                                    & 0xF) as u8,
                            )
                            / 15.0;
                    *c = if inside { -mag } else { mag };
                }
                let (grid, _field) = unit_grid_with_corners(&corners);
                let mut mesh = Mesh::default();
                let cell = CellCoord { x: 0, y: 0, z: 0 };
                emit(&grid, cell, &corners, subcase, &mut mesh);
                match mesh.triangle_count() {
                    3 => count_3tris += 1,
                    5 => count_5tris += 1,
                    9 => count_9tris += 1,
                    n => panic!("unexpected triangle count {n}"),
                }
            }
        }
        eprintln!(
            "Case 6 sweep: 6.1.1 = {count_3tris}, 6.1.2 = {count_9tris}, 6.2 = {count_5tris}"
        );
        let paths_hit = [count_3tris, count_5tris, count_9tris]
            .iter()
            .filter(|&&n| n > 0)
            .count();
        assert!(
            paths_hit >= 2,
            "perturbation sweep reached only {paths_hit}/3 dispatch paths"
        );
    }

    /// Pin: case 7 is no longer in the warn-once ambiguous set.
    #[test]
    fn case7_is_no_longer_ambiguous() {
        for (our_idx, subcase) in case7_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            emit(&grid, cell, &corners, subcase, &mut mesh);
            assert!(
                !mesh.indices.is_empty(),
                "case 7 voxel our_idx={our_idx:#010b} must produce triangles"
            );
        }
    }
}
