//! Case 12 (Chernyaev) / Lewiner base case 13 — "3+1 + diagonal".
//!
//! The configuration: three same-sign corners adjacent on one face plus
//! a fourth same-sign corner diagonally opposite. There are 24 voxel
//! configurations (`base_case == 13` rows in [`super::cases::CASES`]).
//!
//! # Dispatch
//!
//! Identical structure to [`super::case10`]: 2 face tests + 1 interior
//! test, 5 sub-tilings. **Different interior decider**: case 12 uses
//! [`interior_decider_per_edge`] (the per-edge formula, same primitive
//! as Chernyaev 6/7/13) rather than the body-diagonal formula.
//!
//! | (face1, face2) | Subcase | Triangles | Needs c-vertex? |
//! |----------------|---------|-----------|-----------------|
//! | (true, true)   | 12.1.1_ | 4         | no              |
//! | (true, false)  | 12.2    | 8         | yes             |
//! | (false, true)  | 12.2_   | 8         | yes             |
//! | (false, false) | interior decides 12.1.1 vs 12.1.2 | 4 or 8 | 12.1.2 only |
//!
//! # Sign / index conventions
//!
//! - **`TEST_12` columns**:
//!   - 0/1: face indices for [`face_decider`].
//!   - 2: per-edge interior sentinel `s` (always `7`).
//!   - 3: reference edge for [`interior_decider_per_edge`] (already 0-based, range 0..=11).
//!
//! # Port
//!
//! Tables transcribed from `JuliaGeometry/MarchingCubes.jl` `src/lut.jl`:
//! - `test12` at lines 1346-1372
//! - `tiling12_1_1` at lines 1383-1418
//! - `tiling12_1_1_` at lines 1420-1455
//! - `tiling12_1_2` at lines 1457-1492
//! - `tiling12_2` at lines 1494-1529
//! - `tiling12_2_` at lines 1531-1566
//!
//! Dispatch mirrors `case == 13` in `src/MarchingCubes.jl`, lines 279-296.

use super::super::{CellCoord, Mesh};
use super::decider::{face_decider, interior_decider_per_edge};
use super::emit_triangles;
use crate::grid::GridSpec;

/// Per-cfg 4-tuple: `[face_a, face_b, interior_s, ref_edge]`.
///
/// Source: Lewiner `test12[24][4]`. Col 2 is the per-edge sentinel
/// (always `7`). Col 3 is the reference edge — **already 0-based** in
/// the source.
#[rustfmt::skip]
const TEST_12: [[i8; 4]; 24] = [
    [  4,   3,   7,  11],
    [  3,   2,   7,  10],
    [  2,   6,   7,   5],
    [  6,   4,   7,   7],
    [  2,   1,   7,   9],
    [  5,   2,   7,   1],
    [  5,   3,   7,   2],
    [  5,   1,   7,   0],
    [  5,   4,   7,   3],
    [  6,   3,   7,   6],
    [  1,   6,   7,   4],
    [  1,   4,   7,   8],
    [  4,   1,   7,   8],
    [  6,   1,   7,   4],
    [  3,   6,   7,   6],
    [  4,   5,   7,   3],
    [  1,   5,   7,   0],
    [  3,   5,   7,   2],
    [  2,   5,   7,   1],
    [  1,   2,   7,   9],
    [  4,   6,   7,   7],
    [  6,   2,   7,   5],
    [  2,   3,   7,  10],
    [  3,   4,   7,  11],
];

/// Triangulation for subcase 12.1.1 (4 triangles, no c-vertex).
///
/// Source: Lewiner `tiling12_1_1[24]`.
#[rustfmt::skip]
const TILING_12_1_1: [[u8; 12]; 24] = [
    [  7,   6,  11,  10,   3,   2,   3,  10,   8,   9,   8,  10],
    [  6,   5,  10,   9,   2,   1,   2,   9,  11,   8,  11,   9],
    [ 10,   6,   5,   7,   9,   4,   9,   7,   1,   3,   1,   7],
    [  7,   6,  11,   4,   8,   5,   3,   5,   8,   5,   3,   1],
    [  5,   4,   9,   8,   1,   0,   1,   8,  10,  11,  10,   8],
    [  1,   2,  10,   0,   9,   3,   5,   3,   9,   3,   5,   7],
    [ 10,   1,   2,   0,  11,   3,  11,   0,   6,   4,   6,   0],
    [  8,   3,   0,   2,   9,   1,   9,   2,   4,   6,   4,   2],
    [  3,   0,   8,   2,  11,   1,   7,   1,  11,   1,   7,   5],
    [  6,   5,  10,   7,  11,   4,   2,   4,  11,   4,   2,   0],
    [  9,   5,   4,   6,   8,   7,   8,   6,   0,   2,   0,   6],
    [  8,   3,   0,   7,   4,  11,   9,  11,   4,  11,   9,  10],
    [  4,   7,   8,  11,   0,   3,   0,  11,   9,  10,   9,  11],
    [  4,   7,   8,   5,   9,   6,   0,   6,   9,   6,   0,   2],
    [ 11,   7,   6,   4,  10,   5,  10,   4,   2,   0,   2,   4],
    [ 11,   2,   3,   1,   8,   0,   8,   1,   7,   5,   7,   1],
    [  0,   1,   9,   3,   8,   2,   4,   2,   8,   2,   4,   6],
    [  2,   3,  11,   1,  10,   0,   6,   0,  10,   0,   6,   4],
    [  9,   0,   1,   3,  10,   2,  10,   3,   5,   7,   5,   3],
    [  9,   0,   1,   4,   5,   8,  10,   8,   5,   8,  10,  11],
    [  8,   4,   7,   5,  11,   6,  11,   5,   3,   1,   3,   5],
    [  5,   4,   9,   6,  10,   7,   1,   7,  10,   7,   1,   3],
    [ 10,   1,   2,   5,   6,   9,  11,   9,   6,   9,  11,   8],
    [ 11,   2,   3,   6,   7,  10,   8,  10,   7,  10,   8,   9],
];

/// Triangulation for subcase 12.1.1 inverted polarity (4 triangles, no c-vertex).
///
/// Source: Lewiner `tiling12_1_1_[24]`.
#[rustfmt::skip]
const TILING_12_1_1_: [[u8; 12]; 24] = [
    [  3,   2,  11,  10,   7,   6,   7,  10,   8,   9,   8,  10],
    [  2,   1,  10,   9,   6,   5,   6,   9,  11,   8,  11,   9],
    [  9,   4,   5,   7,  10,   6,  10,   7,   1,   3,   1,   7],
    [  7,   4,   8,   6,  11,   5,   3,   5,  11,   5,   3,   1],
    [  1,   0,   9,   8,   5,   4,   5,   8,  10,  11,  10,   8],
    [  1,   0,   9,   2,  10,   3,   5,   3,  10,   3,   5,   7],
    [ 11,   3,   2,   0,  10,   1,  10,   0,   6,   4,   6,   0],
    [  9,   1,   0,   2,   8,   3,   8,   2,   4,   6,   4,   2],
    [  3,   2,  11,   0,   8,   1,   7,   1,   8,   1,   7,   5],
    [  6,   7,  11,   5,  10,   4,   2,   4,  10,   4,   2,   0],
    [  8,   7,   4,   6,   9,   5,   9,   6,   0,   2,   0,   6],
    [  8,   7,   4,   3,   0,  11,   9,  11,   0,  11,   9,  10],
    [  0,   3,   8,  11,   4,   7,   4,  11,   9,  10,   9,  11],
    [  4,   5,   9,   7,   8,   6,   0,   6,   8,   6,   0,   2],
    [ 10,   5,   6,   4,  11,   7,  11,   4,   2,   0,   2,   4],
    [  8,   0,   3,   1,  11,   2,  11,   1,   7,   5,   7,   1],
    [  0,   3,   8,   1,   9,   2,   4,   2,   9,   2,   4,   6],
    [  2,   1,  10,   3,  11,   0,   6,   0,  11,   0,   6,   4],
    [ 10,   2,   1,   3,   9,   0,   9,   3,   5,   7,   5,   3],
    [  9,   4,   5,   0,   1,   8,  10,   8,   1,   8,  10,  11],
    [ 11,   6,   7,   5,   8,   4,   8,   5,   3,   1,   3,   5],
    [  5,   6,  10,   4,   9,   7,   1,   7,   9,   7,   1,   3],
    [ 10,   5,   6,   1,   2,   9,  11,   9,   2,   9,  11,   8],
    [ 11,   6,   7,   2,   3,  10,   8,  10,   3,  10,   8,   9],
];

/// Triangulation for subcase 12.1.2 (8 triangles, no c-vertex).
///
/// Source: Lewiner `tiling12_1_2[24]`.
#[rustfmt::skip]
const TILING_12_1_2: [[u8; 24]; 24] = [
    [  7,   3,  11,   3,   7,   8,   9,   8,   7,   6,   9,   7,   9,   6,  10,   2,  10,   6,  11,   2,   6,   2,  11,   3],
    [  6,   2,  10,   2,   6,  11,   8,  11,   6,   5,   8,   6,   8,   5,   9,   1,   9,   5,  10,   1,   5,   1,  10,   2],
    [ 10,   9,   5,   9,  10,   1,   3,   1,  10,   6,   3,  10,   3,   6,   7,   4,   7,   6,   5,   4,   6,   4,   5,   9],
    [  7,   8,  11,   3,  11,   8,  11,   3,   1,  11,   1,   6,   5,   6,   1,   6,   5,   4,   6,   4,   7,   8,   7,   4],
    [  5,   1,   9,   1,   5,  10,  11,  10,   5,   4,  11,   5,  11,   4,   8,   0,   8,   4,   9,   0,   4,   0,   9,   1],
    [  1,   9,  10,   5,  10,   9,  10,   5,   7,  10,   7,   2,   3,   2,   7,   2,   3,   0,   2,   0,   1,   9,   1,   0],
    [ 10,  11,   2,  11,  10,   6,   4,   6,  10,   1,   4,  10,   4,   1,   0,   3,   0,   1,   2,   3,   1,   3,   2,  11],
    [  8,   9,   0,   9,   8,   4,   6,   4,   8,   3,   6,   8,   6,   3,   2,   1,   2,   3,   0,   1,   3,   1,   0,   9],
    [  3,  11,   8,   7,   8,  11,   8,   7,   5,   8,   5,   0,   1,   0,   5,   0,   1,   2,   0,   2,   3,  11,   3,   2],
    [  6,  11,  10,   2,  10,  11,  10,   2,   0,  10,   0,   5,   4,   5,   0,   5,   4,   7,   5,   7,   6,  11,   6,   7],
    [  9,   8,   4,   8,   9,   0,   2,   0,   9,   5,   2,   9,   2,   5,   6,   7,   6,   5,   4,   7,   5,   7,   4,   8],
    [  8,   4,   0,   9,   0,   4,   0,   9,  10,   0,  10,   3,  11,   3,  10,   3,  11,   7,   3,   7,   8,   4,   8,   7],
    [  4,   0,   8,   0,   4,   9,  10,   9,   4,   7,  10,   4,  10,   7,  11,   3,  11,   7,   8,   3,   7,   3,   8,   0],
    [  4,   9,   8,   0,   8,   9,   8,   0,   2,   8,   2,   7,   6,   7,   2,   7,   6,   5,   7,   5,   4,   9,   4,   5],
    [ 11,  10,   6,  10,  11,   2,   0,   2,  11,   7,   0,  11,   0,   7,   4,   5,   4,   7,   6,   5,   7,   5,   6,  10],
    [ 11,   8,   3,   8,  11,   7,   5,   7,  11,   2,   5,  11,   5,   2,   1,   0,   1,   2,   3,   0,   2,   0,   3,   8],
    [  0,   8,   9,   4,   9,   8,   9,   4,   6,   9,   6,   1,   2,   1,   6,   1,   2,   3,   1,   3,   0,   8,   0,   3],
    [  2,  10,  11,   6,  11,  10,  11,   6,   4,  11,   4,   3,   0,   3,   4,   3,   0,   1,   3,   1,   2,  10,   2,   1],
    [  9,  10,   1,  10,   9,   5,   7,   5,   9,   0,   7,   9,   7,   0,   3,   2,   3,   0,   1,   2,   0,   2,   1,  10],
    [  9,   5,   1,  10,   1,   5,   1,  10,  11,   1,  11,   0,   8,   0,  11,   0,   8,   4,   0,   4,   9,   5,   9,   4],
    [  8,  11,   7,  11,   8,   3,   1,   3,   8,   4,   1,   8,   1,   4,   5,   6,   5,   4,   7,   6,   4,   6,   7,  11],
    [  5,  10,   9,   1,   9,  10,   9,   1,   3,   9,   3,   4,   7,   4,   3,   4,   7,   6,   4,   6,   5,  10,   5,   6],
    [ 10,   6,   2,  11,   2,   6,   2,  11,   8,   2,   8,   1,   9,   1,   8,   1,   9,   5,   1,   5,  10,   6,  10,   5],
    [ 11,   7,   3,   8,   3,   7,   3,   8,   9,   3,   9,   2,  10,   2,   9,   2,  10,   6,   2,   6,  11,   7,  11,   6],
];

/// Triangulation for subcase 12.2 (8 triangles, uses c-vertex slot 12).
///
/// Source: Lewiner `tiling12_2[24]`.
#[rustfmt::skip]
const TILING_12_2: [[u8; 24]; 24] = [
    [  9,   8,  12,  10,   9,  12,   2,  10,  12,   3,   2,  12,  11,   3,  12,   6,  11,  12,   7,   6,  12,   8,   7,  12],
    [  8,  11,  12,   9,   8,  12,   1,   9,  12,   2,   1,  12,  10,   2,  12,   5,  10,  12,   6,   5,  12,  11,   6,  12],
    [  3,   1,  12,   7,   3,  12,   4,   7,  12,   9,   4,  12,   5,   9,  12,   6,   5,  12,  10,   6,  12,   1,  10,  12],
    [ 12,   3,   1,  12,   1,   5,  12,   5,   6,  12,   6,  11,  12,  11,   7,  12,   7,   4,  12,   4,   8,  12,   8,   3],
    [ 11,  10,  12,   8,  11,  12,   0,   8,  12,   1,   0,  12,   9,   1,  12,   4,   9,  12,   5,   4,  12,  10,   5,  12],
    [ 12,   5,   7,  12,   7,   3,  12,   3,   2,  12,   2,  10,  12,  10,   1,  12,   1,   0,  12,   0,   9,  12,   9,   5],
    [  4,   6,  12,   0,   4,  12,   1,   0,  12,  10,   1,  12,   2,  10,  12,   3,   2,  12,  11,   3,  12,   6,  11,  12],
    [  6,   4,  12,   2,   6,  12,   3,   2,  12,   8,   3,  12,   0,   8,  12,   1,   0,  12,   9,   1,  12,   4,   9,  12],
    [ 12,   7,   5,  12,   5,   1,  12,   1,   0,  12,   0,   8,  12,   8,   3,  12,   3,   2,  12,   2,  11,  12,  11,   7],
    [ 12,   2,   0,  12,   0,   4,  12,   4,   5,  12,   5,  10,  12,  10,   6,  12,   6,   7,  12,   7,  11,  12,  11,   2],
    [  2,   0,  12,   6,   2,  12,   7,   6,  12,   8,   7,  12,   4,   8,  12,   5,   4,  12,   9,   5,  12,   0,   9,  12],
    [ 12,   9,  10,  12,  10,  11,  12,  11,   7,  12,   7,   4,  12,   4,   8,  12,   8,   3,  12,   3,   0,  12,   0,   9],
    [ 10,   9,  12,  11,  10,  12,   7,  11,  12,   4,   7,  12,   8,   4,  12,   3,   8,  12,   0,   3,  12,   9,   0,  12],
    [ 12,   0,   2,  12,   2,   6,  12,   6,   7,  12,   7,   8,  12,   8,   4,  12,   4,   5,  12,   5,   9,  12,   9,   0],
    [  0,   2,  12,   4,   0,  12,   5,   4,  12,  10,   5,  12,   6,  10,  12,   7,   6,  12,  11,   7,  12,   2,  11,  12],
    [  5,   7,  12,   1,   5,  12,   0,   1,  12,   8,   0,  12,   3,   8,  12,   2,   3,  12,  11,   2,  12,   7,  11,  12],
    [ 12,   4,   6,  12,   6,   2,  12,   2,   3,  12,   3,   8,  12,   8,   0,  12,   0,   1,  12,   1,   9,  12,   9,   4],
    [ 12,   6,   4,  12,   4,   0,  12,   0,   1,  12,   1,  10,  12,  10,   2,  12,   2,   3,  12,   3,  11,  12,  11,   6],
    [  7,   5,  12,   3,   7,  12,   2,   3,  12,  10,   2,  12,   1,  10,  12,   0,   1,  12,   9,   0,  12,   5,   9,  12],
    [ 12,  10,  11,  12,  11,   8,  12,   8,   0,  12,   0,   1,  12,   1,   9,  12,   9,   4,  12,   4,   5,  12,   5,  10],
    [  1,   3,  12,   5,   1,  12,   6,   5,  12,  11,   6,  12,   7,  11,  12,   4,   7,  12,   8,   4,  12,   3,   8,  12],
    [ 12,   1,   3,  12,   3,   7,  12,   7,   4,  12,   4,   9,  12,   9,   5,  12,   5,   6,  12,   6,  10,  12,  10,   1],
    [ 12,  11,   8,  12,   8,   9,  12,   9,   1,  12,   1,   2,  12,   2,  10,  12,  10,   5,  12,   5,   6,  12,   6,  11],
    [ 12,   8,   9,  12,   9,  10,  12,  10,   2,  12,   2,   3,  12,   3,  11,  12,  11,   6,  12,   6,   7,  12,   7,   8],
];

/// Triangulation for subcase 12.2 inverted polarity (8 triangles, c-vertex).
///
/// Source: Lewiner `tiling12_2_[24]`.
#[rustfmt::skip]
const TILING_12_2_: [[u8; 24]; 24] = [
    [ 12,   2,  11,  12,  11,   7,  12,   7,   6,  12,   6,  10,  12,  10,   9,  12,   9,   8,  12,   8,   3,  12,   3,   2],
    [ 12,   1,  10,  12,  10,   6,  12,   6,   5,  12,   5,   9,  12,   9,   8,  12,   8,  11,  12,  11,   2,  12,   2,   1],
    [ 12,   4,   5,  12,   5,  10,  12,  10,   6,  12,   6,   7,  12,   7,   3,  12,   3,   1,  12,   1,   9,  12,   9,   4],
    [  7,   6,  12,   8,   7,  12,   4,   8,  12,   5,   4,  12,   1,   5,  12,   3,   1,  12,  11,   3,  12,   6,  11,  12],
    [ 12,   0,   9,  12,   9,   5,  12,   5,   4,  12,   4,   8,  12,   8,  11,  12,  11,  10,  12,  10,   1,  12,   1,   0],
    [  1,   2,  12,   9,   1,  12,   0,   9,  12,   3,   0,  12,   7,   3,  12,   5,   7,  12,  10,   5,  12,   2,  10,  12],
    [ 12,   1,   2,  12,   2,  11,  12,  11,   3,  12,   3,   0,  12,   0,   4,  12,   4,   6,  12,   6,  10,  12,  10,   1],
    [ 12,   3,   0,  12,   0,   9,  12,   9,   1,  12,   1,   2,  12,   2,   6,  12,   6,   4,  12,   4,   8,  12,   8,   3],
    [  3,   0,  12,  11,   3,  12,   2,  11,  12,   1,   2,  12,   5,   1,  12,   7,   5,  12,   8,   7,  12,   0,   8,  12],
    [  6,   5,  12,  11,   6,  12,   7,  11,  12,   4,   7,  12,   0,   4,  12,   2,   0,  12,  10,   2,  12,   5,  10,  12],
    [ 12,   7,   4,  12,   4,   9,  12,   9,   5,  12,   5,   6,  12,   6,   2,  12,   2,   0,  12,   0,   8,  12,   8,   7],
    [  8,   7,  12,   0,   8,  12,   3,   0,  12,  11,   3,  12,  10,  11,  12,   9,  10,  12,   4,   9,  12,   7,   4,  12],
    [ 12,   7,   8,  12,   8,   0,  12,   0,   3,  12,   3,  11,  12,  11,  10,  12,  10,   9,  12,   9,   4,  12,   4,   7],
    [  4,   7,  12,   9,   4,  12,   5,   9,  12,   6,   5,  12,   2,   6,  12,   0,   2,  12,   8,   0,  12,   7,   8,  12],
    [ 12,   5,   6,  12,   6,  11,  12,  11,   7,  12,   7,   4,  12,   4,   0,  12,   0,   2,  12,   2,  10,  12,  10,   5],
    [ 12,   0,   3,  12,   3,  11,  12,  11,   2,  12,   2,   1,  12,   1,   5,  12,   5,   7,  12,   7,   8,  12,   8,   0],
    [  0,   3,  12,   9,   0,  12,   1,   9,  12,   2,   1,  12,   6,   2,  12,   4,   6,  12,   8,   4,  12,   3,   8,  12],
    [  2,   1,  12,  11,   2,  12,   3,  11,  12,   0,   3,  12,   4,   0,  12,   6,   4,  12,  10,   6,  12,   1,  10,  12],
    [ 12,   2,   1,  12,   1,   9,  12,   9,   0,  12,   0,   3,  12,   3,   7,  12,   7,   5,  12,   5,  10,  12,  10,   2],
    [  9,   0,  12,   5,   9,  12,   4,   5,  12,   8,   4,  12,  11,   8,  12,  10,  11,  12,   1,  10,  12,   0,   1,  12],
    [ 12,   6,   7,  12,   7,   8,  12,   8,   4,  12,   4,   5,  12,   5,   1,  12,   1,   3,  12,   3,  11,  12,  11,   6],
    [  5,   4,  12,  10,   5,  12,   6,  10,  12,   7,   6,  12,   3,   7,  12,   1,   3,  12,   9,   1,  12,   4,   9,  12],
    [ 10,   1,  12,   6,  10,  12,   5,   6,  12,   9,   5,  12,   8,   9,  12,  11,   8,  12,   2,  11,  12,   1,   2,  12],
    [ 11,   2,  12,   7,  11,  12,   6,   7,  12,  10,   6,  12,   9,  10,  12,   8,   9,  12,   3,   8,  12,   2,   3,  12],
];

/// Emits Case-12 triangles for the given voxel.
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(
        (1..=24).contains(&subcase),
        "Case 12 subcase out of range: {subcase}"
    );
    let cfg = usize::try_from(subcase - 1).unwrap_or(0);
    let f1 = face_decider(corners, TEST_12[cfg][0]);
    let f2 = face_decider(corners, TEST_12[cfg][1]);
    let tiling: &[u8] = match (f1, f2) {
        (true, true) => &TILING_12_1_1_[cfg],
        (true, false) => &TILING_12_2[cfg],
        (false, true) => &TILING_12_2_[cfg],
        (false, false) => {
            let s = TEST_12[cfg][2];
            let edge = u8::try_from(TEST_12[cfg][3]).unwrap_or(0);
            if interior_decider_per_edge(corners, edge, s) {
                &TILING_12_1_1[cfg]
            } else {
                &TILING_12_1_2[cfg]
            }
        }
    };
    emit_triangles(grid, cell, corners, tiling, mesh);
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
    use crate::Mesh;
    use crate::grid::GridSpec;
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

    fn case13_idx_subcase_pairs() -> Vec<(u8, i8)> {
        let mut hits = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 13 {
                hits.push((our_idx, subcase));
            }
        }
        hits
    }

    #[test]
    fn case13_dispatch_has_twenty_four_configs() {
        assert_eq!(case13_idx_subcase_pairs().len(), 24);
    }

    /// Pin: `TEST_12[cfg][2] == 7` and `TEST_12[cfg][3] ∈ 0..=11` for all cfg.
    #[test]
    fn test_12_columns_2_3_invariants() {
        for (cfg, row) in TEST_12.iter().enumerate() {
            assert_eq!(row[2], 7, "TEST_12[{cfg}][2] = {} (expected 7)", row[2]);
            assert!(
                (0..=11).contains(&row[3]),
                "TEST_12[{cfg}][3] = {} not 0..=11",
                row[3]
            );
        }
    }

    /// Transcription oracle: `test12[1]` (Julia line 1347 → `(4, 3, 7, 11)`).
    #[test]
    fn test_12_first_row_matches_lut() {
        assert_eq!(TEST_12[0], [4, 3, 7, 11]);
    }

    /// Transcription oracle: `tiling12_1_1[1]` (Julia line 1384 → first row).
    /// The published row is `(7, 6, 11, 10, 3, 2, 3, 10, 8, 9, 8, 10)`
    /// → 0-based.
    #[test]
    fn tiling_12_1_1_first_row_matches_lut() {
        assert_eq!(TILING_12_1_1[0], [7, 6, 11, 10, 3, 2, 3, 10, 8, 9, 8, 10]);
    }

    #[test]
    fn tilings_use_only_legal_indices() {
        for cfg in 0..24 {
            for &e in &TILING_12_1_1[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_12_1_1_[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_12_1_2[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_12_2[cfg] {
                assert!(e <= 12);
            }
            for &e in &TILING_12_2_[cfg] {
                assert!(e <= 12);
            }
        }
    }

    #[test]
    fn tiling_dimensions_match_lewiner() {
        assert_eq!(TEST_12.len(), 24);
        assert_eq!(TILING_12_1_1.len(), 24);
        assert_eq!(TILING_12_1_1[0].len(), 12);
        assert_eq!(TILING_12_1_2[0].len(), 24);
        assert_eq!(TILING_12_2[0].len(), 24);
    }

    #[test]
    fn case12_dispatch_runs_on_every_case13_voxel() {
        for (our_idx, subcase) in case13_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            emit(
                &grid,
                CellCoord { x: 0, y: 0, z: 0 },
                &corners,
                subcase,
                &mut mesh,
            );
            let n = mesh.triangle_count();
            assert!(matches!(n, 4 | 8), "expected 4/8, got {n}");
            assert_eq!(mesh.vertices.len(), 3 * n);
            for v in &mesh.vertices {
                assert!(v.is_finite());
                assert!((0.0..=1.0).contains(&v.x));
                assert!((0.0..=1.0).contains(&v.y));
                assert!((0.0..=1.0).contains(&v.z));
            }
        }
    }

    #[test]
    fn case12_triangles_are_nondegenerate() {
        for (our_idx, subcase) in case13_idx_subcase_pairs() {
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
                let (grid, _) = unit_grid_with_corners(&corners);
                let mut mesh = Mesh::default();
                emit(
                    &grid,
                    CellCoord { x: 0, y: 0, z: 0 },
                    &corners,
                    subcase,
                    &mut mesh,
                );
                for tri in &mesh.indices {
                    let v0 = mesh.vertices[tri[0] as usize];
                    let v1 = mesh.vertices[tri[1] as usize];
                    let v2 = mesh.vertices[tri[2] as usize];
                    let area2 = (v1 - v0).cross(v2 - v0).length_squared();
                    assert!(area2 > 1e-12);
                }
            }
        }
    }
}
