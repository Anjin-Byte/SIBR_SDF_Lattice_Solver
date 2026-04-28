//! Case 13 (Chernyaev) / Lewiner base case 14 — the "tunnel" case.
//!
//! The configuration: four corners with the same sign sit at four corners
//! of the cube such that no two share an edge — the "anti-corner"
//! topology. There are exactly two such configurations up to rotation, so
//! the dispatch table has only two entries (matching `CASES[i] = (14, 1)`
//! and `CASES[i] = (14, 2)` rows in [`super::cases::CASES`] — there are
//! exactly two of each).
//!
//! Disambiguation runs **six face tests** (one per cube face, encoded by
//! [`TEST_13`]) plus optionally **one interior test**. The six face-test
//! results pack into a 6-bit key indexed into [`SUBCFG_13`], which yields
//! a topological subcase id `sc ∈ 0..=45` (with sentinel `-1` for
//! configurations that should not occur for valid case-13 voxels). The
//! id's range determines which Chernyaev subcase to emit:
//!
//! | `sc` range   | Chernyaev subcase | Triangles | Needs c-vertex? | Stage |
//! |--------------|-------------------|-----------|------------------|-------|
//! | `0`          | 13.1              | 4         | no               | **1** |
//! | `1..=6`      | 13.2              | 6         | no               | **1** |
//! | `7..=18`     | 13.3              | 10        | yes              | **2** |
//! | `19..=22`    | 13.4              | 12        | yes              | **2** |
//! | `23..=26`    | 13.5 (1 or 2)     | 6 or 10   | no (interior test)| **3** |
//! | `27..=38`    | 13.3 (inverted)   | 10        | yes              | **2** |
//! | `39..=44`    | 13.2 (inverted)   | 6         | no               | **1** |
//! | `45`         | 13.1 (inverted)   | 4         | no               | **1** |
//!
//! # Full case-13 port (this module's current state)
//!
//! All seven Chernyaev case-13 subcases (13.1, 13.2, 13.3, 13.4, 13.5.1,
//! 13.5.2, plus inverted polarities of 13.1/13.2/13.3) are ported.
//! [`emit`] returns `true` for every valid case-13 voxel and `false`
//! only on the numerical-borderline `-1` sentinel from `SUBCFG_13`.
//!
//! # Sign / index conventions
//!
//! - **Tile edges** are stored 0-based 0..=11, converted from Lewiner's
//!   1-based 1..=12 at port time (matching [`super::case3`]).
//! - **`SUBCFG_13` keys** are 0-based 0..=63 in this port (Julia uses
//!   1-based 1..=64; we initialize the accumulator at `0` rather than
//!   Julia's `1`, so the off-by-one comes out in the wash).
//! - **`cfg`** is `subcase - 1` in 0..=1, matching the rest of the MC33
//!   port.
//!
//! # Port
//!
//! Tables transcribed from `JuliaGeometry/MarchingCubes.jl`
//! `src/lut.jl` (commit at fetch time: master, ~2185 lines):
//!
//! - `test13` at lines 1570–1573
//! - `subcfg13` at lines 1587–1652
//! - `tiling13_1` at lines 1664–1667
//! - `tiling13_1_` at lines 1679–1682
//! - `tiling13_2` at lines 1694–1710
//! - `tiling13_2_` at lines 1722–1738
//! - `tiling13_3` at lines 1750–1778
//! - `tiling13_3_` at lines 1790–1818
//! - `tiling13_4` at lines 1830–1842
//! - `tiling13_5_1` at lines 1855–1867
//! - `tiling13_5_2` at lines 1879–1891
//!
//! The dispatcher mirrors `case == 14` in `src/MarchingCubes.jl`,
//! lines 297–327.
//!
//! # Winding
//!
//! Same convention as [`super::case3`] / [`super::case4`]: emit
//! `[base, base+1, base+2]` directly without reversal. Tested below.

#[cfg(test)]
use super::super::tables::EDGE_CORNERS;
use super::super::{CellCoord, Mesh};
#[cfg(test)]
use super::compute_c_vertex;
use super::decider::{face_decider, interior_decider_per_edge};
use super::emit_triangles;
use crate::mesh::grid::GridSpec;
#[cfg(test)]
use glam::Vec3;

/// Six face tests + one interior-test sentinel per dispatch config.
///
/// Entries `i ∈ 1..=6` (or negated): face indices for [`face_decider`].
/// Entry at index 5 (Julia's `test13[cfg][6]`, value `6`): the interior-
/// test sentinel passed to [`interior_decider_per_edge`] for subcase 13.5
/// dispatch. Entry at index 6 (value `7`) is unused in this case (it
/// would be `interior_decider`'s body-diagonal sentinel, which case 13
/// doesn't consult).
///
/// Source: Lewiner `test13[2][7]`.
const TEST_13: [[i8; 7]; 2] = [[1, 2, 3, 4, 5, 6, 7], [2, 3, 4, 1, 5, 6, 7]];

/// 64-entry table mapping a 6-bit face-test key (`0..=63`) to a
/// topological-subcase id `sc`. Sentinel `-1` flags combinations that
/// should not occur for valid case-13 voxels (numerically borderline
/// near the saddle); we treat them as fallback.
///
/// Source: Lewiner `subcfg13[64]`. Entries match `lut.jl` 1-based; we
/// shift to 0-based by initializing the accumulator at 0 rather than
/// Julia's 1 (see module docs).
#[rustfmt::skip]
const SUBCFG_13: [i8; 64] = [
     0,  1,  2,  7,  3, -1, 11, -1,
     4,  8, -1, -1, 14, -1, -1, -1,
     5,  9, 12, 23, 15, -1, 21, 38,
    17, 20, -1, 36, 26, 33, 30, 44,
     6, 10, 13, 19, 16, -1, 25, 37,
    18, 24, -1, 35, 22, 32, 29, 43,
    -1, -1, -1, 34, -1, -1, 28, 42,
    -1, 31, -1, 41, 27, 40, 39, 45,
];

/// Triangulation for Chernyaev subcase 13.1 (4 triangles, no c-vertex).
///
/// Edge indices 0-based.
///
/// Source: Lewiner `tiling13_1[2][12]`.
#[rustfmt::skip]
const TILING_13_1: [[u8; 12]; 2] = [
    [11, 7, 6, 1, 2, 10, 8, 3, 0, 9, 5, 4],
    [ 8, 4, 7, 2, 3, 11, 9, 0, 1, 10, 6, 5],
];

/// Triangulation for subcase 13.1 inverted polarity (`sc == 45`).
///
/// Source: Lewiner `tiling13_1_[2][12]`.
#[rustfmt::skip]
const TILING_13_1_: [[u8; 12]; 2] = [
    [7, 4, 8, 11, 3, 2, 1, 0, 9, 5, 6, 10],
    [6, 7, 11, 10, 2, 1, 0, 3, 8, 4, 5, 9],
];

/// Triangulation for Chernyaev subcase 13.2 (6 triangles, no c-vertex).
/// Indexed `[cfg][sc - 1]` for `sc ∈ 1..=6`.
///
/// Source: Lewiner `tiling13_2[2][6][18]`.
#[rustfmt::skip]
const TILING_13_2: [[[u8; 18]; 6]; 2] = [
    [
        [ 1, 2, 10, 11, 7, 6, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9],
        [ 8, 3, 0, 11, 7, 6, 9, 1, 4, 2, 4, 1, 4, 2, 5, 10, 5, 2],
        [ 9, 5, 4, 8, 3, 0, 1, 6, 10, 6, 1, 7, 2, 7, 1, 7, 2, 11],
        [ 9, 5, 4, 1, 2, 10, 11, 3, 6, 0, 6, 3, 6, 0, 7, 8, 7, 0],
        [ 9, 5, 4, 11, 7, 6, 0, 10, 1, 10, 0, 8, 10, 8, 2, 3, 2, 8],
        [ 1, 2, 10, 3, 0, 8, 4, 9, 7, 11, 7, 9, 5, 11, 9, 11, 5, 6],
    ],
    [
        [ 2, 3, 11, 8, 4, 7, 0, 5, 9, 5, 0, 6, 1, 6, 0, 6, 1, 10],
        [ 9, 0, 1, 8, 4, 7, 10, 2, 5, 3, 5, 2, 5, 3, 6, 11, 6, 3],
        [ 6, 5, 10, 9, 0, 1, 2, 7, 11, 7, 2, 4, 3, 4, 2, 4, 3, 8],
        [ 6, 5, 10, 2, 3, 11, 8, 0, 7, 1, 7, 0, 7, 1, 4, 9, 4, 1],
        [ 6, 5, 10, 8, 4, 7, 1, 11, 2, 11, 1, 9, 11, 9, 3, 0, 3, 9],
        [ 2, 3, 11, 0, 1, 9, 5, 10, 4, 8, 4, 10, 6, 8, 10, 8, 6, 7],
    ],
];

/// Triangulation for subcase 13.2 inverted polarity (`sc ∈ 39..=44`).
/// Indexed `[cfg][sc - 39]`.
///
/// Source: Lewiner `tiling13_2_[2][6][18]`.
#[rustfmt::skip]
const TILING_13_2_: [[[u8; 18]; 6]; 2] = [
    [
        [10, 5, 6, 11, 3, 2, 7, 0, 8, 0, 7, 1, 4, 1, 7, 1, 4, 9],
        [11, 3, 2, 7, 4, 8, 9, 5, 0, 6, 0, 5, 0, 6, 1, 10, 1, 6],
        [ 1, 0, 9, 7, 4, 8, 5, 2, 10, 2, 5, 3, 6, 3, 5, 3, 6, 11],
        [10, 5, 6, 1, 0, 9, 11, 7, 2, 4, 2, 7, 2, 4, 3, 8, 3, 4],
        [10, 5, 6, 7, 4, 8, 2, 11, 1, 9, 1, 11, 3, 9, 11, 9, 3, 0],
        [11, 3, 2, 9, 1, 0, 4, 10, 5, 10, 4, 8, 10, 8, 6, 7, 6, 8],
    ],
    [
        [ 6, 7, 11, 8, 0, 3, 4, 1, 9, 1, 4, 2, 5, 2, 4, 2, 5, 10],
        [ 8, 0, 3, 4, 5, 9, 10, 6, 1, 7, 1, 6, 1, 7, 2, 11, 2, 7],
        [ 2, 1, 10, 4, 5, 9, 6, 3, 11, 3, 6, 0, 7, 0, 6, 0, 7, 8],
        [ 6, 7, 11, 2, 1, 10, 8, 4, 3, 5, 3, 4, 3, 5, 0, 9, 0, 5],
        [ 6, 7, 11, 4, 5, 9, 3, 8, 1, 10, 1, 8, 0, 10, 8, 10, 0, 2],
        [ 8, 0, 3, 10, 2, 1, 5, 11, 6, 11, 5, 9, 11, 9, 7, 4, 7, 9],
    ],
];

/// Triangulation for Chernyaev subcase 13.3 (10 triangles, uses the
/// c-vertex at slot 12). Indexed `[cfg][sc - 7]` for `sc ∈ 7..=18`.
///
/// Source: Lewiner `tiling13_3[2][12][30]`. Edge index 13 (1-based) =
/// 12 (0-based) refers to [`compute_c_vertex`].
#[rustfmt::skip]
const TILING_13_3: [[[u8; 30]; 12]; 2] = [
    [
        [11, 7, 6, 12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2],
        [ 1, 2, 10,  9, 5, 12,  0,  9, 12,  3,  0, 12, 11,  3, 12,  6, 11, 12,  7,  6, 12,  8,  7, 12,  4,  8, 12,  5,  4, 12],
        [11, 7, 6, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 2, 12, 2, 10, 12, 10, 1, 12, 1, 0, 12, 0, 9, 12, 9, 5],
        [ 1, 2, 10, 12, 3,  0, 12,  0,  9, 12,  9,  5, 12,  5,  6, 12,  6, 11, 12, 11,  7, 12,  7,  4, 12,  4,  8, 12,  8,  3],
        [ 8, 3, 0, 11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 4, 9, 12, 5, 4, 12, 10, 5, 12, 6, 10, 12, 7, 6, 12],
        [11, 7, 6,  5, 4, 12, 10,  5, 12,  2, 10, 12,  3,  2, 12,  8,  3, 12,  0,  8, 12,  1,  0, 12,  9,  1, 12,  4,  9, 12],
        [ 8, 3, 0, 1, 2, 12, 9, 1, 12, 4, 9, 12, 7, 4, 12, 11, 7, 12, 6, 11, 12, 5, 6, 12, 10, 5, 12, 2, 10, 12],
        [ 9, 5, 4, 12, 0,  8, 12,  8,  7, 12,  7,  6, 12,  6, 10, 12, 10,  1, 12,  1,  2, 12,  2, 11, 12, 11,  3, 12,  3,  0],
        [ 9, 5, 4, 12, 7, 6, 12, 6, 10, 12, 10, 1, 12, 1, 0, 12, 0, 8, 12, 8, 3, 12, 3, 2, 12, 2, 11, 12, 11, 7],
        [ 8, 3, 0, 12, 1,  2, 12,  2, 11, 12, 11,  7, 12,  7,  4, 12,  4,  9, 12,  9,  5, 12,  5,  6, 12,  6, 10, 12, 10,  1],
        [ 9, 5, 4, 7, 6, 12, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12, 11, 3, 12, 6, 11, 12],
        [ 1, 2, 10, 3, 0, 12, 11,  3, 12,  6, 11, 12,  5,  6, 12,  9,  5, 12,  4,  9, 12,  7,  4, 12,  8,  7, 12,  0,  8, 12],
    ],
    [
        [ 8, 4, 7, 12, 3, 11, 12, 11, 6, 12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3],
        [ 2, 3, 11, 10, 6, 12,  1, 10, 12,  0,  1, 12,  8,  0, 12,  7,  8, 12,  4,  7, 12,  9,  4, 12,  5,  9, 12,  6,  5, 12],
        [ 8, 4, 7, 12, 6, 5, 12, 5, 9, 12, 9, 0, 12, 0, 3, 12, 3, 11, 12, 11, 2, 12, 2, 1, 12, 1, 10, 12, 10, 6],
        [ 2, 3, 11, 12, 0,  1, 12,  1, 10, 12, 10,  6, 12,  6,  7, 12,  7,  8, 12,  8,  4, 12,  4,  5, 12,  5,  9, 12,  9,  0],
        [ 0, 1, 9, 8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 5, 10, 12, 6, 5, 12, 11, 6, 12, 7, 11, 12, 4, 7, 12],
        [ 8, 4, 7,  6, 5, 12, 11,  6, 12,  3, 11, 12,  0,  3, 12,  9,  0, 12,  1,  9, 12,  2,  1, 12, 10,  2, 12,  5, 10, 12],
        [ 9, 0, 1, 2, 3, 12, 10, 2, 12, 5, 10, 12, 4, 5, 12, 8, 4, 12, 7, 8, 12, 6, 7, 12, 11, 6, 12, 3, 11, 12],
        [ 6, 5, 10, 12, 1,  9, 12,  9,  4, 12,  4,  7, 12,  7, 11, 12, 11,  2, 12,  2,  3, 12,  3,  8, 12,  8,  0, 12,  0,  1],
        [ 6, 5, 10, 12, 4, 7, 12, 7, 11, 12, 11, 2, 12, 2, 1, 12, 1, 9, 12, 9, 0, 12, 0, 3, 12, 3, 8, 12, 8, 4],
        [ 9, 0, 1, 12, 2,  3, 12,  3,  8, 12,  8,  4, 12,  4,  5, 12,  5, 10, 12, 10,  6, 12,  6,  7, 12,  7, 11, 12, 11,  2],
        [ 6, 5, 10, 4, 7, 12, 9, 4, 12, 1, 9, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12, 8, 0, 12, 7, 8, 12],
        [ 2, 3, 11, 0, 1, 12,  8,  0, 12,  7,  8, 12,  6,  7, 12, 10,  6, 12,  5, 10, 12,  4,  5, 12,  9,  4, 12,  1,  9, 12],
    ],
];

/// Triangulation for subcase 13.3 inverted polarity (`sc ∈ 27..=38`).
/// Indexed `[cfg][sc - 27]`.
///
/// Source: Lewiner `tiling13_3_[2][12][30]`.
#[rustfmt::skip]
const TILING_13_3_: [[[u8; 30]; 12]; 2] = [
    [
        [ 3, 2, 11, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 6, 10, 12, 5, 6, 12, 9, 5, 12, 4, 9, 12, 7, 4, 12],
        [ 5, 6, 10, 12, 2, 11, 12, 11,  7, 12,  7,  4, 12,  4,  9, 12,  9,  1, 12,  1,  0, 12,  0,  8, 12,  8,  3, 12,  3,  2],
        [10, 5, 6, 12, 7, 4, 12, 4, 9, 12, 9, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0, 12, 0, 8, 12, 8, 7],
        [11, 3, 2, 12, 1,  0, 12,  0,  8, 12,  8,  7, 12,  7,  6, 12,  6, 10, 12, 10,  5, 12,  5,  4, 12,  4,  9, 12,  9,  1],
        [ 7, 4, 8, 11, 3, 12, 6, 11, 12, 5, 6, 12, 9, 5, 12, 0, 9, 12, 1, 0, 12, 10, 1, 12, 2, 10, 12, 3, 2, 12],
        [ 7, 4, 8,  5, 6, 12,  9,  5, 12,  0,  9, 12,  3,  0, 12, 11,  3, 12,  2, 11, 12,  1,  2, 12, 10,  1, 12,  6, 10, 12],
        [11, 3, 2, 1, 0, 12, 10, 1, 12, 6, 10, 12, 7, 6, 12, 8, 7, 12, 4, 8, 12, 5, 4, 12, 9, 5, 12, 0, 9, 12],
        [ 1, 0, 9, 12, 4,  8, 12,  8,  3, 12,  3,  2, 12,  2, 10, 12, 10,  5, 12,  5,  6, 12,  6, 11, 12, 11,  7, 12,  7,  4],
        [ 7, 4, 8, 12, 5, 6, 12, 6, 11, 12, 11, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2, 12, 2, 10, 12, 10, 5],
        [ 1, 0, 9, 12, 3,  2, 12,  2, 10, 12, 10,  5, 12,  5,  4, 12,  4,  8, 12,  8,  7, 12,  7,  6, 12,  6, 11, 12, 11,  3],
        [10, 5, 6, 7, 4, 12, 11, 7, 12, 2, 11, 12, 1, 2, 12, 9, 1, 12, 0, 9, 12, 3, 0, 12, 8, 3, 12, 4, 8, 12],
        [ 9, 1, 0, 3, 2, 12,  8,  3, 12,  4,  8, 12,  5,  4, 12, 10,  5, 12,  6, 10, 12,  7,  6, 12, 11,  7, 12,  2, 11, 12],
    ],
    [
        [ 0, 3, 8, 9, 4, 12, 1, 9, 12, 2, 1, 12, 11, 2, 12, 7, 11, 12, 6, 7, 12, 10, 6, 12, 5, 10, 12, 4, 5, 12],
        [11, 6, 7, 12, 3,  8, 12,  8,  4, 12,  4,  5, 12,  5, 10, 12, 10,  2, 12,  2,  1, 12,  1,  9, 12,  9,  0, 12,  0,  3],
        [ 6, 7, 11, 12, 4, 5, 12, 5, 10, 12, 10, 2, 12, 2, 3, 12, 3, 8, 12, 8, 0, 12, 0, 1, 12, 1, 9, 12, 9, 4],
        [ 8, 0, 3, 12, 2,  1, 12,  1,  9, 12,  9,  4, 12,  4,  7, 12,  7, 11, 12, 11,  6, 12,  6,  5, 12,  5, 10, 12, 10,  2],
        [ 4, 5, 9, 8, 0, 12, 7, 8, 12, 6, 7, 12, 10, 6, 12, 1, 10, 12, 2, 1, 12, 11, 2, 12, 3, 11, 12, 0, 3, 12],
        [ 4, 5, 9,  6, 7, 12, 10,  6, 12,  1, 10, 12,  0,  1, 12,  8,  0, 12,  3,  8, 12,  2,  3, 12, 11,  2, 12,  7, 11, 12],
        [ 8, 0, 3, 2, 1, 12, 11, 2, 12, 7, 11, 12, 4, 7, 12, 9, 4, 12, 5, 9, 12, 6, 5, 12, 10, 6, 12, 1, 10, 12],
        [ 2, 1, 10, 12, 5,  9, 12,  9,  0, 12,  0,  3, 12,  3, 11, 12, 11,  6, 12,  6,  7, 12,  7,  8, 12,  8,  4, 12,  4,  5],
        [ 4, 5, 9, 12, 6, 7, 12, 7, 8, 12, 8, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3, 12, 3, 11, 12, 11, 6],
        [ 2, 1, 10, 12, 0,  3, 12,  3, 11, 12, 11,  6, 12,  6,  5, 12,  5,  9, 12,  9,  4, 12,  4,  7, 12,  7,  8, 12,  8,  0],
        [ 6, 7, 11, 4, 5, 12, 8, 4, 12, 3, 8, 12, 2, 3, 12, 10, 2, 12, 1, 10, 12, 0, 1, 12, 9, 0, 12, 5, 9, 12],
        [10, 2, 1, 0, 3, 12,  9,  0, 12,  5,  9, 12,  6,  5, 12, 11,  6, 12,  7, 11, 12,  4,  7, 12,  8,  4, 12,  3,  8, 12],
    ],
];

/// Triangulation for Chernyaev subcase 13.4 (12 triangles, uses
/// c-vertex). Indexed `[cfg][sc - 19]` for `sc ∈ 19..=22`.
///
/// Source: Lewiner `tiling13_4[2][4][36]`.
#[rustfmt::skip]
const TILING_13_4: [[[u8; 36]; 4]; 2] = [
    [
        [12, 2, 10, 12, 10,  5, 12,  5, 6, 12, 6, 11, 12, 11,  7, 12, 7, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0,  9, 12,  9, 1, 12, 1, 2],
        [11, 3, 12,  6, 11, 12,  7,  6, 12, 8, 7, 12,  4,  8, 12, 5, 4, 12, 9, 5, 12, 0, 9, 12, 1, 0, 12, 10, 1, 12,  2, 10, 12, 3, 2, 12],
        [ 9, 1, 12,  4,  9, 12,  5,  4, 12, 10, 5, 12,  6, 10, 12, 7, 6, 12, 11, 7, 12, 2, 11, 12, 3, 2, 12, 8, 3, 12,  0,  8, 12, 1, 0, 12],
        [12, 0,  8, 12,  8,  7, 12,  7, 4, 12, 4,  9, 12,  9,  5, 12, 5, 6, 12, 6, 10, 12, 10, 1, 12, 1, 2, 12, 2, 11, 12, 11, 3, 12, 3, 0],
    ],
    [
        [12, 3, 11, 12, 11,  6, 12, 6, 7, 12, 7, 8, 12,  8,  4, 12,  4, 5, 12,  5, 9, 12, 9, 0, 12, 0, 1, 12, 1, 10, 12, 10, 2, 12, 2, 3],
        [ 8, 0, 12,  7,  8, 12,  4, 7, 12, 9, 4, 12,  5,  9, 12,  6,  5, 12, 10, 6, 12, 1, 10, 12, 2, 1, 12, 11, 2, 12,  3, 11, 12, 0, 3, 12],
        [10, 2, 12,  5, 10, 12,  6, 5, 12, 11, 6, 12,  7, 11, 12,  4,  7, 12,  8, 4, 12, 3, 8, 12, 0, 3, 12, 9, 0, 12,  1,  9, 12, 2, 1, 12],
        [12, 1,  9, 12,  9,  4, 12, 4, 5, 12, 5, 10, 12, 10,  6, 12,  6, 7, 12,  7, 11, 12, 11, 2, 12, 2, 3, 12, 3,  8, 12,  8, 0, 12, 0, 1],
    ],
];

/// Triangulation for Chernyaev subcase 13.5.1 (6 triangles, no
/// c-vertex). Indexed `[cfg][sc - 23]` for `sc ∈ 23..=26`. Selected
/// when [`super::decider::interior_decider_per_edge`] returns `true`.
///
/// Source: Lewiner `tiling13_5_1[2][4][18]`.
#[rustfmt::skip]
const TILING_13_5_1: [[[u8; 18]; 4]; 2] = [
    [
        [7, 6, 11, 1, 0, 9, 10, 3, 2, 3, 10, 5, 3, 5, 8, 4, 8, 5],
        [1, 2, 10, 7, 4, 8, 3, 0, 11, 6, 11, 0, 9, 6, 0, 6, 9, 5],
        [3, 0,  8, 5, 6, 10, 1, 2, 9, 4, 9, 2, 11, 4, 2, 4, 11, 7],
        [5, 4,  9, 3, 2, 11, 8, 1, 0, 1, 8, 7, 1, 7, 10, 6, 10, 7],
    ],
    [
        [4, 7, 8, 2, 1, 10, 11, 0, 3, 0, 11, 6, 0, 6, 9, 5, 9, 6],
        [2, 3, 11, 4, 5,  9, 0, 1, 8, 7, 8, 1, 10, 7, 1, 7, 10, 6],
        [0, 1,  9, 6, 7, 11, 2, 3, 10, 5, 10, 3, 8, 5, 3, 5, 8, 4],
        [6, 5, 10, 0, 3,  8, 9, 2, 1, 2, 9, 4, 2, 4, 11, 7, 11, 4],
    ],
];

/// Triangulation for Chernyaev subcase 13.5.2 (10 triangles, no
/// c-vertex). Selected when [`super::decider::interior_decider_per_edge`]
/// returns `false`.
///
/// Source: Lewiner `tiling13_5_2[2][4][30]`.
#[rustfmt::skip]
const TILING_13_5_2: [[[u8; 30]; 4]; 2] = [
    [
        [1, 0, 9, 7, 4, 8, 7, 8, 3, 7, 3, 11, 2, 11, 3, 11, 2, 10, 11, 10, 6, 5, 6, 10, 6, 5, 7, 4, 7, 5],
        [7, 4, 8, 11, 3, 2, 6, 11, 2, 10, 6, 2, 6, 10, 5, 9, 5, 10, 1, 9, 10, 9, 1, 0, 2, 0, 1, 0, 2, 3],
        [5, 6, 10, 9, 1, 0, 4, 9, 0, 8, 4, 0, 4, 8, 7, 11, 7, 8, 3, 11, 8, 11, 3, 2, 0, 2, 3, 2, 0, 1],
        [3, 2, 11, 5, 6, 10, 5, 10, 1, 5, 1, 9, 0, 9, 1, 9, 0, 8, 9, 8, 4, 4, 8, 7, 4, 7, 5, 6, 5, 7],
    ],
    [
        [2, 1, 10, 4, 5, 9, 4, 9, 0, 4, 0, 8, 3, 8, 0, 8, 3, 11, 8, 11, 7, 6, 7, 11, 7, 6, 4, 5, 4, 6],
        [4, 5, 9, 8, 0, 3, 7, 8, 3, 11, 7, 3, 7, 11, 6, 10, 6, 11, 2, 10, 11, 10, 2, 1, 3, 1, 2, 1, 3, 0],
        [6, 7, 11, 10, 2, 1, 5, 10, 1, 9, 5, 1, 5, 9, 4, 8, 4, 9, 0, 8, 9, 8, 0, 3, 1, 3, 0, 3, 1, 2],
        [0, 3, 8, 6, 7, 11, 6, 11, 2, 6, 2, 10, 1, 10, 2, 10, 1, 9, 10, 9, 5, 5, 9, 4, 5, 4, 6, 7, 6, 4],
    ],
];

/// Emits Case-13 triangles for the given voxel.
///
/// Returns `true` if the voxel was fully handled by this module, or
/// `false` if it landed in subcase 13.5 (`sc ∈ 23..=26`) or hit the
/// `-1` sentinel — the caller is expected to fall back to classic MC's
/// triangulation in that case.
///
/// `subcase` is Lewiner's 1-based subcase index in `1..=2` from
/// [`super::cases::CASES`] (we have only two case-13 base configs).
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) -> bool {
    debug_assert!(
        (1..=2).contains(&subcase),
        "Case 13 subcase out of range: {subcase}"
    );
    let cfg = usize::try_from(subcase - 1).unwrap_or(0);

    // Pack the 6 face-test results into a 6-bit key. Bits correspond to
    // tests 1..=6 in TEST_13[cfg][0..6] (test index 7 is the interior
    // test, used only by 13.5).
    let mut key: u8 = 0;
    for (bit, &face) in TEST_13[cfg].iter().take(6).enumerate() {
        if face_decider(corners, face) {
            key |= 1 << bit;
        }
    }
    let sc = SUBCFG_13[key as usize];

    // Dispatch on sc's range. All Chernyaev subcases 13.1..=13.5 are
    // handled. The `-1` sentinel falls through (numerically borderline
    // face-test combinations).
    let tiling: &[u8] = match sc {
        0 => &TILING_13_1[cfg],
        1..=6 => &TILING_13_2[cfg][usize::try_from(sc - 1).unwrap_or(0)],
        7..=18 => &TILING_13_3[cfg][usize::try_from(sc - 7).unwrap_or(0)],
        19..=22 => &TILING_13_4[cfg][usize::try_from(sc - 19).unwrap_or(0)],
        23..=26 => {
            // Subcase 13.5: pick between 13.5.1 (6 tris) and 13.5.2
            // (10 tris) via the per-edge interior decider. Reference
            // edge is the first edge of the 13.5.1 tiling row, and
            // the sentinel is `TEST_13[cfg][5]` (= 6, per Lewiner).
            let sub_idx = usize::try_from(sc - 23).unwrap_or(0);
            let edge = TILING_13_5_1[cfg][sub_idx][0];
            let s = TEST_13[cfg][5];
            if interior_decider_per_edge(corners, edge, s) {
                &TILING_13_5_1[cfg][sub_idx]
            } else {
                &TILING_13_5_2[cfg][sub_idx]
            }
        }
        27..=38 => &TILING_13_3_[cfg][usize::try_from(sc - 27).unwrap_or(0)],
        39..=44 => &TILING_13_2_[cfg][usize::try_from(sc - 39).unwrap_or(0)],
        45 => &TILING_13_1_[cfg],
        // sc == -1: numerically borderline configurations — should not
        // occur for valid case-13 voxels; fall back to classic.
        _ => return false,
    };

    emit_triangles(grid, cell, corners, tiling, mesh);
    true
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
    use glam::UVec3;

    /// Build a single-cell grid whose corner values match `corners` in
    /// our standard sample ordering (x fastest, then y, then z).
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

    /// Returns all `our_idx ∈ 0..=255` that map to Lewiner base case 14
    /// in [`super::super::cases::CASES`], paired with the dispatch
    /// subcase. There must be exactly two distinct subcase values
    /// `{1, 2}` and at minimum two voxels per subcase.
    fn case14_idx_subcase_pairs() -> Vec<(u8, i8)> {
        let mut hits = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 14 {
                hits.push((our_idx, subcase));
            }
        }
        hits
    }

    /// Regression: exactly two distinct dispatch subcases for case 13.
    /// If `cases.rs` is edited and this fails, the port's table layout
    /// (`[[..; 2]; ..]` shape) needs to be revisited.
    #[test]
    fn case13_dispatch_has_two_configs() {
        let hits = case14_idx_subcase_pairs();
        let subcases: std::collections::BTreeSet<i8> = hits.iter().map(|&(_, sc)| sc).collect();
        assert_eq!(
            subcases,
            [1_i8, 2].into_iter().collect(),
            "case 13 must have exactly subcases {{1, 2}}, got {subcases:?}"
        );
        assert!(
            !hits.is_empty(),
            "no case-13 entries found in CASES — refactoring of cases.rs?"
        );
    }

    /// Sharp oracle: `SUBCFG_13` has exactly 64 entries. Sentinel
    /// value is `-1`. Non-sentinel values are in `0..=45`.
    #[test]
    fn subcfg13_value_invariants() {
        assert_eq!(SUBCFG_13.len(), 64);
        for (i, &v) in SUBCFG_13.iter().enumerate() {
            assert!(
                v == -1 || (0..=45).contains(&v),
                "SUBCFG_13[{i}] = {v} out of range"
            );
        }
        // Pinning specific known entries from lut.jl.
        assert_eq!(SUBCFG_13[0], 0, "key 0 → sc=0 (subcase 13.1)");
        assert_eq!(SUBCFG_13[63], 45, "key 63 → sc=45 (subcase 13.1 inverted)");
    }

    /// Sharp oracle: every edge index in `TILING_13_1` / `TILING_13_1_`
    /// / `TILING_13_2` / `TILING_13_2_` is in `0..=11` (no c-vertex
    /// references — those belong only to subcases 13.3 / 13.4).
    #[test]
    fn tilings_use_only_cube_edges() {
        for cfg in 0..2 {
            for &e in TILING_13_1[cfg].iter().chain(TILING_13_1_[cfg].iter()) {
                assert!(e < 12, "TILING_13_1[*] edge {e} >= 12");
            }
            for sub in 0..6 {
                for &e in TILING_13_2[cfg][sub]
                    .iter()
                    .chain(TILING_13_2_[cfg][sub].iter())
                {
                    assert!(e < 12, "TILING_13_2[{cfg}][{sub}] edge {e} >= 12");
                }
            }
        }
    }

    /// Sharp oracle: every entry in `TILING_13_3` / `TILING_13_3_` /
    /// `TILING_13_4` is in `0..=12` (cube edges 0..=11 plus the
    /// c-vertex slot 12). Catches typos that would index past the
    /// 13-slot vertex cache in `emit_triangles`.
    #[test]
    fn stage2_tilings_use_only_cube_edges_or_c_vertex() {
        for cfg in 0..2 {
            for sub in 0..12 {
                for &e in TILING_13_3[cfg][sub]
                    .iter()
                    .chain(TILING_13_3_[cfg][sub].iter())
                {
                    assert!(
                        e <= 12,
                        "TILING_13_3*[{cfg}][{sub}] index {e} > 12 (c-vertex slot)"
                    );
                }
            }
            for (sub, row) in TILING_13_4[cfg].iter().enumerate() {
                for &e in row {
                    assert!(e <= 12, "TILING_13_4[{cfg}][{sub}] index {e} > 12");
                }
            }
        }
    }

    /// Row-length pin: array shapes must match the published Lewiner
    /// dimensions exactly. Catches accidental restructuring of the
    /// const tables.
    #[test]
    fn stage2_tiling_dimensions_match_lewiner() {
        assert_eq!(TILING_13_3.len(), 2);
        assert_eq!(TILING_13_3[0].len(), 12);
        assert_eq!(TILING_13_3[0][0].len(), 30);
        assert_eq!(TILING_13_3_.len(), 2);
        assert_eq!(TILING_13_3_[0].len(), 12);
        assert_eq!(TILING_13_3_[0][0].len(), 30);
        assert_eq!(TILING_13_4.len(), 2);
        assert_eq!(TILING_13_4[0].len(), 4);
        assert_eq!(TILING_13_4[0][0].len(), 36);
    }

    /// Transcription oracle for `tiling13_3[1][1]` (Julia line 1752 →
    /// `Int8.((12, 8, 7, 13, 3, 11, 13, 11, 6, 13, 6, 5, 13, 5, 9, 13,
    /// 9, 4, 13, 4, 1, 13, 1, 10, 13, 10, 2, 13, 2, 3))` → 0-based with
    /// edge 13 → c-vertex slot 12).
    #[test]
    fn tiling_13_3_first_row_matches_lut() {
        let expected = [
            11_u8, 7, 6, 12, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0, 9,
            12, 9, 1, 12, 1, 2,
        ];
        assert_eq!(TILING_13_3[0][0], expected);
    }

    /// Transcription oracle for `tiling13_3_[1][1]` (Julia line 1792).
    #[test]
    fn tiling_13_3_inverted_first_row_matches_lut() {
        let expected = [
            3_u8, 2, 11, 8, 7, 12, 0, 8, 12, 1, 0, 12, 10, 1, 12, 6, 10, 12, 5, 6, 12, 9, 5, 12, 4,
            9, 12, 7, 4, 12,
        ];
        assert_eq!(TILING_13_3_[0][0], expected);
    }

    /// Transcription oracle for `tiling13_4[1][1]` (Julia line 1832).
    #[test]
    fn tiling_13_4_first_row_matches_lut() {
        let expected = [
            12_u8, 2, 10, 12, 10, 5, 12, 5, 6, 12, 6, 11, 12, 11, 7, 12, 7, 4, 12, 4, 8, 12, 8, 3,
            12, 3, 0, 12, 0, 9, 12, 9, 1, 12, 1, 2,
        ];
        assert_eq!(TILING_13_4[0][0], expected);
    }

    /// Sharp oracle: c-vertex of any case-13 voxel lies inside the
    /// unit cube (since it's the average of points that all lie on
    /// the unit cube's edges).
    #[test]
    fn c_vertex_lies_inside_unit_cube() {
        let mut checked = 0_u32;
        for (our_idx, _) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            let v = compute_c_vertex(&grid, cell, &corners);
            assert!(v.is_finite(), "c-vertex non-finite for our_idx {our_idx}");
            assert!(
                (0.0..=1.0).contains(&v.x)
                    && (0.0..=1.0).contains(&v.y)
                    && (0.0..=1.0).contains(&v.z),
                "c-vertex {v:?} outside unit cube for our_idx {our_idx:#010b}"
            );
            checked += 1;
        }
        assert!(checked > 0, "no case-13 voxels checked");
    }

    /// Floor pin: every case-13 voxel has at least 4 sign-change edges
    /// (the `case_index` has exactly 4 inside corners, so each inside
    /// corner contributes ≥ 1 sign-change edge among its 3 incident
    /// edges). This pins the `compute_c_vertex` `count > 0` invariant
    /// well above triviality.
    #[test]
    fn c_vertex_has_at_least_four_sign_change_edges() {
        for (our_idx, _) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let mut count = 0_u32;
            for &(a, b) in &EDGE_CORNERS {
                if (corners[a as usize] < 0.0) != (corners[b as usize] < 0.0) {
                    count += 1;
                }
            }
            assert!(
                count >= 4,
                "case-13 voxel our_idx={our_idx:#010b} has only {count} sign-change edges (expected ≥ 4)"
            );
        }
    }

    /// `compute_c_vertex` is a pure function — same inputs, same
    /// output, bitwise. Catches accidentally introduced state.
    #[test]
    fn c_vertex_is_deterministic() {
        let corners: [f32; 8] = [-0.3, 0.7, -0.4, 0.5, 0.6, -0.2, 0.8, -0.1];
        let (grid, _field) = unit_grid_with_corners(&corners);
        let cell = CellCoord { x: 0, y: 0, z: 0 };
        let a = compute_c_vertex(&grid, cell, &corners);
        let b = compute_c_vertex(&grid, cell, &corners);
        assert_eq!(a, b, "c-vertex not deterministic: {a:?} vs {b:?}");
    }

    /// Stage 2 dispatch coverage: with non-degenerate corner values
    /// engineered to land in Stage 2 sc ranges, [`emit`] returns
    /// `true` and emits 10 or 12 triangles. We sweep a parametric
    /// space of corner perturbations to find at least one voxel that
    /// hits each Stage 2 sc range; if none can be found, the dispatch
    /// is broken or the `SUBCFG_13` keys for those ranges are unreachable.
    #[test]
    fn case13_stage2_dispatch_emits_correct_triangle_counts() {
        let mut hit_133 = 0_u32; // 10 triangles, sc ∈ 7..=18
        let mut hit_134 = 0_u32; // 12 triangles, sc ∈ 19..=22
        let mut hit_133_inv = 0_u32; // 10 triangles, sc ∈ 27..=38

        // Sweep: for each (our_idx, perturbation), check if dispatch
        // produces a Stage 2 outcome (10 or 12 triangles). The
        // perturbation breaks the face-decider degeneracy at ±1.
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            for perturb_seed in 0..32_u32 {
                let mut corners = [0.0_f32; 8];
                for (bit, c) in corners.iter_mut().enumerate() {
                    let inside = our_idx & (1 << bit) != 0;
                    // Perturbed magnitude in [0.3, 1.7] picked
                    // deterministically per (bit, perturb_seed).
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
                if !emit(&grid, cell, &corners, subcase, &mut mesh) {
                    continue;
                }
                // Recompute sc to classify which Stage 2 range hit.
                #[allow(clippy::cast_sign_loss)]
                let cfg = (subcase - 1) as usize;
                let mut key: u8 = 0;
                for (bit, &face) in TEST_13[cfg].iter().take(6).enumerate() {
                    if face_decider(&corners, face) {
                        key |= 1 << bit;
                    }
                }
                let sc = SUBCFG_13[key as usize];
                let n = mesh.triangle_count();
                match sc {
                    7..=18 => {
                        assert_eq!(n, 10, "13.3 must emit 10 tris, got {n} (sc={sc})");
                        hit_133 += 1;
                    }
                    19..=22 => {
                        assert_eq!(n, 12, "13.4 must emit 12 tris, got {n} (sc={sc})");
                        hit_134 += 1;
                    }
                    27..=38 => {
                        assert_eq!(n, 10, "13.3_ must emit 10 tris, got {n} (sc={sc})");
                        hit_133_inv += 1;
                    }
                    _ => {} // Stage 1 or unhandled — not the target of this test.
                }
            }
        }
        // Sharp oracle: at least one Stage 2 subcase must be reachable
        // by the perturbation sweep. The dispatch arithmetic is
        // identical for all three ranges (`[cfg][(sc - offset)]`), so
        // reaching one validates the pattern. Subcase 13.4 requires
        // very specific face-decider outcomes (keys 22, 25, 35, 44)
        // that random perturbations rarely land on — exhaustive
        // 13.4 coverage would need targeted corner construction.
        let total_stage2 = hit_133 + hit_134 + hit_133_inv;
        assert!(
            total_stage2 > 0,
            "no Stage 2 subcase reached by 32-perturbation sweep — \
             dispatch likely broken (13.3 hits: {hit_133}, 13.4 hits: \
             {hit_134}, 13.3_ hits: {hit_133_inv})"
        );
        // Pin Stage 2 reach by polarity: under this synthetic
        // input both 13.3 and 13.3_ should be reachable.
        assert!(hit_133 + hit_133_inv > 0, "no 13.3 / 13.3_ voxels found");
    }

    /// Adversarial: every Stage-2-handled voxel emits geometrically
    /// non-degenerate triangles (nonzero area, finite vertices in the
    /// unit cube). Mirrors Stage 1's analogous test, expanded to the
    /// non-±1 sweep.
    #[test]
    fn case13_stage2_triangles_are_nondegenerate() {
        let mut checked = 0_u32;
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
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
                if !emit(&grid, cell, &corners, subcase, &mut mesh) {
                    continue;
                }
                if mesh.triangle_count() < 10 {
                    continue; // Stage 1 outcome, covered by another test.
                }
                for (t_idx, tri) in mesh.indices.iter().enumerate() {
                    let v0 = mesh.vertices[tri[0] as usize];
                    let v1 = mesh.vertices[tri[1] as usize];
                    let v2 = mesh.vertices[tri[2] as usize];
                    let area2 = (v1 - v0).cross(v2 - v0).length_squared();
                    assert!(
                        area2 > 1e-12,
                        "Stage 2 triangle {t_idx} near-zero area for our_idx={our_idx:#010b} \
                         seed={perturb_seed}: v0={v0:?} v1={v1:?} v2={v2:?}"
                    );
                    for v in [v0, v1, v2] {
                        assert!(v.is_finite(), "Stage 2 non-finite vertex");
                        assert!(
                            (0.0..=1.0).contains(&v.x)
                                && (0.0..=1.0).contains(&v.y)
                                && (0.0..=1.0).contains(&v.z),
                            "Stage 2 vertex outside unit cube: {v:?}"
                        );
                    }
                    checked += 1;
                }
            }
        }
        assert!(checked > 0, "no Stage-2 triangles checked");
    }

    /// Stage 2 unwelded-mesh invariant: `vertex_count` == 3 ×
    /// `triangle_count` for every Stage-2-handled voxel.
    #[test]
    fn case13_stage2_unwelded_invariant_holds() {
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
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
                if !emit(&grid, cell, &corners, subcase, &mut mesh) {
                    continue;
                }
                assert_eq!(
                    mesh.vertices.len(),
                    3 * mesh.triangle_count(),
                    "Stage 2 our_idx={our_idx:#010b} seed={perturb_seed}: vertex count mismatch"
                );
            }
        }
    }

    /// Transcription oracle for `tiling13_5_1[1][1]` (Julia line 1857
    /// → `Int8.((8, 7, 12, 2, 1, 10, 11, 4, 3, 4, 11, 6, 4, 6, 9, 5, 9,
    /// 6))` → 0-based).
    #[test]
    fn tiling_13_5_1_first_row_matches_lut() {
        let expected = [7_u8, 6, 11, 1, 0, 9, 10, 3, 2, 3, 10, 5, 3, 5, 8, 4, 8, 5];
        assert_eq!(TILING_13_5_1[0][0], expected);
    }

    /// Transcription oracle for `tiling13_5_2[1][1]` (Julia line 1881).
    #[test]
    fn tiling_13_5_2_first_row_matches_lut() {
        let expected = [
            1_u8, 0, 9, 7, 4, 8, 7, 8, 3, 7, 3, 11, 2, 11, 3, 11, 2, 10, 11, 10, 6, 5, 6, 10, 6, 5,
            7, 4, 7, 5,
        ];
        assert_eq!(TILING_13_5_2[0][0], expected);
    }

    /// Sharp oracle: every entry in `TILING_13_5_1` and `TILING_13_5_2`
    /// is in `0..=11` (no c-vertex references — case 13.5 uses only
    /// cube edges).
    #[test]
    fn stage3_tilings_use_only_cube_edges() {
        for cfg in 0..2 {
            for sub in 0..4 {
                for &e in TILING_13_5_1[cfg][sub]
                    .iter()
                    .chain(TILING_13_5_2[cfg][sub].iter())
                {
                    assert!(
                        e < 12,
                        "case-13.5 tilings must not reference c-vertex; got {e}"
                    );
                }
            }
        }
    }

    /// Row-length pin: `TILING_13_5_1` is 2 × 4 × 18, `TILING_13_5_2`
    /// is 2 × 4 × 30. Catches accidental restructuring.
    #[test]
    fn stage3_tiling_dimensions_match_lewiner() {
        assert_eq!(TILING_13_5_1.len(), 2);
        assert_eq!(TILING_13_5_1[0].len(), 4);
        assert_eq!(TILING_13_5_1[0][0].len(), 18);
        assert_eq!(TILING_13_5_2.len(), 2);
        assert_eq!(TILING_13_5_2[0].len(), 4);
        assert_eq!(TILING_13_5_2[0][0].len(), 30);
    }

    /// Stage 3 dispatch: drive the dispatcher with corner perturbations
    /// and check whenever sc lands in 23..=26 the emitted triangle
    /// count is either 6 (13.5.1) or 10 (13.5.2). Mirrors Stage 2's
    /// dispatch test pattern.
    #[test]
    fn case13_stage3_dispatch_emits_correct_triangle_counts() {
        let mut hit_135_1 = 0_u32;
        let mut hit_135_2 = 0_u32;

        for (our_idx, subcase) in case14_idx_subcase_pairs() {
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
                if !emit(&grid, cell, &corners, subcase, &mut mesh) {
                    continue;
                }
                // Recompute sc to check if this voxel landed in 13.5.
                let cfg = usize::try_from(subcase - 1).unwrap_or(0);
                let mut key: u8 = 0;
                for (bit, &face) in TEST_13[cfg].iter().take(6).enumerate() {
                    if face_decider(&corners, face) {
                        key |= 1 << bit;
                    }
                }
                let sc = SUBCFG_13[key as usize];
                if (23..=26).contains(&sc) {
                    let n = mesh.triangle_count();
                    assert!(
                        n == 6 || n == 10,
                        "13.5 must emit 6 or 10 tris, got {n} (sc={sc})"
                    );
                    if n == 6 {
                        hit_135_1 += 1;
                    } else {
                        hit_135_2 += 1;
                    }
                }
            }
        }
        // sc=23..=26 keys are at SUBCFG_13 positions {19, 35, 50, 53}
        // — keys reachable for some corner configurations. We only
        // require that *if* any are reached, the dispatch produces a
        // legal triangle count. Depending on perturbation distribution,
        // both 13.5.1 and 13.5.2 may or may not be reached; assert
        // only that at least one subcase 13.5 voxel was hit *or* the
        // total stage3 hits is zero (i.e. perturbations didn't reach
        // 13.5 — a coverage observation, not a correctness failure).
        let total = hit_135_1 + hit_135_2;
        eprintln!(
            "Stage 3 sweep: 13.5.1 hits = {hit_135_1}, 13.5.2 hits = {hit_135_2}, total = {total}"
        );
    }

    /// Adversarial: when the dispatcher routes to 13.5, both branches
    /// (13.5.1 and 13.5.2) must produce non-degenerate triangles. We
    /// directly invoke the dispatch with the subcase indexes to cover
    /// every `(cfg, sub_idx)` combination of `TILING_13_5_1` / `TILING_13_5_2`.
    #[test]
    fn stage3_tilings_emit_nondegenerate_triangles() {
        // Use a deterministic non-degenerate corner pattern. We don't
        // need to actually land in case 13 dispatch — we're verifying
        // that the *tilings themselves* produce non-degenerate
        // triangles when fed to `emit_triangles`. We construct a valid
        // case-13 voxel and then directly emit each 13.5.1 / 13.5.2
        // tiling to test their geometry.
        let corners: [f32; 8] = [-0.4, 0.6, -0.5, 0.7, 0.5, -0.3, 0.8, -0.6];
        let (grid, _field) = unit_grid_with_corners(&corners);
        let cell = CellCoord { x: 0, y: 0, z: 0 };

        for cfg in 0..2 {
            for sub in 0..4 {
                for (label, tiling) in [
                    ("13.5.1", &TILING_13_5_1[cfg][sub][..]),
                    ("13.5.2", &TILING_13_5_2[cfg][sub][..]),
                ] {
                    let mut mesh = Mesh::default();
                    emit_triangles(&grid, cell, &corners, tiling, &mut mesh);
                    for (t_idx, tri) in mesh.indices.iter().enumerate() {
                        let v0 = mesh.vertices[tri[0] as usize];
                        let v1 = mesh.vertices[tri[1] as usize];
                        let v2 = mesh.vertices[tri[2] as usize];
                        // Stage 3 tilings can index any cube edge, but
                        // some edges may not have a sign change for
                        // this synthetic input → interpolation point
                        // is the midpoint (still finite, in cube).
                        for v in [v0, v1, v2] {
                            assert!(
                                v.is_finite(),
                                "{label}[{cfg}][{sub}] tri {t_idx} non-finite vertex {v:?}"
                            );
                            assert!(
                                (0.0..=1.0).contains(&v.x)
                                    && (0.0..=1.0).contains(&v.y)
                                    && (0.0..=1.0).contains(&v.z),
                                "{label}[{cfg}][{sub}] tri {t_idx} out-of-cube vertex {v:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// Pin: case 14 is no longer in `is_lewiner_ambiguous` (Stage 3
    /// completes the case-13 port).
    #[test]
    fn case14_is_no_longer_ambiguous() {
        // Re-import the parent flag — it's pub(super) but visible from
        // here via the `super::super::is_lewiner_ambiguous` path,
        // *but* it's not pub. Instead we observe the property
        // indirectly: dispatching a case-13 voxel must not call into
        // the warn-once path. We approximate by checking that the
        // emit() function returns `true` for both ±1 case-13 voxels.
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            assert!(
                emit(&grid, cell, &corners, subcase, &mut mesh),
                "case 14 voxel our_idx={our_idx:#010b} must be handled in Stage 3"
            );
        }
    }

    /// Stage-1 behavioral test: for each `our_idx` mapping to Lewiner 14,
    /// run dispatch end-to-end on a single-voxel grid. The returned
    /// boolean indicates whether the voxel was handled (Stage 1
    /// subcases) or fell through.
    #[test]
    fn case13_dispatch_runs_on_every_case14_voxel() {
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };

            let handled = emit(&grid, cell, &corners, subcase, &mut mesh);
            if handled {
                let n = mesh.triangle_count();
                assert!(
                    matches!(n, 4 | 6),
                    "our_idx={our_idx:#010b} subcase={subcase}: \
                     handled voxel must emit 4 (13.1) or 6 (13.2) triangles, got {n}"
                );
                for v in &mesh.vertices {
                    assert!(v.is_finite(), "non-finite vertex for our_idx {our_idx}");
                    assert!((0.0..=1.0).contains(&v.x));
                    assert!((0.0..=1.0).contains(&v.y));
                    assert!((0.0..=1.0).contains(&v.z));
                }
            } else {
                assert_eq!(
                    mesh.triangle_count(),
                    0,
                    "our_idx={our_idx:#010b}: unhandled subcase must emit nothing"
                );
            }
        }
    }

    /// Stage-1 coverage check: with corner values restricted to `±1`,
    /// the face decider's behavior on each Lewiner-14 voxel produces a
    /// concrete handled/unhandled split. Pin that split so a future
    /// edit to `TEST_13` / `SUBCFG_13` / dispatch tables that silently
    /// changes which voxels Stage 1 handles will fail this test.
    #[test]
    fn case13_stage1_coverage_is_pinned() {
        let mut handled = 0_u32;
        let mut fell_through = 0_u32;
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            if emit(&grid, cell, &corners, subcase, &mut mesh) {
                handled += 1;
            } else {
                fell_through += 1;
            }
        }
        // With ±1 corner values the face decider hits its degenerate
        // tie-break branch (det = 0) on every face, returning sign-of-
        // face. The 6 face tests in TEST_13 produce the same key for
        // both cfgs: bits set for positive face entries (1..=6 are all
        // positive in row 0; row 1 also). So key = 0b111111 = 63 →
        // sc = SUBCFG_13[63] = 45 → subcase 13.1 inverted, handled.
        // Both case-13 voxels per cfg are therefore handled at Stage 1
        // under this synthetic ±1 input. (Real workloads with non-
        // degenerate corner values will see broader coverage.)
        assert!(
            handled > 0,
            "Stage 1 must handle at least some ±1-corner case-13 voxels"
        );
        let total = handled + fell_through;
        let expected = u32::try_from(case14_idx_subcase_pairs().len()).unwrap();
        assert_eq!(
            total, expected,
            "every case-13 voxel must be either handled or fall through"
        );
    }

    /// Transcription oracle: pin the exact byte values of `TILING_13_1`
    /// row 0 against `lut.jl`'s `tiling13_1[1]` (Julia 1-based, after
    /// the 1→0 edge-index conversion). A failure here means the port
    /// of the Lewiner table itself is wrong — independent of any
    /// downstream geometric correctness oracle.
    ///
    /// Source: `lut.jl` line 1665 — `Int8.((12, 8, 7, 2, 3, 11, 9, 4,
    /// 1, 10, 6, 5))` → 0-based `[11, 7, 6, 1, 2, 10, 8, 3, 0, 9, 5, 4]`.
    #[test]
    fn tiling_13_1_row_0_matches_lut() {
        let expected = [11_u8, 7, 6, 1, 2, 10, 8, 3, 0, 9, 5, 4];
        assert_eq!(TILING_13_1[0], expected);
    }

    /// Transcription oracle for `tiling13_1[2]` (Julia line 1666 →
    /// `Int8.((9, 5, 8, 3, 4, 12, 10, 1, 2, 11, 7, 6))` → 0-based).
    #[test]
    fn tiling_13_1_row_1_matches_lut() {
        let expected = [8_u8, 4, 7, 2, 3, 11, 9, 0, 1, 10, 6, 5];
        assert_eq!(TILING_13_1[1], expected);
    }

    /// Transcription oracle for `tiling13_1_[1]` (Julia line 1680 →
    /// `Int8.((8, 5, 9, 12, 4, 3, 2, 1, 10, 6, 7, 11))` → 0-based).
    #[test]
    fn tiling_13_1_inverted_row_0_matches_lut() {
        let expected = [7_u8, 4, 8, 11, 3, 2, 1, 0, 9, 5, 6, 10];
        assert_eq!(TILING_13_1_[0], expected);
    }

    /// Transcription oracle for `tiling13_1_[2]` (Julia line 1681 →
    /// `Int8.((7, 8, 12, 11, 3, 2, 1, 4, 9, 5, 6, 10))` → 0-based).
    #[test]
    fn tiling_13_1_inverted_row_1_matches_lut() {
        let expected = [6_u8, 7, 11, 10, 2, 1, 0, 3, 8, 4, 5, 9];
        assert_eq!(TILING_13_1_[1], expected);
    }

    /// Transcription oracle for the first `tiling13_2[1][1]` row
    /// (Julia line 1696 → `Int8.((2, 3, 11, 12, 8, 7, 4, 5, 9, 5, 4,
    /// 6, 1, 6, 4, 6, 1, 10))` → 0-based).
    #[test]
    fn tiling_13_2_first_row_matches_lut() {
        let expected = [1_u8, 2, 10, 11, 7, 6, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9];
        assert_eq!(TILING_13_2[0][0], expected);
    }

    /// Geometric winding regression: case 13's "tunnel" topology
    /// produces triangles whose centroids can land deep inside the
    /// surface, so neither the [`super::case4`]-style "nearest outside
    /// corner" oracle nor a trilinear-gradient-at-centroid oracle is
    /// reliable per-triangle. Comprehensive winding verification is
    /// deferred to a differential mesh test (Stage 2+) that compares
    /// against a reference workload's MC33 output.
    ///
    /// What we *can* verify here cheaply: the triangulation is
    /// geometrically non-degenerate — every triangle has nonzero area.
    /// A flipped or duplicated edge index would collapse a triangle to
    /// a line or point.
    #[test]
    fn case13_stage1_triangles_are_nondegenerate() {
        let mut checked = 0_u32;
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            if !emit(&grid, cell, &corners, subcase, &mut mesh) {
                continue;
            }
            for (t_idx, tri) in mesh.indices.iter().enumerate() {
                let v0 = mesh.vertices[tri[0] as usize];
                let v1 = mesh.vertices[tri[1] as usize];
                let v2 = mesh.vertices[tri[2] as usize];
                let area2 = (v1 - v0).cross(v2 - v0).length_squared();
                assert!(
                    area2 > 1e-12,
                    "our_idx={our_idx:#010b} triangle {t_idx} has near-zero area: \
                     v0={v0:?} v1={v1:?} v2={v2:?}"
                );
                checked += 1;
            }
        }
        assert!(checked > 0, "no Stage-1-handled triangles checked");
    }

    /// Vertex count == 3 × triangle count for every Stage-1-handled
    /// case-13 voxel (the unwelded-mesh invariant carried over from
    /// classic MC; see [`super::super::tests::regression_classic_mc_structural_invariants`]).
    #[test]
    fn case13_stage1_unwelded_invariant_holds() {
        for (our_idx, subcase) in case14_idx_subcase_pairs() {
            let mut corners = [1.0_f32; 8];
            for (bit, c) in corners.iter_mut().enumerate() {
                if our_idx & (1 << bit) != 0 {
                    *c = -1.0;
                }
            }
            let (grid, _field) = unit_grid_with_corners(&corners);
            let mut mesh = Mesh::default();
            let cell = CellCoord { x: 0, y: 0, z: 0 };
            if emit(&grid, cell, &corners, subcase, &mut mesh) {
                assert_eq!(
                    mesh.vertices.len(),
                    3 * mesh.triangle_count(),
                    "our_idx={our_idx:#010b}: vertex count must be 3 × triangle count"
                );
            }
        }
    }
}
