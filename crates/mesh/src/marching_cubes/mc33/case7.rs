//! Case 7 (Chernyaev) / Lewiner base case 8 — "three non-adjacent corners".
//!
//! The configuration: three corners with the same sign such that no two
//! share a face — the "three diagonals" topology. There are 16 voxel
//! configurations up to rotation (`base_case == 8` rows in
//! [`super::cases::CASES`]).
//!
//! # Dispatch
//!
//! Three face tests pack into a 3-bit `subcfg` ∈ `1..=8` (initialized
//! at 1, matching Julia's convention). The arms are *not* monotonically
//! mapped — the 7.2 sub-tilings live at `subcfg ∈ {2, 3, 5}` (skipping
//! 4), which goes to 7.3 instead. `subcfg == 8` triggers a per-edge
//! interior test that selects between 7.4.1 and 7.4.2.
//!
//! | `subcfg` | Subcase   | Triangles | Needs c-vertex? |
//! |----------|-----------|-----------|-----------------|
//! | `1`      | 7.1       | 3         | no              |
//! | `2`      | 7.2 (a)   | 5         | no              |
//! | `3`      | 7.2 (b)   | 5         | no              |
//! | `4`      | 7.3 (a)   | 9         | yes             |
//! | `5`      | 7.2 (c)   | 5         | no              |
//! | `6`      | 7.3 (b)   | 9         | yes             |
//! | `7`      | 7.3 (c)   | 9         | yes             |
//! | `8`      | 7.4.1 / .2 | 5 / 9    | sometimes       |
//!
//! # Sign / index conventions
//!
//! - **Tile edges** are stored 0-based 0..=11; edge index 12 (1-based 13)
//!   refers to the c-vertex computed by [`super::compute_c_vertex`].
//! - **`TEST_7` columns**:
//!   - 0..=2: face indices (1..=6 with sign-encoded inversion) for [`face_decider`].
//!   - 3: interior-test sentinel `s` for `subcfg == 8` (always ±7).
//!   - 4: per-edge reference index `0..=11` for [`interior_decider_per_edge`]
//!     (already 0-based in the source data).
//!
//! # Port
//!
//! Tables transcribed from `JuliaGeometry/MarchingCubes.jl` `src/lut.jl`:
//! - `test7` at lines 878-895
//! - `tiling7_1` at lines 907-924
//! - `tiling7_2` at lines 936-1002
//! - `tiling7_3` at lines 1014-1080
//! - `tiling7_4_1` at lines 1092-1109
//! - `tiling7_4_2` at lines 1121-1138
//!
//! Dispatch mirrors `case == 8` in `src/MarchingCubes.jl`, lines 234-254.

use super::super::{CellCoord, Mesh};
use super::decider::{face_decider, interior_decider_per_edge};
use super::emit_triangles;
use crate::mesh::grid::GridSpec;

/// Per-cfg 5-tuple: `[face_a, face_b, face_c, interior_s, ref_edge]`.
///
/// Source: Lewiner `test7[16][5]`. Columns 0..=2 are face indices for
/// [`face_decider`] (signed for inversion). Column 3 is the interior-
/// test `s` argument (±7). Column 4 is the reference edge for the
/// per-edge interior decider — **already 0-based** in the source.
#[rustfmt::skip]
const TEST_7: [[i8; 5]; 16] = [
    [  1,   2,   5,   7,   1],
    [  3,   4,   5,   7,   3],
    [  4,   1,   6,   7,   4],
    [  4,   1,   5,   7,   0],
    [  2,   3,   5,   7,   2],
    [  1,   2,   6,   7,   5],
    [  2,   3,   6,   7,   6],
    [  3,   4,   6,   7,   7],
    [ -3,  -4,  -6,  -7,   7],
    [ -2,  -3,  -6,  -7,   6],
    [ -1,  -2,  -6,  -7,   5],
    [ -2,  -3,  -5,  -7,   2],
    [ -4,  -1,  -5,  -7,   0],
    [ -4,  -1,  -6,  -7,   4],
    [ -3,  -4,  -5,  -7,   3],
    [ -1,  -2,  -5,  -7,   1],
];

/// Triangulation for Chernyaev subcase 7.1 (3 triangles, no c-vertex).
/// Indexed `[cfg]` for `subcfg == 1`.
///
/// Source: Lewiner `tiling7_1[16][9]`.
#[rustfmt::skip]
const TILING_7_1: [[u8; 9]; 16] = [
    [  9,   5,   4,  10,   1,   2,   8,   3,   0],
    [ 11,   7,   6,   8,   3,   0,  10,   1,   2],
    [  3,   0,   8,   5,   4,   9,   7,   6,  11],
    [  8,   4,   7,   9,   0,   1,  11,   2,   3],
    [ 10,   6,   5,  11,   2,   3,   9,   0,   1],
    [  0,   1,   9,   6,   5,  10,   4,   7,   8],
    [  1,   2,  10,   7,   6,  11,   5,   4,   9],
    [  2,   3,  11,   4,   7,   8,   6,   5,  10],
    [ 11,   3,   2,   8,   7,   4,  10,   5,   6],
    [ 10,   2,   1,  11,   6,   7,   9,   4,   5],
    [  9,   1,   0,  10,   5,   6,   8,   7,   4],
    [  5,   6,  10,   3,   2,  11,   1,   0,   9],
    [  7,   4,   8,   1,   0,   9,   3,   2,  11],
    [  8,   0,   3,   9,   4,   5,  11,   6,   7],
    [  6,   7,  11,   0,   3,   8,   2,   1,  10],
    [  4,   5,   9,   2,   1,  10,   0,   3,   8],
];

/// Triangulation for Chernyaev subcase 7.2 (5 triangles, no c-vertex).
/// Indexed `[cfg][sub]` where `sub ∈ {0, 1, 2}` per the dispatch table.
///
/// **Non-monotonic dispatch**: `subcfg=2 → sub=0`, `subcfg=3 → sub=1`,
/// `subcfg=5 → sub=2` (`subcfg=4` skips 7.2 entirely and goes to 7.3).
///
/// Source: Lewiner `tiling7_2[16][3][15]`.
#[rustfmt::skip]
const TILING_7_2: [[[u8; 15]; 3]; 16] = [
    [
        [  1,   2,  10,   3,   4,   8,   4,   3,   5,   0,   5,   3,   5,   0,   9],
        [  3,   0,   8,   9,   1,   4,   2,   4,   1,   4,   2,   5,  10,   5,   2],
        [  9,   5,   4,   0,  10,   1,  10,   0,   8,  10,   8,   2,   3,   2,   8],
    ],
    [
        [  3,   0,   8,   1,   6,  10,   6,   1,   7,   2,   7,   1,   7,   2,  11],
        [  1,   2,  10,  11,   3,   6,   0,   6,   3,   6,   0,   7,   8,   7,   0],
        [ 11,   7,   6,   2,   8,   3,   8,   2,  10,   8,  10,   0,   1,   0,  10],
    ],
    [
        [  9,   5,   4,  11,   3,   6,   0,   6,   3,   6,   0,   7,   8,   7,   0],
        [ 11,   7,   6,   3,   4,   8,   4,   3,   5,   0,   5,   3,   5,   0,   9],
        [  3,   0,   8,   4,   9,   7,  11,   7,   9,   5,  11,   9,  11,   5,   6],
    ],
    [
        [  0,   1,   9,   2,   7,  11,   7,   2,   4,   3,   4,   2,   4,   3,   8],
        [  2,   3,  11,   8,   0,   7,   1,   7,   0,   7,   1,   4,   9,   4,   1],
        [  8,   4,   7,   3,   9,   0,   9,   3,  11,   9,  11,   1,   2,   1,  11],
    ],
    [
        [  2,   3,  11,   0,   5,   9,   5,   0,   6,   1,   6,   0,   6,   1,  10],
        [  0,   1,   9,  10,   2,   5,   3,   5,   2,   5,   3,   6,  11,   6,   3],
        [  6,   5,  10,   1,  11,   2,  11,   1,   9,  11,   9,   3,   0,   3,   9],
    ],
    [
        [  6,   5,  10,   8,   0,   7,   1,   7,   0,   7,   1,   4,   9,   4,   1],
        [  8,   4,   7,   0,   5,   9,   5,   0,   6,   1,   6,   0,   6,   1,  10],
        [  0,   1,   9,   5,  10,   4,   8,   4,  10,   6,   8,  10,   8,   6,   7],
    ],
    [
        [ 11,   7,   6,   9,   1,   4,   2,   4,   1,   4,   2,   5,  10,   5,   2],
        [  9,   5,   4,   1,   6,  10,   6,   1,   7,   2,   7,   1,   7,   2,  11],
        [  1,   2,  10,   6,  11,   5,   9,   5,  11,   7,   9,  11,   9,   7,   4],
    ],
    [
        [  8,   4,   7,  10,   2,   5,   3,   5,   2,   5,   3,   6,  11,   6,   3],
        [  6,   5,  10,   2,   7,  11,   7,   2,   4,   3,   4,   2,   4,   3,   8],
        [  2,   3,  11,   7,   8,   6,  10,   6,   8,   4,  10,   8,  10,   4,   5],
    ],
    [
        [  7,   4,   8,   5,   2,  10,   2,   5,   3,   6,   3,   5,   3,   6,  11],
        [ 10,   5,   6,  11,   7,   2,   4,   2,   7,   2,   4,   3,   8,   3,   4],
        [ 11,   3,   2,   6,   8,   7,   8,   6,  10,   8,  10,   4,   5,   4,  10],
    ],
    [
        [  6,   7,  11,   4,   1,   9,   1,   4,   2,   5,   2,   4,   2,   5,  10],
        [  4,   5,   9,  10,   6,   1,   7,   1,   6,   1,   7,   2,  11,   2,   7],
        [ 10,   2,   1,   5,  11,   6,  11,   5,   9,  11,   9,   7,   4,   7,   9],
    ],
    [
        [ 10,   5,   6,   7,   0,   8,   0,   7,   1,   4,   1,   7,   1,   4,   9],
        [  7,   4,   8,   9,   5,   0,   6,   0,   5,   0,   6,   1,  10,   1,   6],
        [  9,   1,   0,   4,  10,   5,  10,   4,   8,  10,   8,   6,   7,   6,   8],
    ],
    [
        [ 11,   3,   2,   9,   5,   0,   6,   0,   5,   0,   6,   1,  10,   1,   6],
        [  9,   1,   0,   5,   2,  10,   2,   5,   3,   6,   3,   5,   3,   6,  11],
        [ 10,   5,   6,   2,  11,   1,   9,   1,  11,   3,   9,  11,   9,   3,   0],
    ],
    [
        [  9,   1,   0,  11,   7,   2,   4,   2,   7,   2,   4,   3,   8,   3,   4],
        [ 11,   3,   2,   7,   0,   8,   0,   7,   1,   4,   1,   7,   1,   4,   9],
        [  7,   4,   8,   0,   9,   3,  11,   3,   9,   1,  11,   9,  11,   1,   2],
    ],
    [
        [  4,   5,   9,   6,   3,  11,   3,   6,   0,   7,   0,   6,   0,   7,   8],
        [  6,   7,  11,   8,   4,   3,   5,   3,   4,   3,   5,   0,   9,   0,   5],
        [  8,   0,   3,   7,   9,   4,   9,   7,  11,   9,  11,   5,   6,   5,  11],
    ],
    [
        [  8,   0,   3,  10,   6,   1,   7,   1,   6,   1,   7,   2,  11,   2,   7],
        [ 10,   2,   1,   6,   3,  11,   3,   6,   0,   7,   0,   6,   0,   7,   8],
        [  6,   7,  11,   3,   8,   2,  10,   2,   8,   0,  10,   8,  10,   0,   1],
    ],
    [
        [ 10,   2,   1,   8,   4,   3,   5,   3,   4,   3,   5,   0,   9,   0,   5],
        [  8,   0,   3,   4,   1,   9,   1,   4,   2,   5,   2,   4,   2,   5,  10],
        [  4,   5,   9,   1,  10,   0,   8,   0,  10,   2,   8,  10,   8,   2,   3],
    ],
];

/// Triangulation for Chernyaev subcase 7.3 (9 triangles, uses c-vertex).
/// Indexed `[cfg][sub]` where `sub ∈ {0, 1, 2}` for `subcfg ∈ {4, 6, 7}`.
///
/// Source: Lewiner `tiling7_3[16][3][27]`.
#[rustfmt::skip]
const TILING_7_3: [[[u8; 27]; 3]; 16] = [
    [
        [ 12,   2,  10,  12,  10,   5,  12,   5,   4,  12,   4,   8,  12,   8,   3,  12,   3,   0,  12,   0,   9,  12,   9,   1,  12,   1,   2],
        [ 12,   5,   4,  12,   4,   8,  12,   8,   3,  12,   3,   2,  12,   2,  10,  12,  10,   1,  12,   1,   0,  12,   0,   9,  12,   9,   5],
        [  5,   4,  12,  10,   5,  12,   2,  10,  12,   3,   2,  12,   8,   3,  12,   0,   8,  12,   1,   0,  12,   9,   1,  12,   4,   9,  12],
    ],
    [
        [ 12,   0,   8,  12,   8,   7,  12,   7,   6,  12,   6,  10,  12,  10,   1,  12,   1,   2,  12,   2,  11,  12,  11,   3,  12,   3,   0],
        [ 12,   7,   6,  12,   6,  10,  12,  10,   1,  12,   1,   0,  12,   0,   8,  12,   8,   3,  12,   3,   2,  12,   2,  11,  12,  11,   7],
        [  7,   6,  12,   8,   7,  12,   0,   8,  12,   1,   0,  12,  10,   1,  12,   2,  10,  12,   3,   2,  12,  11,   3,  12,   6,  11,  12],
    ],
    [
        [  9,   5,  12,   0,   9,  12,   3,   0,  12,  11,   3,  12,   6,  11,  12,   7,   6,  12,   8,   7,  12,   4,   8,  12,   5,   4,  12],
        [  3,   0,  12,  11,   3,  12,   6,  11,  12,   5,   6,  12,   9,   5,  12,   4,   9,  12,   7,   4,  12,   8,   7,  12,   0,   8,  12],
        [ 12,   3,   0,  12,   0,   9,  12,   9,   5,  12,   5,   6,  12,   6,  11,  12,  11,   7,  12,   7,   4,  12,   4,   8,  12,   8,   3],
    ],
    [
        [ 12,   1,   9,  12,   9,   4,  12,   4,   7,  12,   7,  11,  12,  11,   2,  12,   2,   3,  12,   3,   8,  12,   8,   0,  12,   0,   1],
        [ 12,   4,   7,  12,   7,  11,  12,  11,   2,  12,   2,   1,  12,   1,   9,  12,   9,   0,  12,   0,   3,  12,   3,   8,  12,   8,   4],
        [  4,   7,  12,   9,   4,  12,   1,   9,  12,   2,   1,  12,  11,   2,  12,   3,  11,  12,   0,   3,  12,   8,   0,  12,   7,   8,  12],
    ],
    [
        [ 12,   3,  11,  12,  11,   6,  12,   6,   5,  12,   5,   9,  12,   9,   0,  12,   0,   1,  12,   1,  10,  12,  10,   2,  12,   2,   3],
        [ 12,   6,   5,  12,   5,   9,  12,   9,   0,  12,   0,   3,  12,   3,  11,  12,  11,   2,  12,   2,   1,  12,   1,  10,  12,  10,   6],
        [  6,   5,  12,  11,   6,  12,   3,  11,  12,   0,   3,  12,   9,   0,  12,   1,   9,  12,   2,   1,  12,  10,   2,  12,   5,  10,  12],
    ],
    [
        [ 10,   6,  12,   1,  10,  12,   0,   1,  12,   8,   0,  12,   7,   8,  12,   4,   7,  12,   9,   4,  12,   5,   9,  12,   6,   5,  12],
        [  0,   1,  12,   8,   0,  12,   7,   8,  12,   6,   7,  12,  10,   6,  12,   5,  10,  12,   4,   5,  12,   9,   4,  12,   1,   9,  12],
        [ 12,   0,   1,  12,   1,  10,  12,  10,   6,  12,   6,   7,  12,   7,   8,  12,   8,   4,  12,   4,   5,  12,   5,   9,  12,   9,   0],
    ],
    [
        [ 11,   7,  12,   2,  11,  12,   1,   2,  12,   9,   1,  12,   4,   9,  12,   5,   4,  12,  10,   5,  12,   6,  10,  12,   7,   6,  12],
        [  1,   2,  12,   9,   1,  12,   4,   9,  12,   7,   4,  12,  11,   7,  12,   6,  11,  12,   5,   6,  12,  10,   5,  12,   2,  10,  12],
        [ 12,   1,   2,  12,   2,  11,  12,  11,   7,  12,   7,   4,  12,   4,   9,  12,   9,   5,  12,   5,   6,  12,   6,  10,  12,  10,   1],
    ],
    [
        [  8,   4,  12,   3,   8,  12,   2,   3,  12,  10,   2,  12,   5,  10,  12,   6,   5,  12,  11,   6,  12,   7,  11,  12,   4,   7,  12],
        [  2,   3,  12,  10,   2,  12,   5,  10,  12,   4,   5,  12,   8,   4,  12,   7,   8,  12,   6,   7,  12,  11,   6,  12,   3,  11,  12],
        [ 12,   2,   3,  12,   3,   8,  12,   8,   4,  12,   4,   5,  12,   5,  10,  12,  10,   6,  12,   6,   7,  12,   7,  11,  12,  11,   2],
    ],
    [
        [ 12,   4,   8,  12,   8,   3,  12,   3,   2,  12,   2,  10,  12,  10,   5,  12,   5,   6,  12,   6,  11,  12,  11,   7,  12,   7,   4],
        [ 12,   3,   2,  12,   2,  10,  12,  10,   5,  12,   5,   4,  12,   4,   8,  12,   8,   7,  12,   7,   6,  12,   6,  11,  12,  11,   3],
        [  3,   2,  12,   8,   3,  12,   4,   8,  12,   5,   4,  12,  10,   5,  12,   6,  10,  12,   7,   6,  12,  11,   7,  12,   2,  11,  12],
    ],
    [
        [ 12,   7,  11,  12,  11,   2,  12,   2,   1,  12,   1,   9,  12,   9,   4,  12,   4,   5,  12,   5,  10,  12,  10,   6,  12,   6,   7],
        [ 12,   2,   1,  12,   1,   9,  12,   9,   4,  12,   4,   7,  12,   7,  11,  12,  11,   6,  12,   6,   5,  12,   5,  10,  12,  10,   2],
        [  2,   1,  12,  11,   2,  12,   7,  11,  12,   4,   7,  12,   9,   4,  12,   5,   9,  12,   6,   5,  12,  10,   6,  12,   1,  10,  12],
    ],
    [
        [ 12,   6,  10,  12,  10,   1,  12,   1,   0,  12,   0,   8,  12,   8,   7,  12,   7,   4,  12,   4,   9,  12,   9,   5,  12,   5,   6],
        [ 12,   1,   0,  12,   0,   8,  12,   8,   7,  12,   7,   6,  12,   6,  10,  12,  10,   5,  12,   5,   4,  12,   4,   9,  12,   9,   1],
        [  1,   0,  12,  10,   1,  12,   6,  10,  12,   7,   6,  12,   8,   7,  12,   4,   8,  12,   5,   4,  12,   9,   5,  12,   0,   9,  12],
    ],
    [
        [ 11,   3,  12,   6,  11,  12,   5,   6,  12,   9,   5,  12,   0,   9,  12,   1,   0,  12,  10,   1,  12,   2,  10,  12,   3,   2,  12],
        [  5,   6,  12,   9,   5,  12,   0,   9,  12,   3,   0,  12,  11,   3,  12,   2,  11,  12,   1,   2,  12,  10,   1,  12,   6,  10,  12],
        [ 12,   5,   6,  12,   6,  11,  12,  11,   3,  12,   3,   0,  12,   0,   9,  12,   9,   1,  12,   1,   2,  12,   2,  10,  12,  10,   5],
    ],
    [
        [  9,   1,  12,   4,   9,  12,   7,   4,  12,  11,   7,  12,   2,  11,  12,   3,   2,  12,   8,   3,  12,   0,   8,  12,   1,   0,  12],
        [  7,   4,  12,  11,   7,  12,   2,  11,  12,   1,   2,  12,   9,   1,  12,   0,   9,  12,   3,   0,  12,   8,   3,  12,   4,   8,  12],
        [ 12,   7,   4,  12,   4,   9,  12,   9,   1,  12,   1,   2,  12,   2,  11,  12,  11,   3,  12,   3,   0,  12,   0,   8,  12,   8,   7],
    ],
    [
        [ 12,   5,   9,  12,   9,   0,  12,   0,   3,  12,   3,  11,  12,  11,   6,  12,   6,   7,  12,   7,   8,  12,   8,   4,  12,   4,   5],
        [ 12,   0,   3,  12,   3,  11,  12,  11,   6,  12,   6,   5,  12,   5,   9,  12,   9,   4,  12,   4,   7,  12,   7,   8,  12,   8,   0],
        [  0,   3,  12,   9,   0,  12,   5,   9,  12,   6,   5,  12,  11,   6,  12,   7,  11,  12,   4,   7,  12,   8,   4,  12,   3,   8,  12],
    ],
    [
        [  8,   0,  12,   7,   8,  12,   6,   7,  12,  10,   6,  12,   1,  10,  12,   2,   1,  12,  11,   2,  12,   3,  11,  12,   0,   3,  12],
        [  6,   7,  12,  10,   6,  12,   1,  10,  12,   0,   1,  12,   8,   0,  12,   3,   8,  12,   2,   3,  12,  11,   2,  12,   7,  11,  12],
        [ 12,   6,   7,  12,   7,   8,  12,   8,   0,  12,   0,   1,  12,   1,  10,  12,  10,   2,  12,   2,   3,  12,   3,  11,  12,  11,   6],
    ],
    [
        [ 10,   2,  12,   5,  10,  12,   4,   5,  12,   8,   4,  12,   3,   8,  12,   0,   3,  12,   9,   0,  12,   1,   9,  12,   2,   1,  12],
        [  4,   5,  12,   8,   4,  12,   3,   8,  12,   2,   3,  12,  10,   2,  12,   1,  10,  12,   0,   1,  12,   9,   0,  12,   5,   9,  12],
        [ 12,   4,   5,  12,   5,  10,  12,  10,   2,  12,   2,   3,  12,   3,   8,  12,   8,   0,  12,   0,   1,  12,   1,   9,  12,   9,   4],
    ],
];

/// Triangulation for Chernyaev subcase 7.4.1 (5 triangles, no c-vertex).
/// Selected when [`interior_decider_per_edge`] returns `true` at `subcfg == 8`.
///
/// Source: Lewiner `tiling7_4_1[16][15]`.
#[rustfmt::skip]
const TILING_7_4_1: [[u8; 15]; 16] = [
    [  3,   4,   8,   4,   3,  10,   2,  10,   3,   4,  10,   5,   9,   1,   0],
    [  1,   6,  10,   6,   1,   8,   0,   8,   1,   6,   8,   7,  11,   3,   2],
    [ 11,   3,   6,   9,   6,   3,   6,   9,   5,   0,   9,   3,   7,   4,   8],
    [  2,   7,  11,   7,   2,   9,   1,   9,   2,   7,   9,   4,   8,   0,   3],
    [  0,   5,   9,   5,   0,  11,   3,  11,   0,   5,  11,   6,  10,   2,   1],
    [  8,   0,   7,  10,   7,   0,   7,  10,   6,   1,  10,   0,   4,   5,   9],
    [  9,   1,   4,  11,   4,   1,   4,  11,   7,   2,  11,   1,   5,   6,  10],
    [ 10,   2,   5,   8,   5,   2,   5,   8,   4,   3,   8,   2,   6,   7,  11],
    [  5,   2,  10,   2,   5,   8,   4,   8,   5,   2,   8,   3,  11,   7,   6],
    [  4,   1,   9,   1,   4,  11,   7,  11,   4,   1,  11,   2,  10,   6,   5],
    [  7,   0,   8,   0,   7,  10,   6,  10,   7,   0,  10,   1,   9,   5,   4],
    [  9,   5,   0,  11,   0,   5,   0,  11,   3,   6,  11,   5,   1,   2,  10],
    [ 11,   7,   2,   9,   2,   7,   2,   9,   1,   4,   9,   7,   3,   0,   8],
    [  6,   3,  11,   3,   6,   9,   5,   9,   6,   3,   9,   0,   8,   4,   7],
    [ 10,   6,   1,   8,   1,   6,   1,   8,   0,   7,   8,   6,   2,   3,  11],
    [  8,   4,   3,  10,   3,   4,   3,  10,   2,   5,  10,   4,   0,   1,   9],
];

/// Triangulation for Chernyaev subcase 7.4.2 (9 triangles, uses c-vertex).
/// Selected when [`interior_decider_per_edge`] returns `false` at `subcfg == 8`.
///
/// Source: Lewiner `tiling7_4_2[16][27]`.
#[rustfmt::skip]
const TILING_7_4_2: [[u8; 27]; 16] = [
    [  9,   4,   8,   4,   9,   5,  10,   5,   9,   1,  10,   9,  10,   1,   2,   0,   2,   1,   2,   0,   3,   8,   3,   0,   9,   8,   0],
    [ 11,   6,  10,   6,  11,   7,   8,   7,  11,   3,   8,  11,   8,   3,   0,   2,   0,   3,   0,   2,   1,  10,   1,   2,  11,  10,   2],
    [ 11,   3,   8,   0,   8,   3,   8,   0,   9,   8,   9,   4,   5,   4,   9,   4,   5,   7,   6,   7,   5,   7,   6,  11,   7,  11,   8],
    [  8,   7,  11,   7,   8,   4,   9,   4,   8,   0,   9,   8,   9,   0,   1,   3,   1,   0,   1,   3,   2,  11,   2,   3,   8,  11,   3],
    [ 10,   5,   9,   5,  10,   6,  11,   6,  10,   2,  11,  10,  11,   2,   3,   1,   3,   2,   3,   1,   0,   9,   0,   1,  10,   9,   1],
    [  8,   0,   9,   1,   9,   0,   9,   1,  10,   9,  10,   5,   6,   5,  10,   5,   6,   4,   7,   4,   6,   4,   7,   8,   4,   8,   9],
    [  9,   1,  10,   2,  10,   1,  10,   2,  11,  10,  11,   6,   7,   6,  11,   6,   7,   5,   4,   5,   7,   5,   4,   9,   5,   9,  10],
    [ 10,   2,  11,   3,  11,   2,  11,   3,   8,  11,   8,   7,   4,   7,   8,   7,   4,   6,   5,   6,   4,   6,   5,  10,   6,  10,  11],
    [ 11,   2,  10,   2,  11,   3,   8,   3,  11,   7,   8,  11,   8,   7,   4,   6,   4,   7,   4,   6,   5,  10,   5,   6,  11,  10,   6],
    [ 10,   1,   9,   1,  10,   2,  11,   2,  10,   6,  11,  10,  11,   6,   7,   5,   7,   6,   7,   5,   4,   9,   4,   5,  10,   9,   5],
    [  9,   0,   8,   0,   9,   1,  10,   1,   9,   5,  10,   9,  10,   5,   6,   4,   6,   5,   6,   4,   7,   8,   7,   4,   9,   8,   4],
    [  9,   5,  10,   6,  10,   5,  10,   6,  11,  10,  11,   2,   3,   2,  11,   2,   3,   1,   0,   1,   3,   1,   0,   9,   1,   9,  10],
    [ 11,   7,   8,   4,   8,   7,   8,   4,   9,   8,   9,   0,   1,   0,   9,   0,   1,   3,   2,   3,   1,   3,   2,  11,   3,  11,   8],
    [  8,   3,  11,   3,   8,   0,   9,   0,   8,   4,   9,   8,   9,   4,   5,   7,   5,   4,   5,   7,   6,  11,   6,   7,   8,  11,   7],
    [ 10,   6,  11,   7,  11,   6,  11,   7,   8,  11,   8,   3,   0,   3,   8,   3,   0,   2,   1,   2,   0,   2,   1,  10,   2,  10,  11],
    [  8,   4,   9,   5,   9,   4,   9,   5,  10,   9,  10,   1,   2,   1,  10,   1,   2,   0,   3,   0,   2,   0,   3,   8,   0,   8,   9],
];

/// Emits Case-7 triangles for the given voxel.
///
/// `subcase` is Lewiner's 1-based subcase index in `1..=16` from
/// [`super::cases::CASES`].
pub(super) fn emit(
    grid: &GridSpec,
    cell: CellCoord,
    corners: &[f32; 8],
    subcase: i8,
    mesh: &mut Mesh,
) {
    debug_assert!(
        (1..=16).contains(&subcase),
        "Case 7 subcase out of range: {subcase}"
    );
    let cfg = usize::try_from(subcase - 1).unwrap_or(0);

    // Pack 3 face-decider results into bits 0..2 of `subcfg`,
    // initialized at 1 (matching Julia). Final value ∈ 1..=8.
    let mut subcfg: u8 = 1;
    for (bit, &face) in TEST_7[cfg].iter().take(3).enumerate() {
        if face_decider(corners, face) {
            subcfg += 1 << bit;
        }
    }

    let tiling: &[u8] = match subcfg {
        1 => &TILING_7_1[cfg],
        2 => &TILING_7_2[cfg][0],
        3 => &TILING_7_2[cfg][1],
        4 => &TILING_7_3[cfg][0],
        // Non-monotonic: subcfg=5 reuses 7.2 sub-row 2 (the third sub-tiling),
        // not 7.3 sub-row 1. Lewiner's published dispatch.
        5 => &TILING_7_2[cfg][2],
        6 => &TILING_7_3[cfg][1],
        7 => &TILING_7_3[cfg][2],
        8 => {
            // Interior test: s = column 3 (±7), edge = column 4 (0..=7,
            // already 0-based).
            let s = TEST_7[cfg][3];
            let edge = TEST_7[cfg][4];
            let edge_u8 = u8::try_from(edge).unwrap_or(0);
            if interior_decider_per_edge(corners, edge_u8, s) {
                &TILING_7_4_2[cfg][..]
            } else {
                &TILING_7_4_1[cfg][..]
            }
        }
        _ => unreachable!("subcfg packed from 3 bits + 1, must be 1..=8"),
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

    fn case8_idx_subcase_pairs() -> Vec<(u8, i8)> {
        let mut hits = Vec::new();
        for our_idx in 0..=255_u8 {
            let lewiner_idx = !our_idx;
            let (base_case, subcase) = super::super::cases::CASES[lewiner_idx as usize];
            if base_case == 8 {
                hits.push((our_idx, subcase));
            }
        }
        hits
    }

    /// Pin: exactly 16 voxels map to Lewiner base case 8.
    #[test]
    fn case8_dispatch_has_sixteen_configs() {
        assert_eq!(case8_idx_subcase_pairs().len(), 16);
    }

    /// Transcription oracle for `test7[1]` (Julia line 879 → 0-based,
    /// columns 0..2 unchanged faces, column 3 unchanged sentinel,
    /// column 4 already 0-based in source).
    #[test]
    fn test_7_first_row_matches_lut() {
        assert_eq!(TEST_7[0], [1, 2, 5, 7, 1]);
    }

    /// Transcription oracle for `tiling7_1[1]` (Julia line 908 →
    /// `Int8.((10, 6, 5, 11, 2, 3, 9, 4, 1))` → 0-based).
    #[test]
    fn tiling_7_1_first_row_matches_lut() {
        assert_eq!(TILING_7_1[0], [9, 5, 4, 10, 1, 2, 8, 3, 0]);
    }

    /// Transcription oracle for `tiling7_2[1][1]` (Julia line 938).
    #[test]
    fn tiling_7_2_first_row_matches_lut() {
        let expected = [1_u8, 2, 10, 3, 4, 8, 4, 3, 5, 0, 5, 3, 5, 0, 9];
        assert_eq!(TILING_7_2[0][0], expected);
    }

    /// Transcription oracle for `tiling7_3[1][1]` (Julia line 1016).
    #[test]
    fn tiling_7_3_first_row_matches_lut() {
        let expected = [
            12_u8, 2, 10, 12, 10, 5, 12, 5, 4, 12, 4, 8, 12, 8, 3, 12, 3, 0, 12, 0, 9, 12, 9, 1,
            12, 1, 2,
        ];
        assert_eq!(TILING_7_3[0][0], expected);
    }

    /// Transcription oracle for `tiling7_4_1[1]` (Julia line 1093).
    #[test]
    fn tiling_7_4_1_first_row_matches_lut() {
        let expected = [3_u8, 4, 8, 4, 3, 10, 2, 10, 3, 4, 10, 5, 9, 1, 0];
        assert_eq!(TILING_7_4_1[0], expected);
    }

    /// Transcription oracle for `tiling7_4_2[1]` (Julia line 1122).
    #[test]
    fn tiling_7_4_2_first_row_matches_lut() {
        let expected = [
            9_u8, 4, 8, 4, 9, 5, 10, 5, 9, 1, 10, 9, 10, 1, 2, 0, 2, 1, 2, 0, 3, 8, 3, 0, 9, 8, 0,
        ];
        assert_eq!(TILING_7_4_2[0], expected);
    }

    /// Sharp oracle: every entry in every tiling is `≤ 12` (cube edges
    /// 0..=11 plus c-vertex slot 12).
    #[test]
    fn tilings_use_only_legal_indices() {
        for cfg in 0..16 {
            for &e in &TILING_7_1[cfg] {
                assert!(e <= 12, "TILING_7_1[{cfg}] index {e}");
            }
            for sub in 0..3 {
                for &e in &TILING_7_2[cfg][sub] {
                    assert!(e <= 12, "TILING_7_2[{cfg}][{sub}] index {e}");
                }
                for &e in &TILING_7_3[cfg][sub] {
                    assert!(e <= 12, "TILING_7_3[{cfg}][{sub}] index {e}");
                }
            }
            for &e in &TILING_7_4_1[cfg] {
                assert!(e <= 12, "TILING_7_4_1[{cfg}] index {e}");
            }
            for &e in &TILING_7_4_2[cfg] {
                assert!(e <= 12, "TILING_7_4_2[{cfg}] index {e}");
            }
        }
    }

    /// Row-length pin.
    #[test]
    fn tiling_dimensions_match_lewiner() {
        assert_eq!(TEST_7.len(), 16);
        assert_eq!(TEST_7[0].len(), 5);
        assert_eq!(TILING_7_1.len(), 16);
        assert_eq!(TILING_7_1[0].len(), 9);
        assert_eq!(TILING_7_2.len(), 16);
        assert_eq!(TILING_7_2[0].len(), 3);
        assert_eq!(TILING_7_2[0][0].len(), 15);
        assert_eq!(TILING_7_3.len(), 16);
        assert_eq!(TILING_7_3[0].len(), 3);
        assert_eq!(TILING_7_3[0][0].len(), 27);
        assert_eq!(TILING_7_4_1.len(), 16);
        assert_eq!(TILING_7_4_1[0].len(), 15);
        assert_eq!(TILING_7_4_2.len(), 16);
        assert_eq!(TILING_7_4_2[0].len(), 27);
    }

    /// Behavioral: every Lewiner-8 voxel under ±1 corner inputs runs
    /// `emit` to completion, produces a non-empty mesh, and the mesh
    /// satisfies the unwelded invariant.
    #[test]
    fn case7_dispatch_runs_on_every_case8_voxel() {
        for (our_idx, subcase) in case8_idx_subcase_pairs() {
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

    /// Edge: explicitly drive `subcfg == 5` and verify it routes to
    /// `TILING_7_2[cfg][2]` (the third sub-row, *not* the second).
    /// Catches a naive `subcfg - 2` indexing bug.
    #[test]
    fn case7_subcfg_5_uses_third_subrow_of_7_2() {
        // subcfg=5 means face tests {1: yes, 2: no, 3: yes} → bits 0 and 2 set.
        // For cfg 0, TEST_7[0] = [1, 2, 5, 7, 1]. We need face_decider to
        // return true for face 1 and 5, false for face 2.
        //
        // Constructing such a voxel directly is fiddly because face_decider
        // depends on the four corner values per face. Instead, we sweep
        // perturbations and check at least one lands in subcfg=5 — then
        // assert it produces 5 triangles AND the first triangle's first
        // edge index matches the first edge of TILING_7_2[cfg][2].
        let mut found = false;
        for (our_idx, subcase) in case8_idx_subcase_pairs() {
            for perturb_seed in 0..256_u32 {
                let mut corners = [0.0_f32; 8];
                for (bit, c) in corners.iter_mut().enumerate() {
                    let inside = our_idx & (1 << bit) != 0;
                    let mag = 0.2
                        + 1.6
                            * f32::from(
                                ((perturb_seed.wrapping_mul(13) ^ u32::try_from(bit).unwrap_or(0))
                                    & 0xF) as u8,
                            )
                            / 15.0;
                    *c = if inside { -mag } else { mag };
                }
                let cfg = usize::try_from(subcase - 1).unwrap_or(0);
                let mut subcfg: u8 = 1;
                for (bit, &face) in TEST_7[cfg].iter().take(3).enumerate() {
                    if face_decider(&corners, face) {
                        subcfg += 1 << bit;
                    }
                }
                if subcfg != 5 {
                    continue;
                }
                let (grid, _field) = unit_grid_with_corners(&corners);
                let mut mesh = Mesh::default();
                let cell = CellCoord { x: 0, y: 0, z: 0 };
                emit(&grid, cell, &corners, subcase, &mut mesh);
                assert_eq!(
                    mesh.triangle_count(),
                    5,
                    "subcfg=5 must produce 5 tris (7.2 sub-row 2)"
                );
                // Verify the first triangle's first vertex came from the
                // first edge of TILING_7_2[cfg][2] (sub-row 2). We can
                // do this by re-running with the explicit tiling and
                // checking vertex equality on the first vertex.
                let expected_edge = TILING_7_2[cfg][2][0];
                let (a, b) = super::super::super::tables::EDGE_CORNERS[expected_edge as usize];
                let expected_v0 = super::super::super::interpolate_edge(
                    &grid,
                    cell,
                    a,
                    b,
                    corners[a as usize],
                    corners[b as usize],
                );
                assert_eq!(
                    mesh.vertices[0], expected_v0,
                    "subcfg=5 first vertex must come from TILING_7_2[cfg][2][0] = edge {expected_edge}"
                );
                found = true;
                break;
            }
            if found {
                break;
            }
        }
        assert!(
            found,
            "perturbation sweep failed to hit subcfg=5; widen the sweep"
        );
    }

    /// Adversarial: every triangle is geometrically non-degenerate.
    #[test]
    fn case7_triangles_are_nondegenerate() {
        for (our_idx, subcase) in case8_idx_subcase_pairs() {
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

    /// Pin: case 8 is no longer in the warn-once ambiguous set.
    #[test]
    fn case8_is_no_longer_ambiguous() {
        // Indirect verification: emit returns and produces non-empty mesh.
        for (our_idx, subcase) in case8_idx_subcase_pairs() {
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
                "case 8 voxel our_idx={our_idx:#010b} must produce triangles"
            );
        }
    }
}
