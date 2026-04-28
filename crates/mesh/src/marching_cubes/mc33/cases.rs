//! Lewiner 2003 `cases[256]` base-case lookup table.
//!
//! Each entry `CASES[i] = (base_case, subcase)` for Lewiner's `case_index` `i`
//! (where bit `p` set iff corner `p+1` is *outside* — the opposite of our
//! internal convention; see [`super::lewiner_index`] for the conversion).
//!
//! `base_case` values run `1..=15` using Lewiner's numbering, which is
//! Chernyaev's case number plus one:
//!
//! | Lewiner case | Chernyaev case | Description                | Ambiguous? |
//! |--------------|----------------|----------------------------|------------|
//! | 1            | 0              | Empty (no surface)         | —          |
//! | 2            | 1              | Single corner              | No         |
//! | 3            | 2              | Two adjacent corners       | No         |
//! | 4            | 3              | Two diag corners on a face | **Yes (face)** |
//! | 5            | 4              | Two diag corners on body   | **Yes (body)** |
//! | 6            | 5              | Three adjacent corners     | No         |
//! | 7            | 6              | Three corners + diag       | **Yes (face+body)** |
//! | 8            | 7              | Three non-adjacent         | **Yes (face+body)** |
//! | 9            | 8              | Four corners on a face     | No         |
//! | 10           | 9              | Four corners, "zigzag"     | No         |
//! | 11           | 10             | Two opposite diagonal pairs| **Yes** |
//! | 12           | 11             | Cross-face pattern         | No         |
//! | 13           | 12             | 3+1 + diag                 | **Yes** |
//! | 14           | 13             | Tunnel configuration       | **Yes (body, 7 subcases)** |
//! | 15           | 14             | Four "flower" corners      | No         |
//!
//! `subcase` is 1-based within each base case (matching Lewiner's / Julia's
//! convention); dispatchers subtract 1 when indexing into the tiling
//! tables.
//!
//! Source: Lewiner, Lopes, Vieira & Tavares, "Efficient Implementation of
//! Marching Cubes' Cases with Topological Guarantees," JGT 8(2), 2003,
//! `LookUpTable.h`, file dated 12/08/2002. Transcribed via `JuliaGeometry`'s
//! `MarchingCubes.jl` (`src/lut.jl`).

/// `CASES[lewiner_index] = (base_case, subcase)`. Both values are 1-based.
pub(crate) const CASES: [(i8, i8); 256] = [
    (1, 0),
    (2, 1),
    (2, 2),
    (3, 1),
    (2, 3),
    (4, 1),
    (3, 4),
    (6, 1),
    (2, 4),
    (3, 2),
    (4, 4),
    (6, 2),
    (3, 6),
    (6, 5),
    (6, 10),
    (9, 1),
    (2, 5),
    (3, 3),
    (4, 5),
    (6, 3),
    (5, 3),
    (7, 3),
    (7, 10),
    (12, 1),
    (4, 9),
    (6, 6),
    (8, 4),
    (10, 2),
    (7, 17),
    (15, 4),
    (13, 13),
    (6, 25),
    (2, 6),
    (4, 2),
    (3, 5),
    (6, 4),
    (4, 7),
    (8, 1),
    (6, 11),
    (10, 1),
    (5, 4),
    (7, 5),
    (7, 12),
    (15, 2),
    (7, 18),
    (13, 5),
    (12, 7),
    (6, 26),
    (3, 9),
    (6, 8),
    (6, 13),
    (9, 2),
    (7, 19),
    (13, 6),
    (15, 8),
    (6, 29),
    (7, 22),
    (12, 5),
    (13, 16),
    (6, 31),
    (11, 6),
    (7, 33),
    (7, 40),
    (3, 13),
    (2, 7),
    (5, 1),
    (4, 6),
    (7, 1),
    (3, 7),
    (7, 4),
    (6, 12),
    (15, 1),
    (4, 10),
    (7, 6),
    (8, 5),
    (13, 2),
    (6, 15),
    (12, 4),
    (10, 5),
    (6, 27),
    (4, 11),
    (7, 7),
    (8, 6),
    (13, 3),
    (7, 20),
    (11, 2),
    (13, 14),
    (7, 25),
    (8, 8),
    (13, 10),
    (14, 2),
    (8, 10),
    (13, 21),
    (7, 34),
    (8, 14),
    (4, 13),
    (3, 11),
    (7, 8),
    (6, 14),
    (12, 3),
    (6, 17),
    (13, 8),
    (9, 4),
    (6, 30),
    (7, 23),
    (11, 3),
    (13, 18),
    (7, 28),
    (15, 10),
    (7, 35),
    (6, 40),
    (3, 15),
    (6, 21),
    (15, 6),
    (10, 6),
    (6, 33),
    (12, 11),
    (7, 36),
    (6, 42),
    (3, 17),
    (13, 24),
    (7, 38),
    (8, 15),
    (4, 17),
    (7, 47),
    (5, 7),
    (4, 22),
    (2, 9),
    (2, 8),
    (4, 3),
    (5, 2),
    (7, 2),
    (4, 8),
    (8, 2),
    (7, 11),
    (13, 1),
    (3, 8),
    (6, 7),
    (7, 13),
    (12, 2),
    (6, 16),
    (10, 3),
    (15, 7),
    (6, 28),
    (3, 10),
    (6, 9),
    (7, 14),
    (15, 3),
    (7, 21),
    (13, 7),
    (11, 4),
    (7, 26),
    (6, 19),
    (9, 3),
    (13, 17),
    (6, 32),
    (12, 10),
    (6, 35),
    (7, 41),
    (3, 14),
    (4, 12),
    (8, 3),
    (7, 15),
    (13, 4),
    (8, 7),
    (14, 1),
    (13, 15),
    (8, 9),
    (7, 24),
    (13, 11),
    (11, 5),
    (7, 29),
    (13, 22),
    (8, 11),
    (7, 42),
    (4, 14),
    (6, 22),
    (10, 4),
    (12, 9),
    (6, 34),
    (13, 23),
    (8, 12),
    (7, 43),
    (4, 15),
    (15, 12),
    (6, 37),
    (7, 45),
    (3, 18),
    (7, 48),
    (4, 19),
    (5, 8),
    (2, 10),
    (3, 12),
    (7, 9),
    (7, 16),
    (11, 1),
    (6, 18),
    (13, 9),
    (12, 8),
    (7, 27),
    (6, 20),
    (15, 5),
    (13, 19),
    (7, 30),
    (9, 5),
    (6, 36),
    (6, 41),
    (3, 16),
    (6, 23),
    (12, 6),
    (13, 20),
    (7, 31),
    (15, 11),
    (7, 37),
    (7, 44),
    (5, 5),
    (10, 8),
    (6, 38),
    (8, 16),
    (4, 18),
    (6, 45),
    (3, 20),
    (4, 23),
    (2, 11),
    (6, 24),
    (13, 12),
    (15, 9),
    (7, 32),
    (10, 7),
    (8, 13),
    (6, 43),
    (4, 16),
    (12, 12),
    (7, 39),
    (7, 46),
    (5, 6),
    (6, 46),
    (4, 20),
    (3, 22),
    (2, 12),
    (9, 6),
    (6, 39),
    (6, 44),
    (3, 19),
    (6, 47),
    (4, 21),
    (3, 23),
    (2, 13),
    (6, 48),
    (3, 21),
    (4, 24),
    (2, 14),
    (3, 24),
    (2, 15),
    (2, 16),
    (1, 0),
];

#[cfg(test)]
#[allow(clippy::unwrap_used, missing_docs)]
mod tests {
    use super::*;

    #[test]
    fn empty_cases_at_both_ends() {
        // Both all-outside (index 0) and all-inside (index 255) are Lewiner
        // "case 1" (empty surface).
        assert_eq!(CASES[0], (1, 0));
        assert_eq!(CASES[255], (1, 0));
    }

    #[test]
    fn single_corner_cases_are_base_case_2() {
        // Any case_index with exactly one bit set is Lewiner's case 2 (one corner).
        for bit in 0..8 {
            let idx = 1_usize << bit;
            assert_eq!(
                CASES[idx].0, 2,
                "case_index {idx:#010b} (single bit {bit}) should be base case 2; got {:?}",
                CASES[idx]
            );
        }
    }

    #[test]
    fn face_ambiguous_diagonal_case_is_base_case_4() {
        // case_index 5 = 0b00000101: Lewiner bits 0 and 2 set = Lewiner
        // corners 1 and 3 outside. These corners sit at (0,0,0) and (1,1,0)
        // in the cube — diagonally opposite on the bottom face (z=0). This is
        // the canonical Chernyaev case 3 = Lewiner case 4 (face-ambiguous).
        assert_eq!(CASES[5].0, 4);
    }

    #[test]
    fn complement_symmetry() {
        // For any case, the all-inside complement (bitwise NOT) must map to
        // the same base case (geometry is sign-symmetric). Verify on a few.
        assert_eq!(CASES[0].0, CASES[255].0);
        assert_eq!(CASES[5].0, CASES[!5_u8 as usize].0);
        assert_eq!(CASES[7].0, CASES[!7_u8 as usize].0);
    }
}
