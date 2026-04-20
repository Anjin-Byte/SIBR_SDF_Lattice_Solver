//! Marching Cubes 33 (Chernyaev 1995 / Lewiner 2003) — per-voxel processing.
//!
//! # Current scope (incrementally being ported)
//!
//! - **Case 3** (face-ambiguous, 2 diag corners on a face): **fully
//!   disambiguated** via [`decider::face_decider`] and
//!   [`case3`]'s alternative triangulation tables.
//! - **Case 4** (body-ambiguous, 2 diag corners on a body diagonal):
//!   **fully disambiguated** via [`decider::interior_decider`] and
//!   [`case4`]'s alternative triangulation tables.
//! - **Unambiguous cases (2, 3, 6, 9, 10, 12, 15 in Lewiner's
//!   numbering)**: handled by [`unambiguous`]'s face-consistent
//!   Lewiner tables. Replaces classic MC's 1987 triangulations, which
//!   are self-manifold per voxel but not consistent across shared faces.
//!   This is what eliminates the ~2–3% non-manifold edge rate on typical
//!   SDF workloads.
//! - **Empty case (Lewiner 1)**: no geometry, early-return.
//! - **Cases 7, 8, 11, 13, 14 (Chernyaev 6/7/10/12/13)**: still ambiguous
//!   and still fall back to classic MC's triangulation with a single
//!   `tracing::warn!` per run. Scheduled for follow-on sessions; see
//!   the [Domain Knowledge note](../../../../SDF_Lattice_Knowledge_Base/Domain%20Knowledge/Isosurface%20Extraction%20Methods.md)
//!   for the full roadmap.
//!
//! # Sign/index convention mismatch with the Lewiner tables
//!
//! Our internal `case_index` convention sets bit `p` iff corner `p` has
//! `f < 0` (inside). Lewiner's convention is the opposite (bit set iff
//! corner is **outside**). To use Lewiner's `cases` lookup table we
//! convert at dispatch time via `lewiner_index(our_idx) = !our_idx`.
//!
//! Edge indices in Lewiner's tables are 1-based `(1..=12)`; ours are
//! 0-based `(0..=11)`. Conversion is done at port time in each case's
//! tiling tables, not at runtime.

mod case3;
mod case4;
mod cases;
mod decider;
mod unambiguous;

use std::sync::atomic::{AtomicBool, Ordering};

use super::{CellCoord, Mesh, case_index, classic, read_corner_values};
use crate::mesh::grid::GridSpec;

/// Rate-limits the "cases 7/8/11/13/14 not yet ported" warning to one
/// emission per process run.
static PARTIAL_PORT_WARNED: AtomicBool = AtomicBool::new(false);

/// Converts our internal `case_index` (bit set iff corner is inside) to
/// Lewiner's (bit set iff corner is outside).
#[inline]
fn lewiner_index(our_idx: u8) -> u8 {
    !our_idx
}

/// Processes a single voxel and appends emitted triangles to `mesh`.
///
/// Dispatches per the Lewiner base-case table:
/// - Case 1 (empty) → no-op.
/// - Cases 2/3/6/9/10/12/15 (unambiguous) → Lewiner face-consistent tables.
/// - Case 4 → Case 3 disambiguator (face decider).
/// - Case 5 → Case 4 disambiguator (interior decider).
/// - Cases 7/8/11/13/14 (still ambiguous) → warn-once + classic fallback.
pub(crate) fn process_cell(grid: &GridSpec, field: &[f32], cell: CellCoord, mesh: &mut Mesh) {
    let corners = read_corner_values(grid, field, cell);
    let our_idx = case_index(&corners);
    let lewiner_idx = lewiner_index(our_idx);
    let (base_case, subcase) = cases::CASES[lewiner_idx as usize];

    // Empty case — no surface crosses this voxel.
    if base_case == 1 {
        return;
    }

    // Unambiguous cases use Lewiner's face-consistent tables.
    if matches!(base_case, 2 | 3 | 6 | 9 | 10 | 12 | 15) {
        unambiguous::emit(grid, cell, &corners, base_case, subcase, mesh);
        return;
    }

    // Lewiner base case 4 = Chernyaev case 3 (face-ambiguous).
    if base_case == 4 {
        case3::emit(grid, cell, &corners, subcase, mesh);
        return;
    }
    // Lewiner base case 5 = Chernyaev case 4 (body-ambiguous).
    if base_case == 5 {
        case4::emit(grid, cell, &corners, subcase, mesh);
        return;
    }

    // Remaining ambiguous cases — Lewiner 7, 8, 11, 13, 14 — are not yet
    // ported. Emit a single warning and delegate to classic.
    if is_lewiner_ambiguous(base_case) && !PARTIAL_PORT_WARNED.swap(true, Ordering::Relaxed) {
        tracing::warn!(
            "MC33 partial port: Cases 3, 4, and all unambiguous cases use Lewiner's \
             face-consistent tables. Lewiner cases 7/8/11/13/14 \
             (Chernyaev 6/7/10/12/13) still fall back to classic's triangulation \
             for this run. Output is at least as good as --extraction-method classic. \
             This warning will not repeat."
        );
    }

    classic::process_cell(grid, field, cell, mesh);
}

/// Returns `true` for Lewiner base cases that are ambiguous and not yet
/// handled by MC33 in this session.
#[inline]
fn is_lewiner_ambiguous(base_case: i8) -> bool {
    // Lewiner cases 4 (Chernyaev 3) and 5 (Chernyaev 4) are handled above
    // and return early, so we don't include them here.
    matches!(base_case, 7 | 8 | 11 | 13 | 14)
}

/// Test-only: returns the Lewiner base case for a given internal
/// `case_index`. Used by manifoldness diagnostics in the parent module.
#[cfg(test)]
pub(super) fn lewiner_base_case(our_idx: u8) -> i8 {
    cases::CASES[lewiner_index(our_idx) as usize].0
}

#[cfg(test)]
#[allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::float_cmp,
    missing_docs
)]
mod tests {
    use super::*;

    #[test]
    fn lewiner_index_is_bitwise_complement() {
        assert_eq!(lewiner_index(0), 255);
        assert_eq!(lewiner_index(255), 0);
        assert_eq!(lewiner_index(5), 250);
        // Round-trip involution.
        for i in 0..=255_u8 {
            assert_eq!(lewiner_index(lewiner_index(i)), i);
        }
    }

    #[test]
    fn cases_lookup_at_empty_configurations() {
        // our_idx 0 (all outside) → lewiner_idx 255 → (1, 0) empty.
        let (base, _) = cases::CASES[lewiner_index(0) as usize];
        assert_eq!(base, 1);
        // our_idx 255 (all inside) → lewiner_idx 0 → (1, 0) empty.
        let (base, _) = cases::CASES[lewiner_index(255) as usize];
        assert_eq!(base, 1);
    }

    #[test]
    fn cases_lookup_at_known_case3_configurations() {
        // our_idx 5 = 0b00000101: corners 0 and 2 inside, diagonally
        // opposite on the bottom face. This is the canonical Case 3
        // configuration.
        let (base, _) = cases::CASES[lewiner_index(5) as usize];
        assert_eq!(
            base, 4,
            "our_idx 5 should map to Lewiner base case 4 (Chernyaev 3)"
        );
    }

    #[test]
    fn is_lewiner_ambiguous_identifies_all_pending_cases() {
        // All cases still unhandled by the port should register as ambiguous.
        for c in [7_i8, 8, 11, 13, 14] {
            assert!(is_lewiner_ambiguous(c), "case {c} should be ambiguous");
        }
        // Non-ambiguous cases (including Cases 4 and 5, now handled) should not.
        for c in [1_i8, 2, 3, 4, 5, 6, 9, 10, 12, 15] {
            assert!(!is_lewiner_ambiguous(c), "case {c} should not be ambiguous");
        }
    }

    #[test]
    fn cases_lookup_at_known_case4_configurations() {
        // our_idx 0b10111110 = 190 → bits 1,2,3,4,5,7 inside = 6 corners
        // inside; corners 0 and 6 are outside (opposite ends of the main
        // body diagonal). Canonical Case 4 subcase-1 voxel.
        let (base, sub) = cases::CASES[lewiner_index(0b1011_1110) as usize];
        assert_eq!(
            base, 5,
            "our_idx 0b10111110 should map to Lewiner base case 5 (Chernyaev 4)"
        );
        assert_eq!(sub, 1, "expected subcase 1");
    }
}
