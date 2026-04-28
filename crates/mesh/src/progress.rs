//! Progress-reporting seam for long-running pipeline stages.
//!
//! Stages accept a `&mut impl Progress` during their outer loops and
//! report work completed via [`Progress::set_len`] / [`Progress::inc`] /
//! [`Progress::finish`]. Binaries (the `sibr-lattice` CLI and the
//! `xtask` post-processing tool) implement this trait against a terminal
//! UI library such as `indicatif`; the library itself stays
//! UI-independent.
//!
//! # Why a trait
//!
//! The library has two concrete [`Progress`] implementations: the no-op
//! `()` impl used by tests and any headless caller, and the real
//! indicatif-backed impl that lives in the CLI crate. See the vault's
//! "Pure Core, Effectful Edges" and "Traits as Seams" patterns — this
//! trait is a minimal substitution point, not speculation.
//!
//! # Contract
//!
//! Stages MUST call `set_len` exactly once before the first `inc`, and
//! `finish` exactly once after the last `inc`. Stages with a trivial
//! amount of work (zero iterations, empty mesh) MAY skip `set_len` and
//! `inc` entirely but still MUST call `finish` so the reporting side can
//! release its UI resources.
//!
//! Implementations MUST be tolerant of being called with `total = 0`
//! (empty-work stage) and of `inc` overshooting `total` (the last tick
//! sometimes rounds up).

/// Progress reporter for a single pipeline stage.
///
/// Methods take `&mut self` so implementations can be simple stateful
/// structs without interior mutability. All three methods are infallible
/// — rendering failures are UI concerns and must not propagate into the
/// compute path.
pub trait Progress {
    /// Declare the total number of ticks this stage will emit. Called
    /// once before the outer loop.
    fn set_len(&mut self, total: u64);

    /// Advance progress by `delta` ticks.
    fn inc(&mut self, delta: u64);

    /// Mark the stage complete. Idempotent — callers may invoke it once
    /// per stage and implementations must tolerate double-calls.
    fn finish(&mut self);
}

/// No-op implementation. Used by library tests and any caller that does
/// not want progress reporting.
///
/// Concrete use: `mesh_with_progress(&job, &grid, method, &mut ())`.
impl Progress for () {
    #[inline]
    fn set_len(&mut self, _total: u64) {}
    #[inline]
    fn inc(&mut self, _delta: u64) {}
    #[inline]
    fn finish(&mut self) {}
}

/// Counting spy used by every stage's Level-1 progress test.
///
/// Accumulates calls so the test can assert that the stage declared a
/// length, ticked up to it, and finished. Lives outside the `tests`
/// module so sibling modules' unit tests can import it via
/// `crate::progress::Spy`.
#[cfg(test)]
#[derive(Debug, Default)]
pub(crate) struct Spy {
    pub(crate) set_len_calls: u32,
    pub(crate) total: u64,
    pub(crate) inc_calls: u32,
    pub(crate) inc_sum: u64,
    pub(crate) finish_calls: u32,
}

#[cfg(test)]
impl Progress for Spy {
    fn set_len(&mut self, total: u64) {
        self.set_len_calls += 1;
        self.total = total;
    }
    fn inc(&mut self, delta: u64) {
        self.inc_calls += 1;
        self.inc_sum += delta;
    }
    fn finish(&mut self) {
        self.finish_calls += 1;
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::panic, missing_docs)]
mod tests {
    use super::*;

    #[test]
    fn unit_impl_is_a_zero_cost_noop() {
        let mut p: () = ();
        p.set_len(42);
        p.inc(1);
        p.inc(99);
        p.finish();
        p.finish();
        // The () impl has no state to inspect; this test pins that the
        // impl exists, compiles, and accepts arbitrary calls.
    }

    #[test]
    fn spy_accumulates_calls() {
        let mut s = Spy::default();
        s.set_len(10);
        s.inc(3);
        s.inc(7);
        s.finish();
        assert_eq!(s.set_len_calls, 1);
        assert_eq!(s.total, 10);
        assert_eq!(s.inc_sum, 10);
        assert_eq!(s.finish_calls, 1);
    }
}
