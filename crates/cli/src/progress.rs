//! Terminal progress rendering for the meshing pipeline.
//!
//! This module is the **only** place in the `sibr-lattice` binary that
//! mentions `indicatif`. Elsewhere the rest of the CLI (and the whole
//! `lattice-gen` library) sees only the [`lattice_gen::Progress`] trait.
//!
//! # What lives here
//!
//! - [`Pipeline`] — owns a shared [`MultiProgress`] and the two
//!   [`ProgressStyle`]s used by all bars. Constructed once in
//!   [`crate::run::run`] and handed out bars via [`Pipeline::stage`]
//!   (determinate) and [`Pipeline::spinner`] (indeterminate).
//! - [`Bar`] — the trait-implementing wrapper returned to callers. It
//!   forwards `set_len` / `inc` / `finish` to its inner `ProgressBar`,
//!   or no-ops when the terminal is non-interactive.
//! - [`LogWriter`] — a [`std::io::Write`] adapter passed to
//!   `tracing-subscriber` so log lines suspend the bars around each
//!   write. This is what keeps the output legible when both bars and
//!   `tracing::info!` lines are active on the same stderr.
//!
//! # TTY gating
//!
//! When stderr is not a terminal (pipe, redirect, CI), [`Pipeline::new`]
//! sets `enabled = false` and every [`Bar`] returned is a no-op. Tracing
//! output flows straight through, unchanged. No per-call-site
//! conditionals are needed in [`crate::run`].

use std::io::{self, IsTerminal, Write};
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};

/// A pipeline-wide holder for the shared progress drawing surface and
/// the per-bar styles.
///
/// Construct once at the top of a CLI action, hand out [`Bar`]s per
/// stage, and pass the pipeline's log-writer factory to
/// `tracing-subscriber`. Cloning is cheap — indicatif's `MultiProgress`
/// is `Arc`-backed.
pub struct Pipeline {
    multi: MultiProgress,
    determinate_style: ProgressStyle,
    spinner_style: ProgressStyle,
    enabled: bool,
}

impl Pipeline {
    /// Builds a pipeline that renders if and only if stderr is an
    /// interactive terminal. Piped or redirected runs (CI, log capture)
    /// produce silent bars.
    pub fn new() -> Self {
        let enabled = io::stderr().is_terminal();
        // When disabled, route drawing to a hidden target so any stray
        // `inc` calls after we've handed out no-op bars also stay
        // quiet. Belt and suspenders — the `Bar` wrapper already
        // short-circuits.
        let multi = if enabled {
            MultiProgress::new()
        } else {
            MultiProgress::with_draw_target(ProgressDrawTarget::hidden())
        };

        let determinate_style = ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:<12} [{bar:30.cyan/blue}] {pos:>7}/{len:<7} {elapsed_precise} eta {eta_precise} {wide_msg}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar())
        .progress_chars("=>-");

        let spinner_style = ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:<12} {elapsed_precise} {wide_msg}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_spinner());

        Self {
            multi,
            determinate_style,
            spinner_style,
            enabled,
        }
    }

    /// Creates a determinate bar with the given label. Call
    /// `set_len` before ticking; that happens automatically if you pass
    /// this `Bar` to a `*_with_progress` library function.
    pub fn stage(&self, prefix: &'static str) -> Bar {
        if !self.enabled {
            return Bar(None);
        }
        let bar = self.multi.add(ProgressBar::new(0));
        bar.set_style(self.determinate_style.clone());
        bar.set_prefix(prefix);
        Bar(Some(bar))
    }

    /// Creates an indeterminate spinner with the given label. Suitable
    /// for opaque work (FFI calls, fast bookkeeping passes where a
    /// determinate bar would distract). Ticks itself every 80 ms until
    /// finished.
    pub fn spinner(&self, prefix: &'static str) -> Bar {
        if !self.enabled {
            return Bar(None);
        }
        let bar = self.multi.add(ProgressBar::new_spinner());
        bar.set_style(self.spinner_style.clone());
        bar.set_prefix(prefix);
        bar.enable_steady_tick(Duration::from_millis(80));
        Bar(Some(bar))
    }

    /// Returns a [`LogWriter`] suitable for
    /// `tracing_subscriber::fmt().with_writer(...)`. Every tracing line
    /// written through it suspends the bars so the output stays legible.
    pub fn log_writer(&self) -> LogWriter {
        LogWriter {
            multi: self.multi.clone(),
        }
    }
}

impl Default for Pipeline {
    fn default() -> Self {
        Self::new()
    }
}

/// A single pipeline stage's progress reporter. Implements
/// [`lattice_gen::Progress`] so it can be passed directly to
/// `*_with_progress` library calls.
///
/// `Bar(None)` is the non-TTY / disabled variant — every method is a
/// cheap no-op.
pub struct Bar(Option<ProgressBar>);

impl Bar {
    /// Attach a status message shown after the bar (e.g., "43 MB
    /// written"). No-op when disabled.
    #[allow(dead_code)] // Reserved for later polish; harmless to keep.
    pub fn message(&self, msg: impl Into<std::borrow::Cow<'static, str>>) {
        if let Some(bar) = &self.0 {
            bar.set_message(msg);
        }
    }
}

impl lattice_gen::Progress for Bar {
    fn set_len(&mut self, total: u64) {
        if let Some(bar) = &self.0 {
            bar.set_length(total);
        }
    }
    fn inc(&mut self, delta: u64) {
        if let Some(bar) = &self.0 {
            bar.inc(delta);
        }
    }
    fn finish(&mut self) {
        if let Some(bar) = &self.0 {
            bar.finish();
        }
    }
}

/// A [`std::io::Write`] adapter that suspends the shared
/// [`MultiProgress`] for each write. Used as
/// `tracing_subscriber::fmt().with_writer(...)` so bars and log lines
/// share stderr without colliding.
///
/// Each `write` call is a single tracing event (tracing-subscriber
/// buffers the whole formatted line before calling through), so
/// suspending per write is correct granularity.
#[derive(Clone)]
pub struct LogWriter {
    multi: MultiProgress,
}

impl Write for LogWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // `MultiProgress::suspend` redraws the bars after `f` returns,
        // so writing directly to stderr underneath is safe.
        self.multi.suspend(|| {
            let mut err = io::stderr().lock();
            err.write_all(buf).map(|()| buf.len())
        })
    }

    fn flush(&mut self) -> io::Result<()> {
        self.multi.suspend(|| io::stderr().lock().flush())
    }
}

impl<'a> tracing_subscriber::fmt::MakeWriter<'a> for LogWriter {
    type Writer = LogWriter;
    fn make_writer(&'a self) -> Self::Writer {
        self.clone()
    }
}
