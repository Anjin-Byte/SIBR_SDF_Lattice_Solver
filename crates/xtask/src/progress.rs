//! Minimal progress rendering for the `xtask` post-processing commands.
//!
//! Mirrors the shape of [`crate::progress`](../../cli/src/progress.rs) in
//! the `sibr-lattice` CLI, but scoped to `xtask`'s needs — only spinners
//! (every slow op here is an opaque call into `stl_io` or `meshopt`),
//! and logging goes through plain `eprintln!` rather than `tracing`.
//!
//! All informational lines from [`remesh`](super::remesh) should be
//! emitted via [`Pipeline::println`] so they render above any live
//! spinners instead of tearing through them.

use std::io::{self, IsTerminal};
use std::time::Duration;

use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};

pub struct Pipeline {
    multi: MultiProgress,
    spinner_style: ProgressStyle,
    enabled: bool,
}

impl Pipeline {
    pub fn new() -> Self {
        let enabled = io::stderr().is_terminal();
        let multi = if enabled {
            MultiProgress::new()
        } else {
            MultiProgress::with_draw_target(ProgressDrawTarget::hidden())
        };
        let spinner_style = ProgressStyle::with_template(
            "{spinner:.cyan} {prefix:<10} {elapsed_precise} {wide_msg}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_spinner());
        Self {
            multi,
            spinner_style,
            enabled,
        }
    }

    /// Adds a steady-tick spinner with the given label. Returns a
    /// [`Spinner`] whose [`Spinner::finish`] clears the animation.
    pub fn spinner(&self, prefix: &'static str) -> Spinner {
        if !self.enabled {
            return Spinner { bar: None, prefix };
        }
        let bar = self.multi.add(ProgressBar::new_spinner());
        bar.set_style(self.spinner_style.clone());
        bar.set_prefix(prefix);
        bar.enable_steady_tick(Duration::from_millis(80));
        Spinner {
            bar: Some(bar),
            prefix,
        }
    }

    /// Prints `line` above the live spinners, respecting their draw
    /// target. Use for all informational output from `xtask` commands.
    ///
    /// `MultiProgress::println` silently drops lines when the draw
    /// target is hidden (non-TTY), so route through plain stderr in
    /// that case. When enabled, `MultiProgress::println` handles the
    /// suspend-and-redraw.
    pub fn println(&self, line: impl AsRef<str>) {
        if self.enabled {
            // In TTY mode, let MultiProgress handle it so live
            // spinners redraw cleanly after the line.
            let _ = self.multi.println(line.as_ref());
        } else {
            eprintln!("{}", line.as_ref());
        }
    }
}

impl Default for Pipeline {
    fn default() -> Self {
        Self::new()
    }
}

/// Handle to a steady-tick spinner. Call [`Spinner::finish`] explicitly
/// when the stage is complete, so the final message can be attached.
/// Drop alone does NOT finish the spinner.
pub struct Spinner {
    bar: Option<ProgressBar>,
    /// Kept for the non-TTY fallback path — `finish` prints
    /// `prefix: final_msg` to stderr when no live bar exists.
    prefix: &'static str,
}

impl Spinner {
    /// Completes the spinner, leaving a static final message in place
    /// of the animation. In non-TTY mode, prints `prefix: final_msg`
    /// to stderr so the information is still visible in logs.
    pub fn finish(self, final_msg: impl Into<std::borrow::Cow<'static, str>>) {
        let msg = final_msg.into();
        if let Some(bar) = self.bar {
            bar.finish_with_message(msg);
        } else {
            eprintln!("[xtask] {}: {}", self.prefix, msg);
        }
    }
}
