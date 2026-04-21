//! `sibr-lattice` — command-line interface for SDF-based lattice generation.
//!
//! This crate is a thin orchestration layer over [`lattice_gen`]. It parses
//! arguments, constructs a validated [`lattice_gen::LatticeJob`], runs the
//! CPU reference mesher, and writes the result in a chosen file format.
//!
//! The CLI has **no compute logic of its own** — everything substantive
//! lives in `sdf` and `lattice-gen`. See the project's `cli crate`
//! architecture note.

#![forbid(unsafe_code)]

mod args;
mod run;

use std::process::ExitCode;

fn main() -> ExitCode {
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(
            tracing_subscriber::EnvFilter::try_from_default_env()
                .unwrap_or_else(|_| tracing_subscriber::EnvFilter::new("info")),
        )
        .init();

    let args = <args::Args as clap::Parser>::parse();

    match run::run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            // Print the full error chain to stderr.
            eprintln!("error: {err}");
            let mut source = err.source();
            while let Some(cause) = source {
                eprintln!("  caused by: {cause}");
                source = cause.source();
            }
            ExitCode::FAILURE
        }
    }
}
