//! Mesh file-format exporters.
//!
//! Each format has a dedicated module exposing a `write` function that
//! takes a [`Mesh`](super::Mesh) and a `&mut impl Write`. The CLI (or other
//! consumers) wire up a file handle to the chosen writer.
//!
//! # Supported formats (feature-gated)
//!
//! - [`stl`] — binary STL (the format `PreForm` and other slicers consume).
//!   Behind the `stl` Cargo feature (default-on).
//! - [`obj`] — Wavefront OBJ ASCII (the format every 3D viewer reads).
//!   Behind the `obj` Cargo feature (default-on).
//!
//! The [`Format`] enum and [`Format::from_extension`] are always available
//! regardless of which features are enabled — only the writer modules
//! themselves are gated. Callers attempting to invoke a writer whose
//! feature is disabled get a compile error pointing at the missing
//! module path.
//!
//! Format selection is typically driven by file extension; see the CLI's
//! `run::run` for the mapping.

#[cfg(feature = "obj")]
pub mod obj;
#[cfg(feature = "stl")]
pub mod stl;

/// Which file format to write the mesh in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Format {
    /// Binary STL — most common 3D-printing exchange format.
    Stl,
    /// Wavefront OBJ — ubiquitous in 3D viewers.
    Obj,
}

impl Format {
    /// Attempts to infer a [`Format`] from a file extension (case-insensitive).
    ///
    /// Returns `None` for unrecognized extensions; the CLI reports that as
    /// an actionable error with the list of known extensions.
    pub fn from_extension(ext: &str) -> Option<Self> {
        match ext.to_ascii_lowercase().as_str() {
            "stl" => Some(Self::Stl),
            "obj" => Some(Self::Obj),
            _ => None,
        }
    }
}

#[cfg(test)]
#[allow(clippy::unwrap_used, missing_docs)]
mod tests {
    use super::*;

    #[test]
    fn from_extension_accepts_common_extensions() {
        assert_eq!(Format::from_extension("stl"), Some(Format::Stl));
        assert_eq!(Format::from_extension("obj"), Some(Format::Obj));
    }

    #[test]
    fn from_extension_is_case_insensitive() {
        assert_eq!(Format::from_extension("STL"), Some(Format::Stl));
        assert_eq!(Format::from_extension("Obj"), Some(Format::Obj));
    }

    #[test]
    fn from_extension_rejects_unknown() {
        assert_eq!(Format::from_extension("ply"), None);
        assert_eq!(Format::from_extension(""), None);
        assert_eq!(Format::from_extension("txt"), None);
    }
}
