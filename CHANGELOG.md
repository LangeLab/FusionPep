<!-- markdownlint-disable MD024 -->

# FusionPep changelog

All notable changes to FusionPep are documented in this file. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and FusionPep follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html): a major version changes the command-line interface, R functions, or output files incompatibly; a minor version adds backward-compatible features; a patch version fixes bugs.

Record changes under `[Unreleased]` as they are merged. To release, rename that section to `[X.Y.Z] - YYYY-MM-DD` and push a `vX.Y.Z` tag; the section for that version becomes the GitHub release notes.

## [Unreleased]

First public version, planned as 0.1.0.

### Added

- Peptide mapping of normalized sequences to a proposed fusion protein and its two supplied parents, retaining every occurrence, overlapping matches, input row identifiers, and original metadata. I/L equivalence is the default; `--exact-il` keeps I and L distinct.
- Junction review that evaluates each supplied junction separately and reports both flank lengths against the junction's minimum flank requirement. Missing junction metadata is reported as not assessed rather than as a failed flank rule.
- Residue coverage as the union of inclusive occurrence intervals, and global and local alignments of the fusion with each parent.
- Outputs: CSV tables, an RDS result, warnings, 600 dpi PNG and vector PDF figures, and a self-contained HTML report with section navigation and file manifests. Saved outputs identify inputs by file name and MD5 hash, not by local directory.
- Command-line runner `run_fusion_mapper.R` that works from any working directory, and reusable R functions such as `run_fusion_analysis()`.
- Reproducible setup: `setup_renv.R` installs exactly the versions pinned in `renv.lock` into a project-local library without rewriting the lockfile.
- Continuous integration with lint, tests on Linux, macOS, and Windows, a 90% coverage floor, workflow security checks, and tag-driven GitHub releases.
