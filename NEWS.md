# fusionpep (development version)

FusionPep has not had a tagged release. This section describes the current development state and becomes the notes for the first release.

## Analysis

- Maps normalized peptide sequences to a proposed fusion protein and its two supplied parents, retaining every occurrence, overlapping matches, input row identifiers, and original metadata. I/L equivalence is the default; `--exact-il` keeps I and L distinct.
- Evaluates each supplied junction separately and reports both flank lengths against the junction's minimum flank requirement. Missing junction metadata is reported as not assessed rather than as a failed flank rule.
- Calculates residue coverage as the union of inclusive occurrence intervals and aligns the fusion globally and locally with each parent.
- Writes CSV tables, an RDS result, warnings, 600 dpi PNG and vector PDF figures, and a self-contained HTML report with section navigation and file manifests.

## Reliability

- Setup installs exactly the versions pinned in `renv.lock` and no longer rewrites the lockfile. Bioconductor versions superseded within their release are installed from the release archive. `pak` is bootstrapped from r-lib's prebuilt binaries, so setup needs no system libraries before its first install.
- `testthat`, `covr`, and `lintr` are pinned in `renv.lock` with the analysis dependencies.
- The runner and project check work from any working directory without a global `renv` installation. The runner treats Windows drive and UNC paths as absolute.
- Saved outputs identify inputs by file name and MD5 hash instead of local directory paths. The run manifest keys are `fasta_file`, `peptide_csv_file`, and `junction_csv_file`.
- Junction CSVs read every column as text: `NA` is a valid identifier or inserted sequence (Asn-Ala), and identifiers keep leading zeros. Non-finite and out-of-range coordinates are rejected.
