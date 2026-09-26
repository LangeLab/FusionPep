# Troubleshooting

Start with the reported file, column, or junction ID. FusionPep rejects inconsistent inputs and records excluded peptide rows; changing a biological input to suppress an error can change the question being analyzed.

## Setup or package checks fail

Confirm that R meets the project's minimum version and can reach the network: starting R in the project directory downloads the locked `renv` version. Run `Rscript setup_renv.R` from the project root to populate the local library. This command installs the locked dependency versions without changing `renv.lock`; it is separate from analysis.

The runner checks packages in the project-local library against `renv.lock`. An installed global copy does not satisfy that check. A version mismatch means the local library needs to be reconciled with the lockfile, not that the check should be removed.

`--help` currently runs after the package checks. If help fails because a dependency is missing, complete setup first. See [Getting started](Getting-Started) for prerequisites and the current unpinned test dependency.

## A file cannot be found

The CLI resolves relative paths from the project root. Calling the script from another directory does not change that rule. Use an absolute path when the intended file is elsewhere and quote an argument containing spaces.

R function calls use the R session's working directory. Check `getwd()` when the same relative path works in the CLI but fails in your R script.

## FASTA identifiers or sequences are rejected

The default identifiers are `Fusion`, `ParentA`, and `ParentB`. Each must select exactly one different record. Keep accession details after the first whitespace or use `sequence_ids` through the [R interface](R-Interface).

Protein records must be non-empty and ungapped. Stop symbols and unsupported letters are rejected. Confirm that the supplied sequence is a protein and that the intended isoform or fusion reference was selected. Additional FASTA records are not used as a background proteome.

## A peptide column or annotation is rejected

Set `--peptide-column=NAME` to the actual sequence column. The reader needs one exact or unique case-insensitive match; duplicate or ambiguous column names cannot identify the intended field.

For excluded sequences, read `normalized_peptides.csv`. It distinguishes missing or blank values, unsupported annotations, invalid characters, and sequences that become empty after normalization. Retain the original export and correct the format using the exporting tool's documented peptide notation. Do not assume that every bracketed string is a supported modification.

## Junction coordinates are rejected

Check the coordinates against the supplied fusion sequence. The left boundary is the last upstream-parent residue; the right boundary is the first downstream-parent residue. Both are one-based.

The inserted sequence must exactly match the residues between those positions. A direct join has adjacent positions and an empty inserted sequence. Supply both corresponding parent positions or neither; supplied boundary residues must agree exactly with those parent positions.

If your new FASTA appears to be checked against the wrong junction, inspect the command. Omitting `--junctions` retains the bundled example junction file. Supply your own junction path explicitly. [Input formats](Input-Formats) lists all required fields.

## A peptide maps, but is not a junction candidate

Check three separate conditions: the occurrence crosses the supplied boundary, both flanks meet the minimum, and the peptide is absent from both supplied parents under the matching settings.

A fusion-only match can lie away from the junction. A crossing match can have a short flank. A match that passes both flanks can still occur in a supplied parent. Missing junction metadata gives an unassessed state. [Interpreting results](Interpreting-Results) explains each class.

If I/L-equivalent matching changes the expected parent sharing, inspect `match_basis` and `matched_subsequence`. Run a separate analysis with `--exact-il` only when that is the intended comparison.

## Coverage or row counts differ from expectations

All valid peptide rows contribute, including controls and rows with low upstream confidence scores. Select a restricted input table before analysis if that is the intended dataset.

Input rows, distinct matching sequences, occurrences, and junction evaluations are different counts. Repeated input peptides can share a matching sequence; one sequence can occur at several positions; one occurrence can cross more than one junction. Coverage counts each reference residue once even when several peptides overlap it.

An unassessed or excluded row also differs from a valid peptide with no match. Inspect the normalization and occurrence tables before interpreting a count difference as a missing result.

## PDF export or report links fail

The PDF figures use Cairo graphics. Before writing any output, FusionPep checks that a Cairo PDF can be written and otherwise stops with "PDF figure export needs R's Cairo graphics". On macOS, CRAN's R loads Cairo from [XQuartz](https://www.xquartz.org): install it and restart R. `capabilities("cairo")` can report `TRUE` even when the library fails to load.

The report embeds PNG figures, while PDF and table links point to other files in the output directory. Share the entire directory with its `figures/` subdirectory to retain those links.

If detail panels disappear when printing, open the panels you want before printing. Closed panels are intentionally omitted. If a figure or HTML table displays fewer rows than the CSV, read the stated display limit and total count.
