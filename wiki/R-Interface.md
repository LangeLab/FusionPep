# R interface

Run these examples from the FusionPep project root after completing [setup](Getting-Started). The project-local R environment must be active. FusionPep is currently a script-based R project; load its functions with `source()`.

## Analyze the bundled inputs

```r
source("R/fusion_mapper.R")
source("report/fusion_report.R")

result <- run_fusion_analysis(
  sequence_file = "input/sequences.fasta",
  peptide_file = "input/peptides.csv",
  junction_file = "input/fusion_junctions.csv"
)

result$peptide_summary
result$junction_evaluations
result$coverage
```

`run_fusion_analysis()` reads and validates inputs, normalizes peptides, maps sequences, evaluates supplied junctions, calculates coverage, and aligns the fusion with each parent. It returns a `fusion_peptide_mapping_result` list without writing output files.

### Analysis arguments

The public analysis function accepts:

| Argument | Meaning | Default |
| --- | --- | --- |
| `sequence_file` | Protein FASTA path | Required |
| `peptide_file` | Peptide CSV path | Required |
| `sequence_ids` | Named mapping from roles to FASTA identifiers | `Fusion`, `ParentA`, `ParentB` |
| `peptide_column` | Peptide sequence column | `"peptide"` |
| `il_equivalent` | Match I and L equivalently | `TRUE` |
| `junction_file` | Explicit junction CSV path | `NULL` |
| `gap_opening` | Pairwise alignment gap-opening penalty | `10` |
| `gap_extension` | Pairwise alignment gap-extension penalty | `0.5` |
| `substitution_matrix` | Requested alignment matrix | `"BLOSUM62"` |
| `store_alignment_objects` | Retain raw alignment S4 objects | `FALSE` |

For custom FASTA identifiers, keep the role names `Fusion`, `ParentA`, and `ParentB` in `sequence_ids` and set their values to the respective identifiers in your file. Junction CSV role fields still use those canonical role names.

With `junction_file = NULL`, mapping and alignment proceed and fusion matches have an unassessed junction state. Unlike the CLI, this function does not default to the bundled junction file.

I/L equivalence controls peptide matching. The alignment output records exact and I/L-equivalent identity separately. Alignment summaries and the run manifest record the matrix actually used; inspect those values when interpreting alignments containing non-standard amino-acid symbols.

### Returned data

The result retains `config`, selected reference records, the original peptide table, normalized input rows, validated junctions, mapped occurrences, junction evaluations, peptide summaries, coverage, uncovered intervals, alignment tables, warnings, and the run manifest.

Raw alignment objects are omitted by default to reduce retained memory. Complete alignment columns, regions, and summaries remain available. Set `store_alignment_objects = TRUE` only when your downstream R code needs the raw S4 objects.

## Write outputs separately

Using the `result` created above:

```r
output_dir <- file.path(tempdir(), "fusionpep-example")
paths <- write_fusion_outputs(result, output_dir)
paths$report <- write_fusion_report(result, output_dir, paths)
paths$report
```

This example uses the R session's temporary directory. Choose a persistent output directory to retain the files after the session ends. Reusing a directory replaces generated files with matching names.

`write_fusion_outputs()` writes CSV, RDS, warning, and figure files and returns their paths. `write_fusion_report()` consumes the result and those paths, writes HTML, and returns the report path. Figure PDFs require Cairo graphics. [Reports and figures](Reports-and-Figures) describes the saved files.

## Reuse mapping or alignment without figure export

The lower-level functions let a script request only the calculations it needs:

- `read_fusion_sequences()` returns the selected sequences, headers, requested identifiers, source path, and input hash.
- `read_peptide_table()` returns the original input, normalized rows, selected column, source path, and input hash.
- `read_fusion_junctions()` returns validated junction rows and their source information; pass the selected sequences for coordinate checks.
- `map_peptides_to_references()` takes selected sequences, normalized rows, optional validated junction rows, and the I/L setting. It returns occurrence, evaluation, summary, coverage, and uncovered-interval tables.
- `align_fusion_to_parents()` returns the pairwise alignment columns, regions, summaries, and optional raw objects.

Use the readers to prepare inputs for these lower-level calls. Keep the same `il_equivalent` setting when normalizing and mapping. The mapping and alignment functions do not write files or create figures.

R function paths resolve relative to the current R working directory. The CLI's project-root path resolution does not apply to arbitrary R calls.
