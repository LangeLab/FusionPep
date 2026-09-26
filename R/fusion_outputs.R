# Output files for FusionPep.
#
# Writes the analysis result, tables, and figures into an output directory,
# staging the files first so an interrupted run does not leave a partial set.

normalize_output_directory <- function(output_dir) {
  if (!is.character(output_dir) ||
      length(output_dir) != 1L ||
      is.na(output_dir) ||
      !nzchar(output_dir)) {
    stopf("output_dir must be one non-empty path.")
  }
  output_dir <- normalizePath(path.expand(output_dir), mustWork = FALSE)
  if (file.exists(output_dir) && !dir.exists(output_dir)) {
    stopf("output_dir is an existing file, not a directory: %s", output_dir)
  }
  output_dir
}

fusion_output_paths <- function(output_dir) {
  output_dir <- normalize_output_directory(output_dir)
  figures_dir <- file.path(output_dir, "figures")
  c(
    input_peptides = file.path(output_dir, "input_peptides.csv"),
    fusion_junctions = file.path(output_dir, "fusion_junctions.csv"),
    junction_evaluations = file.path(output_dir, "junction_evaluations.csv"),
    normalized_peptides = file.path(output_dir, "normalized_peptides.csv"),
    peptide_hits = file.path(output_dir, "peptide_hits.csv"),
    peptide_summary = file.path(output_dir, "peptide_summary.csv"),
    coverage_summary = file.path(output_dir, "coverage_summary.csv"),
    uncovered_regions = file.path(output_dir, "uncovered_regions.csv"),
    alignment_columns = file.path(output_dir, "alignment_columns.csv"),
    alignment_regions = file.path(output_dir, "alignment_regions.csv"),
    alignment_summaries = file.path(output_dir, "alignment_summaries.csv"),
    sequence_metadata = file.path(output_dir, "sequence_metadata.csv"),
    run_manifest = file.path(output_dir, "run_manifest.csv"),
    warnings = file.path(output_dir, "warnings.txt"),
    peptide_coverage_plot = file.path(figures_dir, "peptide_coverage.png"),
    junction_evidence_plot = file.path(figures_dir, "junction_evidence.png"),
    alignment_status_plot = file.path(figures_dir, "alignment_status.png"),
    peptide_coverage_pdf = file.path(figures_dir, "peptide_coverage.pdf"),
    junction_evidence_pdf = file.path(figures_dir, "junction_evidence.pdf"),
    alignment_status_pdf = file.path(figures_dir, "alignment_status.pdf"),
    result_rds = file.path(output_dir, "fusion_peptide_mapper_result.rds"),
    report = file.path(output_dir, "fusion_peptide_mapper_report.html")
  )
}

assert_output_paths_do_not_overwrite_inputs <- function(result, output_paths) {
  input_paths <- c(
    result$sequence_input$source_path,
    result$peptide_input$source_path,
    result$config$junction_file
  )
  input_paths <- input_paths[!is.na(input_paths) & nzchar(input_paths)]
  input_paths <- normalizePath(input_paths, mustWork = FALSE)
  target_paths <- normalizePath(unname(output_paths), mustWork = FALSE)
  collisions <- target_paths[target_paths %in% input_paths]
  if (length(collisions) > 0L) {
    stopf(
      "Output would overwrite input file(s): %s. Choose a different output directory.",
      paste(basename(collisions), collapse = ", ")
    )
  }
  invisible(TRUE)
}

write_csv_output <- function(data, path) {
  utils::write.csv(
    data,
    file = path,
    row.names = FALSE,
    na = "",
    quote = TRUE
  )
}

write_fusion_outputs <- function(result, output_dir) {
  if (!inherits(result, "fusion_peptide_mapping_result")) {
    stopf("write_fusion_outputs expects a fusion_peptide_mapping_result.")
  }
  output_dir <- normalize_output_directory(output_dir)

  final_paths <- fusion_output_paths(output_dir)
  assert_output_paths_do_not_overwrite_inputs(result, final_paths)
  dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
  staging_dir <- tempfile(
    pattern = ".fusionpep-",
    tmpdir = dirname(output_dir)
  )
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)
  artifact_paths <- fusion_output_paths(staging_dir)

  figures_dir <- file.path(staging_dir, "figures")
  dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

  sequence_metadata <- data.frame(
    sequence_id = names(result$sequence_input$sequence_text),
    fasta_header = unname(result$sequence_input$source_headers),
    sequence_length = nchar(result$sequence_input$sequence_text),
    stringsAsFactors = FALSE
  )

  csv_paths <- artifact_paths[c(
    "input_peptides",
    "fusion_junctions",
    "junction_evaluations",
    "normalized_peptides",
    "peptide_hits",
    "peptide_summary",
    "coverage_summary",
    "uncovered_regions",
    "alignment_columns",
    "alignment_regions",
    "alignment_summaries",
    "sequence_metadata",
    "run_manifest"
  )]
  write_csv_output(result$peptide_input$input, csv_paths[["input_peptides"]])
  write_csv_output(result$junctions, csv_paths[["fusion_junctions"]])
  write_csv_output(
    result$junction_evaluations,
    csv_paths[["junction_evaluations"]]
  )
  write_csv_output(result$normalized_peptides, csv_paths[["normalized_peptides"]])
  write_csv_output(result$peptide_hits, csv_paths[["peptide_hits"]])
  write_csv_output(result$peptide_summary, csv_paths[["peptide_summary"]])
  write_csv_output(result$coverage, csv_paths[["coverage_summary"]])
  write_csv_output(result$uncovered_regions, csv_paths[["uncovered_regions"]])
  write_csv_output(result$alignment_columns, csv_paths[["alignment_columns"]])
  write_csv_output(result$alignment_regions, csv_paths[["alignment_regions"]])
  write_csv_output(result$alignment_summaries, csv_paths[["alignment_summaries"]])
  write_csv_output(sequence_metadata, csv_paths[["sequence_metadata"]])
  write_csv_output(result$manifest, csv_paths[["run_manifest"]])
  saveRDS(result, artifact_paths[["result_rds"]])

  warnings_path <- artifact_paths[["warnings"]]
  if (length(result$warnings) == 0L) {
    writeLines("No normalization warnings.", warnings_path, useBytes = TRUE)
  } else {
    writeLines(result$warnings, warnings_path, useBytes = TRUE)
  }

  coverage_plot_path <- artifact_paths[["peptide_coverage_plot"]]
  coverage_plot <- plot_peptide_coverage(
    result$peptide_hits,
    result$sequence_input$sequence_text,
    result$junctions
  )
  write_publication_figure(
    coverage_plot, coverage_plot_path,
    artifact_paths[["peptide_coverage_pdf"]], height = 3.7
  )
  rm(coverage_plot)

  junction_plot_path <- artifact_paths[["junction_evidence_plot"]]
  junction_plot <- plot_junction_evidence(
    result$peptide_hits,
    result$sequence_input$sequence_text,
    result$junctions,
    max_candidates = 12L
  )
  junction_counts <- tabulate(
    match(result$junction_evaluations$junction_id, result$junctions$junction_id),
    nbins = nrow(result$junctions)
  )
  junction_height <- max(
    3.5, 1.9 + sum(0.8 + 0.48 * pmin(junction_counts, 12L))
  )
  write_publication_figure(
    junction_plot, junction_plot_path,
    artifact_paths[["junction_evidence_pdf"]], height = junction_height
  )
  rm(junction_plot)

  alignment_plot_path <- artifact_paths[["alignment_status_plot"]]
  alignment_plot <- plot_alignment_status(
    result$alignment_columns,
    result$junctions
  )
  write_publication_figure(
    alignment_plot, alignment_plot_path,
    artifact_paths[["alignment_status_pdf"]], height = 3.3
  )

  promoted_names <- c(
    names(csv_paths),
    "result_rds",
    "warnings",
    "peptide_coverage_plot",
    "junction_evidence_plot",
    "alignment_status_plot",
    "peptide_coverage_pdf",
    "junction_evidence_pdf",
    "alignment_status_pdf"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(
    file.path(output_dir, "figures"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  for (name in promoted_names) {
    source_path <- artifact_paths[[name]]
    target_path <- final_paths[[name]]
    moved <- file.rename(source_path, target_path)
    if (!moved) {
      copied <- file.copy(source_path, target_path, overwrite = TRUE)
      if (!copied) {
        stopf("Could not finalize output artifact: %s", target_path)
      }
      unlink(source_path)
    }
  }

  c(
    as.list(final_paths[names(csv_paths)]),
    result_rds = final_paths[["result_rds"]],
    warnings = final_paths[["warnings"]],
    peptide_coverage_plot = final_paths[["peptide_coverage_plot"]],
    junction_evidence_plot = final_paths[["junction_evidence_plot"]],
    alignment_status_plot = final_paths[["alignment_status_plot"]],
    peptide_coverage_pdf = final_paths[["peptide_coverage_pdf"]],
    junction_evidence_pdf = final_paths[["junction_evidence_pdf"]],
    alignment_status_pdf = final_paths[["alignment_status_pdf"]]
  )
}
