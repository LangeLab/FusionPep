#!/usr/bin/env Rscript

# Deterministic local project gate. This script never installs packages,
# changes renv.lock, or modifies a global library.

script_path <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  if (length(file_argument) == 1L) {
    return(normalizePath(
      sub("^--file=", "", file_argument),
      mustWork = TRUE
    ))
  }
  normalizePath(file.path(getwd(), "check_project.R"), mustWork = TRUE)
}

gate_fail <- function(...) {
  stop(paste0(...), call. = FALSE)
}

project_root <- dirname(script_path())
setwd(project_root)

activation_file <- file.path(project_root, "renv", "activate.R")
if (!file.exists(activation_file)) {
  gate_fail("Missing renv activation file: ", activation_file)
}
# Startup files were read from the caller's directory, not necessarily this
# project; activate the project library explicitly, as the runner does.
Sys.setenv(RENV_PROJECT = project_root, RENV_CONFIG_SANDBOX_ENABLED = "FALSE")
source(activation_file, local = TRUE)
if (!requireNamespace("renv", quietly = TRUE)) {
  gate_fail("The renv package is unavailable; the gate will not install it.")
}

local_library <- renv::paths$library(project = project_root)
required_packages <- c(
  "Biostrings",
  "IRanges",
  "pwalign",
  "ggplot2",
  "htmltools",
  "base64enc",
  "lintr",
  "testthat"
)
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE,
    lib.loc = local_library
  )
]
if (length(missing_packages) > 0L) {
  gate_fail(
    "Missing project-local package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Run setup_renv.R; this gate does not install packages."
  )
}
.libPaths(unique(c(local_library, .libPaths())))

analysis_files <- c(
  list.files(file.path(project_root, "R"), pattern = "\\.R$", full.names = TRUE),
  list.files(
    file.path(project_root, "report"),
    pattern = "\\.R$",
    full.names = TRUE
  )
)
parse_files <- c(
  analysis_files,
  file.path(
    project_root,
    c(
      "run_fusion_mapper.R",
      "setup_renv.R",
      "check_project.R"
    )
  ),
  list.files(
    file.path(project_root, "tests"),
    pattern = "\\.R$",
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    file.path(project_root, ".github", "scripts"),
    pattern = "\\.R$",
    full.names = TRUE
  )
)
parse_files <- unique(parse_files[file.exists(parse_files)])
for (path in parse_files) {
  tryCatch(
    parse(file = path),
    error = function(error) {
      gate_fail("Could not parse ", path, ": ", conditionMessage(error))
    }
  )
}

# Style rules live in .lintr at the project root.
lint_results <- lapply(parse_files, lintr::lint)
lint_count <- sum(lengths(lint_results))
if (lint_count > 0L) {
  for (lints in lint_results[lengths(lint_results) > 0L]) {
    print(lints)
  }
  gate_fail(lint_count, " lint finding(s); see the listing above.")
}

production_files <- c(
  analysis_files,
  file.path(project_root, "run_fusion_mapper.R")
)
forbidden_patterns <- c(
  "install\\.packages\\s*\\(",
  "renv::install\\s*\\(",
  "renv::snapshot\\s*\\(",
  "browser\\s*\\(",
  "debug(?:once)?\\s*\\(",
  "recover\\s*\\(",
  "<<-"
)
for (path in production_files[file.exists(production_files)]) {
  lines <- readLines(path, warn = FALSE)
  for (pattern in forbidden_patterns) {
    matches <- which(grepl(pattern, lines, perl = TRUE))
    if (length(matches) > 0L) {
      gate_fail(
        "Forbidden production pattern '",
        pattern,
        "' in ",
        path,
        " at line(s) ",
        paste(matches, collapse = ", "),
        "."
      )
    }
  }
}

testthat::test_dir(
  file.path(project_root, "tests", "testthat"),
  reporter = "summary"
)

rscript <- file.path(R.home("bin"), "Rscript")
runner_status <- system2(
  rscript,
  args = file.path(project_root, "run_fusion_mapper.R")
)
if (!identical(as.integer(runner_status), 0L)) {
  gate_fail("The real example runner exited with status ", runner_status, ".")
}

for (path in list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(path, local = TRUE)
}
output_dir <- file.path(project_root, "results")
output_paths <- fusion_output_paths(output_dir)
missing_outputs <- names(output_paths)[
  !file.exists(unname(output_paths))
]
if (length(missing_outputs) > 0L) {
  gate_fail(
    "The example run did not produce required output(s): ",
    paste(missing_outputs, collapse = ", ")
  )
}
empty_outputs <- names(output_paths)[
  file.info(unname(output_paths))$size <= 0
]
if (length(empty_outputs) > 0L) {
  gate_fail(
    "The example run produced empty output(s): ",
    paste(empty_outputs, collapse = ", ")
  )
}

read_output_csv <- function(name) {
  utils::read.csv(
    output_paths[[name]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

required_columns <- list(
  peptide_hits = c(
    "sequence_id",
    "matching_peptide",
    "start",
    "end",
    "display_class"
  ),
  peptide_summary = c(
    "matching_peptide",
    "presence_classification",
    "fusion_evidence_class"
  ),
  junction_evaluations = c(
    "hit_row_id",
    "junction_id",
    "junction_left_flank",
    "junction_right_flank",
    "min_flank_aa",
    "junction_flank_pass"
  ),
  coverage_summary = c(
    "sequence_id",
    "covered_residues",
    "coverage_percent"
  ),
  run_manifest = c("key", "value")
)
for (name in names(required_columns)) {
  table <- read_output_csv(name)
  missing_columns <- setdiff(required_columns[[name]], names(table))
  if (length(missing_columns) > 0L) {
    gate_fail(
      "Output ",
      name,
      " is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

result <- readRDS(output_paths[["result_rds"]])
if (!inherits(result, "fusion_peptide_mapping_result")) {
  gate_fail("The result RDS does not contain a fusion_peptide_mapping_result.")
}
required_result_fields <- c(
  "sequence_input",
  "normalized_peptides",
  "junctions",
  "junction_evaluations",
  "peptide_hits",
  "peptide_summary",
  "coverage",
  "alignment_columns",
  "manifest"
)
missing_result_fields <- setdiff(required_result_fields, names(result))
if (length(missing_result_fields) > 0L) {
  gate_fail(
    "The result RDS is missing field(s): ",
    paste(missing_result_fields, collapse = ", ")
  )
}

report_text <- paste(
  readLines(output_paths[["report"]], warn = FALSE),
  collapse = "\n"
)
required_report_text <- c(
  "FusionPep: peptide mapping decision report",
  "Junction decision",
  "Annotated reference sequences",
  "Output manifest: file contents and row counts"
)
missing_report_text <- required_report_text[
  !vapply(
    required_report_text,
    grepl,
    logical(1),
    x = report_text,
    fixed = TRUE
  )
]
if (length(missing_report_text) > 0L) {
  gate_fail(
    "The report is missing required section(s): ",
    paste(missing_report_text, collapse = ", ")
  )
}

message("Project checks passed.")
