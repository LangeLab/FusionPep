#!/usr/bin/env Rscript

# Routine entry point.  Change the configuration block below when replacing
# the input files, or use the documented command-line overrides.

script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) == 1L) {
    return(normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE))
  }
  normalizePath(file.path(getwd(), "run_fusion_mapper.R"), mustWork = TRUE)
}

PROJECT_ROOT <- dirname(script_path())
INPUT_FASTA <- file.path(PROJECT_ROOT, "input", "sequences.fasta")
PEPTIDE_CSV <- file.path(PROJECT_ROOT, "input", "peptides.csv")
JUNCTION_CSV <- file.path(PROJECT_ROOT, "input", "fusion_junctions.csv")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "results")
SEQUENCE_IDS <- c(
  Fusion = "Fusion",
  ParentA = "ParentA",
  ParentB = "ParentB"
)
PEPTIDE_COLUMN <- "peptide"
I_L_EQUIVALENT <- TRUE

print_help <- function() {
  cat(
    paste(
      "FusionPep: peptide mapping and junction review for fusion proteins.",
      "",
      "Usage: Rscript run_fusion_mapper.R [options]",
      "",
      "Options:",
      "  --fasta=PATH            Protein FASTA (default: input/sequences.fasta).",
      "  --peptides=PATH         Peptide CSV (default: input/peptides.csv).",
      "  --junctions=PATH        Junction CSV (default: input/fusion_junctions.csv).",
      "  --output=PATH           Output directory (default: results/).",
      "  --peptide-column=NAME   Peptide sequence column (default: peptide).",
      "  --exact-il              Keep I and L distinct (default: I/L equivalence).",
      "  --help, -h              Show this message.",
      "",
      "FASTA identifiers: Fusion, ParentA, ParentB; other records are not compared.",
      "Relative paths resolve from the project root. Use --name=value syntax.",
      "For your own analysis, supply all three input paths; omitted inputs keep",
      "their bundled defaults. Reusing an output directory replaces generated files.",
      "",
      "Example: Rscript run_fusion_mapper.R --output=results/example",
      "Sequence matches do not establish peptide detection or fusion expression.",
      sep = "\n"
    ),
    "\n",
    sep = ""
  )
}

parse_arguments <- function(arguments) {
  values <- list()
  for (argument in arguments) {
    if (identical(argument, "--help") || identical(argument, "-h")) {
      print_help()
      quit(save = "no", status = 0L)
    }
    if (identical(argument, "--exact-il")) {
      if (identical(values$il_equivalent, FALSE)) {
        stop("I/L matching mode was specified more than once.", call. = FALSE)
      }
      values$il_equivalent <- FALSE
      next
    }
    if (!grepl("^--[^=]+=.+$", argument, perl = TRUE)) {
      stop(
        "Unrecognised argument: ",
        argument,
        ". Use --help for supported options.",
        call. = FALSE
      )
    }
    key_value <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1]]
    key <- key_value[[1]]
    value <- paste(key_value[-1L], collapse = "=")
    key <- switch(
      key,
      fasta = "fasta",
      peptides = "peptides",
      junctions = "junctions",
      output = "output",
      `peptide-column` = "peptide_column",
      stop("Unrecognised option: --", key, call. = FALSE)
    )
    if (!nzchar(value)) {
      stop("Option --", gsub("_", "-", key), " cannot be empty.", call. = FALSE)
    }
    if (!is.null(values[[key]])) {
      stop(
        "Option --",
        gsub("_", "-", key),
        " was specified more than once.",
        call. = FALSE
      )
    }
    values[[key]] <- value
  }
  values
}

resolve_project_path <- function(path, project_root) {
  if (is.null(path)) {
    return(NULL)
  }
  if (grepl("^(/|~[/])", path)) {
    return(normalizePath(path.expand(path), mustWork = FALSE))
  }
  normalizePath(file.path(project_root, path), mustWork = FALSE)
}

activation_file <- file.path(PROJECT_ROOT, "renv", "activate.R")
if (!file.exists(activation_file)) {
  stop(
    paste(
      "This project has not been initialised with renv.",
      "Run Rscript setup_renv.R first.",
      sep = "\n"
    ),
    call. = FALSE
  )
}
Sys.setenv(RENV_PROJECT = PROJECT_ROOT)
local_library <- renv::paths$library(project = PROJECT_ROOT)
normalised_library <- normalizePath(local_library, mustWork = FALSE)
normalised_paths <- normalizePath(.libPaths(), mustWork = FALSE)
if (!normalised_library %in% normalised_paths) {
  source(activation_file, local = TRUE)
}

local_library <- renv::paths$library(project = PROJECT_ROOT)
required_packages <- c(
  "Biostrings",
  "IRanges",
  "pwalign",
  "ggplot2",
  "htmltools",
  "base64enc"
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
  stop(
    paste(
      "The project-local renv library is incomplete.",
      paste0("Missing: ", paste(missing_packages, collapse = ", ")),
      "Run Rscript setup_renv.R; the analysis runner never installs packages.",
      sep = "\n"
    ),
    call. = FALSE
  )
}
.libPaths(unique(c(local_library, .libPaths())))

lockfile_path <- file.path(PROJECT_ROOT, "renv.lock")
if (file.exists(lockfile_path)) {
  lock <- renv::lockfile_read(lockfile_path)
  locked_packages <- lock$Packages
  version_mismatches <- vapply(names(locked_packages), function(package) {
    record <- locked_packages[[package]]
    expected <- if (is.null(record$Version) || length(record$Version) == 0L) {
      NA_character_
    } else {
      as.character(record$Version)
    }
    actual <- tryCatch(
      as.character(utils::packageDescription(
        package,
        lib.loc = local_library,
        fields = "Version"
      )),
      error = function(error) NA_character_
    )
    is.na(expected) || is.na(actual) || !identical(actual, expected)
  }, logical(1))
  if (any(version_mismatches)) {
    mismatched <- names(locked_packages)[version_mismatches]
    stop(
      paste(
        "The project-local library does not match renv.lock.",
        paste("Mismatched or missing:", paste(mismatched, collapse = ", ")),
        "Run Rscript setup_renv.R before analysis.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }
}

source(file.path(PROJECT_ROOT, "R", "fusion_mapper.R"), local = TRUE)
source(file.path(PROJECT_ROOT, "report", "fusion_report.R"), local = TRUE)

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
input_fasta <- resolve_project_path(arguments$fasta %||% INPUT_FASTA, PROJECT_ROOT)
peptide_csv <- resolve_project_path(arguments$peptides %||% PEPTIDE_CSV, PROJECT_ROOT)
junction_csv <- resolve_project_path(arguments$junctions %||% JUNCTION_CSV, PROJECT_ROOT)
output_dir <- resolve_project_path(arguments$output %||% OUTPUT_DIR, PROJECT_ROOT)
peptide_column <- arguments$peptide_column %||% PEPTIDE_COLUMN
il_equivalent <- arguments$il_equivalent %||% I_L_EQUIVALENT

message("Running FusionPep...")
message("  FASTA:    ", input_fasta)
message("  Peptides: ", peptide_csv)
message("  Junctions: ", junction_csv)
message("  Output:   ", output_dir)

result <- run_fusion_analysis(
  sequence_file = input_fasta,
  peptide_file = peptide_csv,
  junction_file = junction_csv,
  sequence_ids = SEQUENCE_IDS,
  peptide_column = peptide_column,
  il_equivalent = il_equivalent
)
artifact_paths <- write_fusion_outputs(result, output_dir)
artifact_paths$report <- write_fusion_report(
  result,
  output_dir,
  artifact_paths
)

message("Completed.")
message("  Report:   ", artifact_paths$report)
message("  Coverage: ", artifact_paths$coverage_summary)
message("  Hits:     ", artifact_paths$peptide_hits)
