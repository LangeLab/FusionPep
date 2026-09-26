#!/usr/bin/env Rscript

# Measure line coverage of the analysis and report sources by the test suite.
# FusionPep is an R project rather than a package, so covr instruments the
# source files directly and runs the test files against them.
#
# Usage: Rscript .github/scripts/coverage.R [minimum-percent] [cobertura-path]

arguments <- commandArgs(trailingOnly = TRUE)
minimum <- as.numeric(if (length(arguments) >= 1L) arguments[[1L]] else "90")
cobertura_path <- if (length(arguments) >= 2L) arguments[[2L]] else "coverage.xml"
if (is.na(minimum) || minimum < 0 || minimum > 100) {
  stop("The minimum coverage must be a percentage.", call. = FALSE)
}

project_root <- normalizePath(".", mustWork = TRUE)
# Relative paths make the report name files as the repository does
# ("R/fusion_mapper.R"), which is how Codecov matches them to the source.
source_files <- c(
  list.files("R", pattern = "[.]R$", full.names = TRUE),
  file.path("report", "fusion_report.R")
)
# The helper file is left out: it would re-source the uninstrumented code.
test_files <- list.files(
  file.path(project_root, "tests", "testthat"),
  pattern = "^test_.*[.]R$",
  full.names = TRUE
)
test_env <- new.env(parent = globalenv())
test_env$project_root <- project_root

coverage <- covr::file_coverage(source_files, test_files, parent_env = test_env)
percent <- covr::percent_coverage(coverage)
lines <- covr::tally_coverage(coverage, by = "line")
by_file <- vapply(split(lines$value, lines$filename), function(hits) 100 * mean(hits > 0), numeric(1))
print(coverage)
covr::to_cobertura(coverage, filename = cobertura_path)

summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary_path)) {
  file_rows <- sprintf(
    "| `%s` | %.2f |",
    names(by_file),
    by_file
  )
  cat(
    "## Test coverage", "",
    sprintf("Overall line coverage: **%.2f%%** (minimum %.0f%%).", percent, minimum), "",
    "| File | Coverage (%) |", "| :--- | ---: |", file_rows, "",
    file = summary_path, sep = "\n", append = TRUE
  )
}

if (percent < minimum) {
  stop(sprintf("Coverage %.2f%% is below the %.0f%% minimum.", percent, minimum), call. = FALSE)
}
cat(sprintf("Coverage %.2f%% meets the %.0f%% minimum.\n", percent, minimum))
