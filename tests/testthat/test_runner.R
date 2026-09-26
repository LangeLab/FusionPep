# Run the CLI as users do: a separate Rscript process started from an
# unrelated directory, with no global or inherited library that could supply
# renv or the analysis packages.
run_cli <- function(arguments) {
  working_dir <- tempfile("fusionpep-cwd-")
  dir.create(working_dir)
  empty_library <- tempfile("fusionpep-no-library-")
  dir.create(empty_library)
  old_dir <- setwd(working_dir)
  on.exit(setwd(old_dir), add = TRUE)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(file.path(project_root, "run_fusion_mapper.R")), shQuote(arguments)),
    stdout = TRUE,
    stderr = TRUE,
    env = c(
      "R_LIBS=",
      paste0("R_LIBS_USER=", empty_library),
      paste0("R_LIBS_SITE=", empty_library),
      "RENV_PROJECT="
    )
  ))
  list(status = attr(output, "status") %||% 0L, output = output)
}

testthat::test_that("the runner works from another directory without a global renv", {
  help <- run_cli("--help")
  testthat::expect_identical(help$status, 0L)
  testthat::expect_true(any(grepl("^Usage: Rscript run_fusion_mapper.R", help$output)))

  output_dir <- tempfile("fusionpep-cli-")
  run <- run_cli(paste0("--output=", output_dir))
  testthat::expect_identical(run$status, 0L, info = paste(run$output, collapse = "\n"))
  testthat::expect_true(file.exists(
    file.path(output_dir, "fusion_peptide_mapper_report.html")
  ))
  testthat::expect_true(file.exists(file.path(output_dir, "peptide_hits.csv")))
})

testthat::test_that("the runner keeps Windows absolute paths instead of joining them to the project", {
  run <- run_cli(c(
    "--fasta=C:/fusionpep-missing/sequences.fasta",
    paste0("--output=", tempfile("fusionpep-cli-"))
  ))
  testthat::expect_false(identical(run$status, 0L))
  # normalizePath() keeps the drive path on every platform; Windows reports it
  # with backslashes.
  testthat::expect_true(any(grepl(
    "FASTA file does not exist: C:[/\\\\]fusionpep-missing[/\\\\]sequences[.]fasta",
    run$output
  )), info = paste(run$output, collapse = "\n"))
})
