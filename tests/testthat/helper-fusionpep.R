# testthat sources helper files before the tests; load the project code once.
project_root <- normalizePath(testthat::test_path("..", ".."), mustWork = TRUE)
for (path in list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  source(path)
}
source(file.path(project_root, "report", "fusion_report.R"))
