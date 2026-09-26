#!/usr/bin/env Rscript

# Report locked R packages that have a newer version in their repository.
# Dependabot does not read renv.lock, so the monthly CI run writes this table
# to the job summary. It does not change renv.lock or fail on findings; a
# dependency update remains a reviewed change (see the Architecture wiki page).

lock <- renv::lockfile_read("renv.lock")
available <- pak::meta_list()
available <- available[!duplicated(available[c("package", "version")]), c("package", "version")]

rows <- lapply(lock$Packages, function(record) {
  versions <- available$version[available$package == record$Package]
  if (length(versions) == 0L) {
    return(NULL)
  }
  newest <- versions[order(package_version(versions), decreasing = TRUE)][[1L]]
  if (package_version(newest) <= package_version(record$Version)) {
    return(NULL)
  }
  data.frame(
    package = record$Package,
    source = record$Source,
    locked = record$Version,
    available = newest,
    stringsAsFactors = FALSE
  )
})
outdated <- do.call(rbind, rows)

lines <- if (is.null(outdated)) {
  c("## R dependencies", "", sprintf("All %d locked packages are current.", length(lock$Packages)))
} else {
  c(
    "## R dependencies", "",
    sprintf("%d of %d locked packages have a newer version:", nrow(outdated), length(lock$Packages)), "",
    "| Package | Source | Locked | Available |", "| :--- | :--- | :--- | :--- |",
    sprintf("| %s | %s | %s | %s |", outdated$package, outdated$source, outdated$locked, outdated$available)
  )
}
writeLines(lines)
summary_path <- Sys.getenv("GITHUB_STEP_SUMMARY")
if (nzchar(summary_path)) {
  cat(lines, file = summary_path, sep = "\n", append = TRUE)
}
