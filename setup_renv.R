#!/usr/bin/env Rscript

# Initialise the project-local renv environment and install dependencies with
# pak.  This script intentionally does not call install.packages(): renv must
# already be available as the one bootstrap prerequisite.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

script_path <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_argument) == 1L) {
    return(normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE))
  }
  normalizePath(file.path(getwd(), "setup_renv.R"), mustWork = TRUE)
}

project_root <- dirname(script_path())
setwd(project_root)

if (!requireNamespace("renv", quietly = TRUE)) {
  stop(
    paste(
      "The renv package is required to bootstrap this project.",
      "It was not found in the current R installation.",
      "Install or provision renv before running setup_renv.R;",
      "this script will not modify a global library.",
      sep = "\n"
    ),
    call. = FALSE
  )
}
if (getRversion() < "4.6.0") {
  stop(
    paste(
      "This project requires R 4.6.0 or newer.",
      "The committed Bioconductor 3.23 lockfile is not compatible with older R releases.",
      sep = "\n"
    ),
    call. = FALSE
  )
}

renv_root <- file.path(project_root, "renv")
local_paths <- c(
  file.path(renv_root, "cache"),
  file.path(renv_root, "sources"),
  file.path(renv_root, "binaries"),
  file.path(renv_root, "cellar"),
  file.path(renv_root, "pak-cache"),
  file.path(renv_root, "pak-downloads")
)
for (path in local_paths) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Keep renv and pak caches inside this project.  The cache is also disabled
# below so installed packages are written directly to the project library.
Sys.setenv(
  RENV_PATHS_CACHE = file.path(renv_root, "cache"),
  RENV_PATHS_SOURCE = file.path(renv_root, "sources"),
  RENV_PATHS_BINARY = file.path(renv_root, "binaries"),
  RENV_PATHS_CELLAR = file.path(renv_root, "cellar"),
  PKG_PACKAGE_CACHE_DIR = file.path(renv_root, "pak-cache"),
  R_PKG_CACHE_DIR = file.path(renv_root, "pak-cache"),
  PKG_CACHE_DIR = file.path(renv_root, "pak-downloads"),
  RENV_CONFIG_SANDBOX_ENABLED = "FALSE",
  RENV_CONFIG_AUTO_SNAPSHOT = "FALSE"
)

activation_file <- file.path(renv_root, "activate.R")
if (!file.exists(activation_file)) {
  message("Initialising project-local renv...")
  renv::init(
    project = project_root,
    bare = TRUE,
    load = FALSE,
    restart = FALSE
  )
}

if (!file.exists(activation_file)) {
  stop("renv initialisation did not create renv/activate.R.", call. = FALSE)
}
source(activation_file, local = TRUE)

# renv::load() may restore startup variables, so reassert every local cache
# path immediately before invoking pak.
Sys.setenv(
  RENV_PATHS_CACHE = file.path(renv_root, "cache"),
  RENV_PATHS_SOURCE = file.path(renv_root, "sources"),
  RENV_PATHS_BINARY = file.path(renv_root, "binaries"),
  RENV_PATHS_CELLAR = file.path(renv_root, "cellar"),
  PKG_PACKAGE_CACHE_DIR = file.path(renv_root, "pak-cache"),
  R_PKG_CACHE_DIR = file.path(renv_root, "pak-cache"),
  PKG_CACHE_DIR = file.path(renv_root, "pak-downloads")
)
options(repos = c(CRAN = "https://cloud.r-project.org"))

renv::settings$use.cache(FALSE, project = project_root)
renv::settings$snapshot.type("explicit", project = project_root)

local_library <- renv::paths$library(project = project_root)
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)

# pak cannot install itself before it exists in the local library. Prefer an
# already available pak bootstrap tool, but install the package into the
# project library rather than using that global library as a target. If pak is
# genuinely unavailable, renv is the unavoidable bootstrap fallback.
if (!requireNamespace("pak", quietly = TRUE, lib.loc = local_library)) {
  message("Bootstrapping pak into the project-local renv library...")
  if (requireNamespace("pak", quietly = TRUE)) {
    pak::pkg_install(
      "pak",
      lib = local_library,
      upgrade = FALSE,
      ask = FALSE
    )
  } else {
    renv::install(
      packages = "pak",
      library = local_library,
      prompt = FALSE,
      project = project_root
    )
  }
}

if (!requireNamespace("pak", quietly = TRUE, lib.loc = local_library)) {
  stop(
    "pak could not be made available in the project-local renv library.",
    call. = FALSE
  )
}

locked_package_specs <- function(lockfile) {
  if (!file.exists(lockfile)) {
    return(NULL)
  }
  lock <- tryCatch(
    renv::lockfile_read(lockfile),
    error = function(error) {
      stop(
        "Could not read existing renv.lock: ",
        conditionMessage(error),
        call. = FALSE
      )
    }
  )
  records <- lock$Packages
  if (is.null(records) || length(records) == 0L) {
    return(NULL)
  }

  specs <- vapply(records, function(record) {
    package <- as.character(record$Package %||% NA_character_)
    version <- as.character(record$Version %||% NA_character_)
    source <- as.character(record$Source %||% "")
    repository <- as.character(record$Repository %||% "")
    if (is.na(package) || is.na(version) || !nzchar(package) || !nzchar(version)) {
      return(NA_character_)
    }
    prefix <- if (
      identical(source, "Bioconductor") ||
        grepl("Bioconductor", repository, fixed = TRUE)
    ) {
      "bioc"
    } else if (source %in% c("Repository", "CRAN", "standard")) {
      "cran"
    } else {
      stop(
        "Unsupported package source in renv.lock for ",
        package,
        ": ",
        source,
        call. = FALSE
      )
    }
    paste0(prefix, "::", package, "@", version)
  }, character(1))
  specs <- specs[!is.na(specs)]
  if (length(specs) == 0L) NULL else unname(specs)
}

analysis_packages <- c(
  "BiocManager",
  "bioc::BiocVersion",
  "bioc::Biostrings",
  "bioc::IRanges",
  "bioc::pwalign",
  "ggplot2",
  "htmltools",
  "base64enc",
  "testthat"
)
lockfile_path <- file.path(project_root, "renv.lock")
locked_specs <- locked_package_specs(lockfile_path)
package_specs <- if (length(locked_specs) > 0L) {
  message("Using pinned package versions from renv.lock.")
  # Explicit snapshots can omit development-only Suggests packages. Keep
  # testthat in the install set because the documented test command depends
  # on it, even when the existing lockfile predates that dependency.
  unique(c(locked_specs, "testthat"))
} else {
  analysis_packages
}

message("Installing project dependencies with pak...")
pak::pkg_install(
  package_specs,
  lib = local_library,
  upgrade = FALSE,
  ask = FALSE
)

message("Writing renv.lock...")
renv::snapshot(
  project = project_root,
  type = "explicit",
  prompt = FALSE
)

message("Project-local renv setup complete.")
message("Library: ", local_library)
message("Lockfile: ", file.path(project_root, "renv.lock"))
