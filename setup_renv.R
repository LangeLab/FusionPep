#!/usr/bin/env Rscript

# Initialise the project-local renv environment and install dependencies with
# pak.  This script intentionally does not call install.packages(): renv must
# already be available as the one bootstrap prerequisite.
#
# When renv.lock exists, setup installs exactly the locked versions and leaves
# the lockfile unchanged.  A lockfile is written only when none exists; changing
# dependencies is a separate, reviewed step (see the Architecture wiki page).

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

# pak cannot install itself before it exists in the local library. Bootstrap
# it from r-lib's prebuilt pak repository: those builds bundle their own
# dependencies, so no compiler or system libraries (such as libcurl headers)
# are needed. pak is the installer, not a project dependency, so renv.lock does
# not pin it.
if (!requireNamespace("pak", quietly = TRUE, lib.loc = local_library)) {
  message("Bootstrapping pak into the project-local renv library...")
  pak_repository <- sprintf(
    "https://r-lib.github.io/p/pak/stable/%s/%s/%s",
    .Platform$pkgType, R.Version()$os, R.Version()$arch
  )
  renv::install(
    packages = "pak",
    repos = c(pak = pak_repository),
    library = local_library,
    prompt = FALSE,
    project = project_root
  )
}

if (!requireNamespace("pak", quietly = TRUE, lib.loc = local_library)) {
  stop(
    "pak could not be made available in the project-local renv library.",
    call. = FALSE
  )
}

# Return install specifications for the locked packages: `pak` specs for pak,
# and `archived` specs for Bioconductor versions that are no longer current.
# pak finds bioc::pkg@version only while that version is current in its
# Bioconductor release; a later patch release moves it to the release's source
# archive.  pak installs those by archive URL, so their dependents build
# against the locked version, and renv then reinstalls them from the same
# archive to record them as Bioconductor installs matching renv.lock.
locked_package_specs <- function(lockfile, available, library) {
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
  bioc_version <- as.character(lock$Bioconductor$Version %||% NA_character_)

  specs <- vapply(records, function(record) {
    package <- as.character(record$Package %||% NA_character_)
    version <- as.character(record$Version %||% NA_character_)
    source <- as.character(record$Source %||% "")
    repository <- as.character(record$Repository %||% "")
    if (is.na(package) || is.na(version) || !nzchar(package) || !nzchar(version)) {
      stop("renv.lock contains a record without a package name or version.", call. = FALSE)
    }
    if (identical(source, "Bioconductor") ||
        grepl("Bioconductor", repository, fixed = TRUE)) {
      current <- any(available$package == package & available$version == version)
      if (!current) {
        return(paste0("archived::", package, "@", version))
      }
      # Reuse the recorded pak reference so the installed package metadata
      # matches renv.lock; the version check after installation still applies.
      return(as.character(record$RemotePkgRef %||% paste0("bioc::", package, "@", version)))
    }
    if (source %in% c("Repository", "CRAN", "standard")) {
      return(paste0("cran::", package, "@", version))
    }
    stop(
      "Unsupported package source in renv.lock for ",
      package,
      ": ",
      source,
      call. = FALSE
    )
  }, character(1))
  specs <- unname(specs)
  archived <- startsWith(specs, "archived::")
  # An archived package already installed by renv at its locked version needs
  # no work; skipping it keeps repeat setups from rebuilding it twice.
  installed_by_renv <- vapply(specs, function(spec) {
    if (!startsWith(spec, "archived::")) {
      return(FALSE)
    }
    package <- sub("^archived::([^@]+)@.*$", "\\1", spec)
    fields <- suppressWarnings(utils::packageDescription(
      package, lib.loc = library, fields = c("Version", "RemoteType")
    ))
    is.list(fields) &&
      identical(fields$Version, sub("^archived::[^@]+@", "", spec)) &&
      is.na(fields$RemoteType)
  }, logical(1))
  specs <- specs[!installed_by_renv]
  archived <- archived[!installed_by_renv]
  archived_packages <- sub("^archived::([^@]+)@.*$", "\\1", specs[archived])
  archived_versions <- sub("^archived::[^@]+@", "", specs[archived])
  if (any(archived) && is.na(bioc_version)) {
    stop("renv.lock does not record its Bioconductor version.", call. = FALSE)
  }
  archive_urls <- sprintf(
    "url::https://bioconductor.org/packages/%s/bioc/src/contrib/Archive/%s/%s_%s.tar.gz",
    bioc_version, archived_packages, archived_packages, archived_versions
  )
  list(
    pak = c(specs[!archived], archive_urls),
    archived = sub("^archived::", "bioc::", specs[archived])
  )
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
locked_specs <- locked_package_specs(lockfile_path, pak::meta_list(), local_library)
if (length(locked_specs) > 0L) {
  message("Installing pinned package versions from renv.lock with pak...")
  pak::pkg_install(locked_specs$pak, lib = local_library, upgrade = FALSE, ask = FALSE)
  if (length(locked_specs$archived) > 0L) {
    message(
      "Reinstalling archived Bioconductor versions with renv: ",
      paste(locked_specs$archived, collapse = ", ")
    )
    # Dependencies are already installed at their locked versions, so renv
    # rebuilds only these packages.
    renv::install(
      locked_specs$archived,
      library = local_library,
      prompt = FALSE,
      project = project_root
    )
  }
  # Compare the recorded version strings directly, as the analysis runner does.
  locked <- renv::lockfile_read(lockfile_path)$Packages
  installed <- vapply(names(locked), function(package) {
    tryCatch(
      as.character(utils::packageDescription(
        package, lib.loc = c(local_library, .Library), fields = "Version"
      )),
      error = function(error) NA_character_
    )
  }, character(1))
  expected <- vapply(locked, function(record) as.character(record$Version), character(1))
  mismatched <- names(locked)[is.na(installed) | installed != expected]
  if (length(mismatched) > 0L) {
    stop(
      "The project-local library does not match renv.lock after installation: ",
      paste(mismatched, collapse = ", "),
      call. = FALSE
    )
  }
} else {
  message("No renv.lock found; installing current dependency versions with pak...")
  pak::pkg_install(analysis_packages, lib = local_library, upgrade = FALSE, ask = FALSE)
  message("Writing renv.lock...")
  renv::snapshot(project = project_root, type = "explicit", prompt = FALSE)
}

message("Project-local renv setup complete.")
message("Library: ", local_library)
message("Lockfile: ", lockfile_path)
