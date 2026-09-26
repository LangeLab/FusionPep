# Core functions for FusionPep: input readers, normalization, matching,
# junction evaluation, coverage, alignment, and the analysis pipeline.
# Figures are drawn in fusion_figures.R and files written in fusion_outputs.R.
#
# The functions in this file deliberately use explicit namespaces.  The
# runner activates renv before sourcing this file and does not install
# packages while analysing a dataset.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

stopf <- function(...) {
  stop(sprintf(...), call. = FALSE)
}

require_data_frame_columns <- function(data, columns, label) {
  if (!is.data.frame(data)) {
    stopf("%s must be a data frame.", label)
  }
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0L) {
    stopf(
      "%s is missing required column(s): %s",
      label,
      paste(missing, collapse = ", ")
    )
  }
  invisible(data)
}

validate_logical_option <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stopf("%s must be one non-missing TRUE/FALSE value.", name)
  }
  invisible(TRUE)
}

validate_nonnegative_number <- function(value, name) {
  if (!is.numeric(value) ||
      length(value) != 1L ||
      is.na(value) ||
      !is.finite(value) ||
      value < 0) {
    stopf("%s must be one finite non-negative number.", name)
  }
  invisible(TRUE)
}

empty_hit_table <- function() {
  data.frame(
    sequence_id = character(),
    matching_peptide = character(),
    normalized_peptide = character(),
    peptide_input_values = character(),
    input_row_ids = character(),
    input_count = integer(),
    start = integer(),
    end = integer(),
    length = integer(),
    matched_subsequence = character(),
    match_basis = character(),
    crosses_junction = logical(),
    junction_assessed = logical(),
    junction_flank_pass = logical(),
    junction_id = character(),
    junction_left_flank = integer(),
    junction_right_flank = integer(),
    junction_min_flank_aa = integer(),
    display_class = character(),
    stringsAsFactors = FALSE
  )
}

empty_alignment_columns <- function() {
  data.frame(
    parent = character(),
    alignment_type = character(),
    alignment_column = integer(),
    fusion_pos = integer(),
    parent_pos = integer(),
    fusion_aa = character(),
    parent_aa = character(),
    raw_match = logical(),
    il_equivalent_match = logical(),
    status = character(),
    stringsAsFactors = FALSE
  )
}

empty_alignment_regions <- function() {
  data.frame(
    parent = character(),
    alignment_type = character(),
    region_id = integer(),
    status = character(),
    fusion_start = integer(),
    fusion_end = integer(),
    parent_start = integer(),
    parent_end = integer(),
    alignment_columns = integer(),
    stringsAsFactors = FALSE
  )
}

empty_uncovered_regions <- function() {
  data.frame(
    sequence_id = character(),
    start = integer(),
    end = integer(),
    length = integer(),
    stringsAsFactors = FALSE
  )
}

empty_interval_table <- function() {
  data.frame(
    start = integer(),
    end = integer(),
    length = integer(),
    stringsAsFactors = FALSE
  )
}

empty_junction_table <- function() {
  data.frame(
    fusion_id = character(),
    junction_id = character(),
    upstream_parent = character(),
    downstream_parent = character(),
    fusion_left_position = integer(),
    fusion_right_position = integer(),
    upstream_parent_position = integer(),
    downstream_parent_position = integer(),
    inserted_sequence = character(),
    min_flank_aa = integer(),
    source = character(),
    notes = character(),
    stringsAsFactors = FALSE
  )
}

empty_junction_evaluations <- function() {
  data.frame(
    hit_row_id = integer(),
    sequence_id = character(),
    matching_peptide = character(),
    junction_id = character(),
    junction_left_flank = integer(),
    junction_right_flank = integer(),
    min_flank_aa = integer(),
    junction_flank_pass = logical(),
    stringsAsFactors = FALSE
  )
}

standard_amino_acids <- strsplit("ACDEFGHIKLMNPQRSTVWYBXZJUO", "", fixed = TRUE)[[1]]
invalid_amino_acid_pattern <- paste0(
  "[^",
  paste(standard_amino_acids, collapse = ""),
  "]"
)
supported_modification_names <- c(
  "oxidation", "carbamidomethyl", "phospho", "acetyl", "methyl",
  "dimethyl", "trimethyl", "deamidated", "glygly", "pyro", "pyro-glu",
  "amidated", "biotin", "sulfo", "nitro", "ubiquitin", "carboxy",
  "propionyl"
)
supported_modification_name_pattern <- paste(
  supported_modification_names,
  collapse = "|"
)

# Sequence validation and input normalization

# Biostrings' PDict matcher is implemented for nucleotide alphabets, not
# AAStringSet. Encode each supported amino-acid symbol as a
# unique three-base DNA word so that a dictionary can still be used for large
# peptide tables. Three bases provide 64 codes, enough for every accepted
# amino-acid symbol. Coordinates are converted back to amino-acid positions
# immediately after matching.
protein_dictionary_codebook <- local({
  alphabet <- unique(standard_amino_acids)
  dna_words <- apply(
    expand.grid(
      c("A", "C", "G", "T"),
      c("A", "C", "G", "T"),
      c("A", "C", "G", "T"),
      stringsAsFactors = FALSE
    ),
    1L,
    paste0,
    collapse = ""
  )
  setNames(dna_words[seq_along(alphabet)], alphabet)
})

encode_protein_for_dictionary <- function(sequence_text) {
  chars <- strsplit(sequence_text, "", fixed = TRUE)[[1]]
  encoded <- unname(protein_dictionary_codebook[chars])
  if (anyNA(encoded)) {
    stopf("Cannot encode unsupported amino-acid symbols for peptide matching.")
  }
  paste0(encoded, collapse = "")
}

validate_sequence_text <- function(sequence_text, sequence_id) {
  if (length(sequence_text) != 1L || is.na(sequence_text) || !nzchar(sequence_text)) {
    stopf("Sequence '%s' is empty.", sequence_id)
  }

  sequence_text <- toupper(sequence_text)
  invalid_positions <- gregexpr(
    invalid_amino_acid_pattern,
    sequence_text,
    perl = TRUE
  )[[1L]]
  if (length(invalid_positions) > 0L && invalid_positions[[1L]] != -1L) {
    invalid <- sort(unique(substring(
      sequence_text,
      invalid_positions,
      invalid_positions
    )))
    stopf(
      "Sequence '%s' contains unsupported character(s): %s. ",
      sequence_id,
      paste(invalid, collapse = ", ")
    )
  }

  invisible(TRUE)
}

read_raw_fasta_records <- function(path) {
  lines <- tryCatch(
    readLines(path, warn = FALSE, encoding = "UTF-8"),
    error = function(error) {
      stopf("Could not read FASTA file '%s': %s", path, conditionMessage(error))
    }
  )
  if (anyNA(lines)) {
    stopf("FASTA file contains invalid text encoding: %s", path)
  }
  if (length(lines) > 0L) {
    lines[[1L]] <- sub("^\uFEFF", "", lines[[1L]])
  }
  nonblank <- which(nzchar(trimws(lines)))
  if (length(nonblank) == 0L) {
    stopf("FASTA file contains no records: %s", path)
  }

  header_rows <- which(startsWith(lines, ">"))
  if (length(header_rows) == 0L || header_rows[[1L]] != nonblank[[1L]]) {
    stopf("FASTA file must begin with a record header: %s", path)
  }

  headers <- substring(lines[header_rows], 2L)
  if (any(!nzchar(trimws(headers)))) {
    stopf("Every FASTA record must have a non-empty header.")
  }
  sequences <- vapply(seq_along(header_rows), function(index) {
    start <- header_rows[[index]] + 1L
    end <- if (index < length(header_rows)) {
      header_rows[[index + 1L]] - 1L
    } else {
      length(lines)
    }
    sequence_lines <- if (start <= end) lines[start:end] else character()
    sequence_lines <- sequence_lines[nzchar(trimws(sequence_lines))]
    gsub("[[:space:]]+", "", paste(sequence_lines, collapse = ""), perl = TRUE)
  }, character(1))

  list(headers = headers, sequences = sequences)
}

read_fusion_sequences <- function(path,
                                  sequence_ids = c(
                                    Fusion = "Fusion",
                                    ParentA = "ParentA",
                                    ParentB = "ParentB"
                                  )) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stopf("FASTA file does not exist: %s", path %||% "<missing>")
  }
  if (is.null(names(sequence_ids)) ||
      !all(c("Fusion", "ParentA", "ParentB") %in% names(sequence_ids))) {
    stopf("sequence_ids must be a named character vector containing Fusion, ParentA, and ParentB.")
  }
  sequence_ids <- sequence_ids[c("Fusion", "ParentA", "ParentB")]
  if (anyNA(sequence_ids) || any(!nzchar(sequence_ids))) {
    stopf("All configured FASTA record names must be non-empty.")
  }

  raw_records <- read_raw_fasta_records(path)
  headers <- raw_records$headers
  if (anyDuplicated(headers)) {
    duplicate_headers <- unique(headers[duplicated(headers)])
    stopf(
      "FASTA contains duplicate record headers: %s",
      paste(duplicate_headers, collapse = ", ")
    )
  }

  first_tokens <- sub("\\s.*$", "", headers, perl = TRUE)
  selected_indices <- integer(length(sequence_ids))
  for (role in names(sequence_ids)) {
    requested <- sequence_ids[[role]]
    candidates <- which(headers == requested | first_tokens == requested)
    if (length(candidates) == 0L) {
      stopf(
        "Could not find FASTA record '%s' for %s. Available records: %s",
        requested,
        role,
        paste(headers, collapse = ", ")
      )
    }
    if (length(candidates) > 1L) {
      stopf(
        "FASTA name '%s' matches multiple records for %s.",
        requested,
        role
      )
    }
    selected_indices[[which(names(sequence_ids) == role)]] <- candidates
  }
  if (anyDuplicated(selected_indices)) {
    stopf("Fusion, ParentA, and ParentB must refer to three different FASTA records.")
  }

  sequence_text <- toupper(raw_records$sequences[selected_indices])
  names(sequence_text) <- names(sequence_ids)
  for (role in names(sequence_text)) {
    validate_sequence_text(sequence_text[[role]], role)
  }

  path <- normalizePath(path, mustWork = TRUE)

  list(
    sequence_text = sequence_text,
    requested_ids = sequence_ids,
    source_path = path,
    source_md5 = unname(as.character(tools::md5sum(path))),
    source_headers = setNames(headers[selected_indices], names(sequence_ids))
  )
}

supported_modification_token <- function(token) {
  token <- gsub("\\s+", "", tolower(trimws(token)), perl = TRUE)
  if (!nzchar(token)) {
    return(FALSE)
  }
  common_name <- grepl(
    paste0(
      "^(?:",
      supported_modification_name_pattern,
      ")(?:\\([^()]+\\))?$"
    ),
    token,
    perl = TRUE
  )
  numeric_shift <- grepl(
    "^[+-]?[0-9]+(?:\\.[0-9]+)?(?:da)?$",
    token,
    perl = TRUE
  )
  shorthand <- grepl("^(?:ox|cam|p|ac)$", token, perl = TRUE)
  unimod_id <- grepl("^unimod:[1-9][0-9]{0,3}$", token, perl = TRUE)
  unimod_name <- grepl(
    paste0(
      "^unimod:(?:",
      supported_modification_name_pattern,
      ")(?:\\([^()]+\\))?$"
    ),
    token,
    perl = TRUE
  )
  isotope_label <- grepl(
    "^label:(?:(?:13c|15n|18o|2h|d|t)\\([0-9]+\\))+$",
    token,
    perl = TRUE
  )
  reporter_label <- grepl(
    "^(?:tmt(?:pro)?|itraq)(?:[0-9]+plex)?$",
    token,
    perl = TRUE
  )
  common_name || numeric_shift || shorthand || unimod_id || unimod_name ||
    isotope_label || reporter_label
}

strip_modification_groups <- function(value) {
  chars <- strsplit(value, "", fixed = TRUE)[[1]]
  if (length(chars) == 0L) {
    return(list(value = value, removed = character(), unsupported = character()))
  }

  openers <- c("[", "(", "{")
  closers <- c("]", ")", "}")
  matching_closer <- setNames(closers, openers)
  # Preallocate token buffers to avoid repeated copying while parsing.
  output <- character(length(chars))
  removed <- character(length(chars))
  unsupported <- character(max(1L, 2L * length(chars)))
  output_count <- 0L
  removed_count <- 0L
  unsupported_count <- 0L
  i <- 1L
  n_chars <- length(chars)

  while (i <= n_chars) {
    current <- chars[[i]]
    if (current == "*") {
      unsupported_count <- unsupported_count + 1L
      unsupported[[unsupported_count]] <- "*"
      i <- i + 1L
      next
    }
    if (!current %in% openers) {
      if (current %in% closers) {
        unsupported_count <- unsupported_count + 1L
        unsupported[[unsupported_count]] <- current
      } else {
        output_count <- output_count + 1L
        output[[output_count]] <- current
      }
      i <- i + 1L
      next
    }

    stack <- character(n_chars - i + 1L)
    stack_depth <- 1L
    stack[[stack_depth]] <- current
    j <- i + 1L
    while (j <= n_chars && stack_depth > 0L) {
      if (chars[[j]] %in% openers) {
        stack_depth <- stack_depth + 1L
        stack[[stack_depth]] <- chars[[j]]
      } else if (chars[[j]] %in% closers) {
        expected <- matching_closer[[stack[[stack_depth]]]]
        if (!identical(chars[[j]], expected)) {
          unsupported_count <- unsupported_count + 1L
          unsupported[[unsupported_count]] <- paste(chars[i:j], collapse = "")
          stack_depth <- 0L
          j <- n_chars + 1L
          break
        }
        stack_depth <- stack_depth - 1L
      }
      j <- j + 1L
    }
    if (stack_depth > 0L) {
      unsupported_count <- unsupported_count + 1L
      unsupported[[unsupported_count]] <- paste(chars[i:n_chars], collapse = "")
      break
    }

    token_end <- j - 1L
    token_indices <- if (token_end - i <= 1L) {
      integer()
    } else {
      (i + 1L):(token_end - 1L)
    }
    token <- if (length(token_indices) == 0L) {
      ""
    } else {
      paste(chars[token_indices], collapse = "")
    }
    if (supported_modification_token(token)) {
      removed_count <- removed_count + 1L
      removed[[removed_count]] <- token
    } else {
      unsupported_count <- unsupported_count + 1L
      unsupported[[unsupported_count]] <- token
    }
    i <- j
  }

  list(
    value = if (output_count == 0L) {
      ""
    } else {
      paste(output[seq_len(output_count)], collapse = "")
    },
    removed = if (removed_count == 0L) {
      character()
    } else {
      removed[seq_len(removed_count)]
    },
    unsupported = if (unsupported_count == 0L) {
      character()
    } else {
      unsupported[seq_len(unsupported_count)]
    }
  )
}

normalize_one_peptide_fields <- function(value, il_equivalent = TRUE) {
  if (length(value) != 1L || is.na(value)) {
    return(list(
      peptide_input = as.character(value %||% NA_character_),
      normalized_peptide = NA_character_,
      matching_peptide = NA_character_,
      normalization_status = "missing",
      normalization_notes = "Missing peptide value."
    ))
  }

  raw <- trimws(as.character(value))
  if (!nzchar(raw)) {
    return(list(
      peptide_input = raw,
      normalized_peptide = NA_character_,
      matching_peptide = NA_character_,
      normalization_status = "blank",
      normalization_notes = "Blank peptide value."
    ))
  }

  notes <- character()
  x <- gsub("\\s+", "", raw, perl = TRUE)
  if (!identical(x, raw)) {
    notes <- c(notes, "whitespace removed")
  }

  if (startsWith(x, "_") || endsWith(x, "_")) {
    x <- sub("^_", "", x)
    x <- sub("_$", "", x)
    notes <- c(notes, "terminal underscores removed")
  }

  flank_match <- regexec("^([A-Za-z_-])\\.(.*)\\.([A-Za-z_-])$", x, perl = TRUE)
  flank_parts <- regmatches(x, flank_match)[[1]]
  if (length(flank_parts) == 4L && nzchar(flank_parts[[3]])) {
    x <- flank_parts[[3]]
    notes <- c(notes, "flanking residues removed")
  }

  stripped <- if (grepl("[*\\[\\](){}]", x, perl = TRUE)) {
    strip_modification_groups(x)
  } else {
    list(value = x, removed = character(), unsupported = character())
  }
  if (length(stripped$removed) > 0L) {
    notes <- c(notes, "supported modification annotations removed")
  }
  if (length(stripped$unsupported) > 0L) {
    return(list(
      peptide_input = raw,
      normalized_peptide = NA_character_,
      matching_peptide = NA_character_,
      normalization_status = "unsupported_annotation",
      normalization_notes = paste(
        c(notes, paste0("unsupported annotation(s): ", paste(stripped$unsupported, collapse = ", "))),
        collapse = "; "
      )
    ))
  }

  x <- toupper(stripped$value)
  if (!nzchar(x)) {
    return(list(
      peptide_input = raw,
      normalized_peptide = NA_character_,
      matching_peptide = NA_character_,
      normalization_status = "empty_after_normalization",
      normalization_notes = paste(
        c(notes, "no amino-acid characters remain"),
        collapse = "; "
      )
    ))
  }

  chars <- strsplit(x, "", fixed = TRUE)[[1]]
  invalid <- sort(unique(chars[!chars %in% standard_amino_acids]))
  if (length(invalid) > 0L) {
    return(list(
      peptide_input = raw,
      normalized_peptide = NA_character_,
      matching_peptide = NA_character_,
      normalization_status = "invalid_characters",
      normalization_notes = paste(
        c(notes, paste0("unsupported character(s): ", paste(invalid, collapse = ", "))),
        collapse = "; "
      )
    ))
  }

  matching <- if (isTRUE(il_equivalent)) chartr("I", "L", x) else x
  if (isTRUE(il_equivalent) && !identical(matching, x)) {
    notes <- c(notes, "I/L-equivalent matching applied")
  }

  list(
    peptide_input = raw,
    normalized_peptide = x,
    matching_peptide = matching,
    normalization_status = "ok",
    normalization_notes = if (length(notes) == 0L) {
      "none"
    } else {
      paste(notes, collapse = "; ")
    }
  )
}

normalize_one_peptide <- function(value, il_equivalent = TRUE) {
  validate_logical_option(il_equivalent, "il_equivalent")
  as.data.frame(
    normalize_one_peptide_fields(value, il_equivalent = il_equivalent),
    stringsAsFactors = FALSE
  )
}

normalize_peptides <- function(values, il_equivalent = TRUE) {
  validate_logical_option(il_equivalent, "il_equivalent")
  if (length(values) == 0L) {
    return(data.frame(
      input_row_id = integer(),
      peptide_input = character(),
      normalized_peptide = character(),
      matching_peptide = character(),
      normalization_status = character(),
      normalization_notes = character(),
      stringsAsFactors = FALSE
    ))
  }
  # Normalize each raw value once. Large PSM tables commonly repeat the same
  # peptide sequence; caching those results avoids repeating the parser and
  # modification-group scan for every duplicate row.
  n_values <- length(values)
  raw_values <- as.character(values)
  unique_values <- unique(raw_values)
  unique_results <- lapply(unique_values, normalize_one_peptide_fields,
                           il_equivalent = il_equivalent)
  cache_index <- match(raw_values, unique_values)
  peptide_input <- vapply(
    unique_results[cache_index],
    `[[`,
    character(1),
    "peptide_input"
  )
  normalized_peptide <- vapply(
    unique_results[cache_index],
    `[[`,
    character(1),
    "normalized_peptide"
  )
  matching_peptide <- vapply(
    unique_results[cache_index],
    `[[`,
    character(1),
    "matching_peptide"
  )
  normalization_status <- vapply(
    unique_results[cache_index],
    `[[`,
    character(1),
    "normalization_status"
  )
  normalization_notes <- vapply(
    unique_results[cache_index],
    `[[`,
    character(1),
    "normalization_notes"
  )
  normalized <- data.frame(
    input_row_id = seq_len(n_values),
    peptide_input = peptide_input,
    normalized_peptide = normalized_peptide,
    matching_peptide = matching_peptide,
    normalization_status = normalization_status,
    normalization_notes = normalization_notes,
    stringsAsFactors = FALSE
  )
  rownames(normalized) <- NULL
  normalized
}

read_peptide_table <- function(path,
                               peptide_column = "peptide",
                               il_equivalent = TRUE) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stopf("Peptide CSV file does not exist: %s", path %||% "<missing>")
  }
  input <- tryCatch(
    utils::read.csv(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c(""),
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(error) {
      stopf("Could not read peptide CSV '%s': %s", path, conditionMessage(error))
    }
  )
  if (ncol(input) == 0L) {
    stopf("Peptide CSV has no columns: %s", path)
  }
  if (is.null(peptide_column) || length(peptide_column) != 1L || !nzchar(peptide_column)) {
    stopf("peptide_column must name one column in the peptide CSV.")
  }

  exact <- which(names(input) == peptide_column)
  case_insensitive <- which(tolower(names(input)) == tolower(peptide_column))
  selected <- if (length(exact) == 1L) exact else case_insensitive
  if (length(selected) != 1L) {
    stopf(
      "Peptide column '%s' was not found uniquely. Available columns: %s",
      peptide_column,
      paste(names(input), collapse = ", ")
    )
  }
  selected_name <- names(input)[[selected]]
  normalized <- normalize_peptides(input[[selected]], il_equivalent = il_equivalent)
  path <- normalizePath(path, mustWork = TRUE)

  list(
    input = input,
    peptide_column = selected_name,
    normalized = normalized,
    source_path = path,
    source_md5 = unname(as.character(tools::md5sum(path)))
  )
}

read_fusion_junctions <- function(path,
                                  sequence_text,
                                  sequence_ids = names(sequence_text)) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) {
    return(list(
      table = empty_junction_table(),
      source_path = NA_character_,
      source_md5 = NA_character_
    ))
  }
  if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
    stopf("Fusion-junction CSV file does not exist: %s", path %||% "<missing>")
  }
  if (is.null(names(sequence_text)) || any(!nzchar(names(sequence_text)))) {
    stopf("sequence_text must be a named character vector before reading junction metadata.")
  }

  # Read every column as text: "NA" is a valid identifier or inserted
  # sequence (Asn-Ala), and type conversion would rewrite IDs such as "01".
  # Only blank cells are missing; the optional coordinates also accept "NA".
  input <- tryCatch(
    utils::read.csv(
      path,
      check.names = FALSE,
      colClasses = "character",
      na.strings = "",
      fileEncoding = "UTF-8-BOM"
    ),
    error = function(error) {
      stopf("Could not read fusion-junction CSV '%s': %s", path, conditionMessage(error))
    }
  )
  required <- c(
    "fusion_id",
    "junction_id",
    "upstream_parent",
    "downstream_parent",
    "fusion_left_position",
    "fusion_right_position",
    "min_flank_aa"
  )
  require_data_frame_columns(input, required, "Fusion-junction CSV")
  if (nrow(input) == 0L) {
    stopf("Fusion-junction CSV contains no junction rows: %s", path)
  }

  for (column in c("source", "notes", "inserted_sequence")) {
    if (!column %in% names(input)) {
      input[[column]] <- ""
    }
  }
  text_columns <- c(
    "fusion_id",
    "junction_id",
    "upstream_parent",
    "downstream_parent"
  )
  for (column in text_columns) {
    values <- trimws(as.character(input[[column]]))
    if (any(is.na(values) | !nzchar(values))) {
      stopf("Fusion-junction column '%s' contains blank values.", column)
    }
    input[[column]] <- values
  }
  if (anyDuplicated(input$junction_id)) {
    stopf("Fusion-junction IDs must be unique.")
  }

  parse_integer_column <- function(column, allow_na = FALSE) {
    raw <- trimws(as.character(input[[column]]))
    if (allow_na) {
      raw[!is.na(raw) & raw == "NA"] <- NA_character_
    }
    numeric_values <- suppressWarnings(as.numeric(raw))
    not_integer <- !is.finite(numeric_values) |
      abs(numeric_values) > .Machine$integer.max |
      numeric_values != floor(numeric_values)
    invalid <- if (allow_na) {
      !is.na(raw) & nzchar(raw) & not_integer
    } else {
      not_integer
    }
    if (any(invalid)) {
      stopf("Fusion-junction column '%s' must contain integer values.", column)
    }
    as.integer(numeric_values)
  }

  input$fusion_left_position <- parse_integer_column("fusion_left_position")
  input$fusion_right_position <- parse_integer_column("fusion_right_position")
  input$min_flank_aa <- parse_integer_column("min_flank_aa")
  input$inserted_sequence <- toupper(trimws(as.character(input$inserted_sequence)))
  input$inserted_sequence[is.na(input$inserted_sequence)] <- ""
  invalid_insertions <- vapply(
    input$inserted_sequence,
    function(value) {
      if (!nzchar(value)) return(FALSE)
      chars <- strsplit(value, "", fixed = TRUE)[[1]]
      any(!chars %in% standard_amino_acids)
    },
    logical(1)
  )
  if (any(invalid_insertions)) {
    stopf("Fusion-junction inserted_sequence values must contain amino-acid symbols.")
  }
  for (column in c("upstream_parent_position", "downstream_parent_position")) {
    if (column %in% names(input)) {
      input[[column]] <- parse_integer_column(column, allow_na = TRUE)
    } else {
      input[[column]] <- NA_integer_
    }
  }

  if (any(!input$fusion_id %in% sequence_ids)) {
    stopf("Fusion-junction fusion_id values must match configured sequence roles.")
  }
  if (any(!input$upstream_parent %in% sequence_ids) ||
      any(!input$downstream_parent %in% sequence_ids)) {
    stopf("Fusion-junction parent values must match configured sequence roles.")
  }
  if (any(input$fusion_id != "Fusion")) {
    stopf("Fusion-junction fusion_id values must be 'Fusion'.")
  }
  parent_roles <- setdiff(sequence_ids, "Fusion")
  if (any(!input$upstream_parent %in% parent_roles) ||
      any(!input$downstream_parent %in% parent_roles)) {
    stopf("Fusion-junction parent values must refer to parent sequences, not Fusion.")
  }
  if (any(input$upstream_parent == input$downstream_parent)) {
    stopf("A fusion junction must have different upstream and downstream parents.")
  }
  if (any(input$fusion_left_position < 1L) ||
      any(input$fusion_right_position <= input$fusion_left_position)) {
    stopf(
      "Fusion junction coordinates must define an ordered boundary: right > left."
    )
  }
  for (row in seq_len(nrow(input))) {
    fusion_id <- input$fusion_id[[row]]
    upstream <- input$upstream_parent[[row]]
    downstream <- input$downstream_parent[[row]]
    left <- input$fusion_left_position[[row]]
    right <- input$fusion_right_position[[row]]
    if (right > nchar(sequence_text[[fusion_id]])) {
      stopf("Fusion junction '%s' lies outside the fusion sequence.", input$junction_id[[row]])
    }
    expected_inserted <- if (right - left > 1L) {
      substr(sequence_text[[fusion_id]], left + 1L, right - 1L)
    } else {
      ""
    }
    if (nchar(input$inserted_sequence[[row]]) != right - left - 1L ||
        !identical(input$inserted_sequence[[row]], expected_inserted)) {
      stopf(
        "Fusion junction '%s' inserted_sequence does not match the fusion sequence.",
        input$junction_id[[row]]
      )
    }
    parent_positions <- c(
      input$upstream_parent_position[[row]],
      input$downstream_parent_position[[row]]
    )
    if (xor(is.na(parent_positions[[1]]), is.na(parent_positions[[2]]))) {
      stopf(
        "Fusion junction '%s' must provide both parent coordinates or neither.",
        input$junction_id[[row]]
      )
    }
    if (!anyNA(parent_positions)) {
      upstream_position <- parent_positions[[1]]
      downstream_position <- parent_positions[[2]]
      if (upstream_position < 1L ||
          upstream_position > nchar(sequence_text[[upstream]]) ||
          downstream_position < 1L ||
          downstream_position > nchar(sequence_text[[downstream]])) {
        stopf(
          "Fusion junction '%s' contains a parent coordinate outside its parent sequence.",
          input$junction_id[[row]]
        )
      }
      fusion_left <- substr(sequence_text[[fusion_id]], left, left)
      fusion_right <- substr(sequence_text[[fusion_id]], right, right)
      parent_left <- substr(sequence_text[[upstream]], upstream_position, upstream_position)
      parent_right <- substr(sequence_text[[downstream]], downstream_position, downstream_position)
      if (!identical(fusion_left, parent_left) || !identical(fusion_right, parent_right)) {
        stopf(
          "Fusion junction '%s' does not agree with the supplied parent coordinates.",
          input$junction_id[[row]]
        )
      }
    }
  }
  if (any(input$min_flank_aa < 1L)) {
    stopf("min_flank_aa must be at least 1 for every fusion junction.")
  }

  input <- input[, c(
    "fusion_id",
    "junction_id",
    "upstream_parent",
    "downstream_parent",
    "fusion_left_position",
    "fusion_right_position",
    "upstream_parent_position",
    "downstream_parent_position",
    "inserted_sequence",
    "min_flank_aa",
    "source",
    "notes"
  )]
  rownames(input) <- NULL
  path <- normalizePath(path, mustWork = TRUE)
  list(
    table = input,
    source_path = path,
    source_md5 = unname(as.character(tools::md5sum(path)))
  )
}

# Peptide matching, junction annotation, and coverage

merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) {
    return(empty_interval_table())
  }
  starts <- as.integer(starts)
  ends <- as.integer(ends)
  valid <- !is.na(starts) & !is.na(ends) & starts <= ends
  starts <- starts[valid]
  ends <- ends[valid]
  if (length(starts) == 0L) {
    return(empty_interval_table())
  }

  order_index <- order(starts, ends)
  starts <- starts[order_index]
  ends <- ends[order_index]
  merged_start <- integer(length(starts))
  merged_end <- integer(length(starts))
  merged_count <- 1L

  current_start <- starts[[1]]
  current_end <- ends[[1]]
  if (length(starts) > 1L) {
    for (i in seq.int(2L, length(starts))) {
      if (starts[[i]] <= current_end + 1L) {
        current_end <- max(current_end, ends[[i]])
      } else {
        merged_start[[merged_count]] <- current_start
        merged_end[[merged_count]] <- current_end
        merged_count <- merged_count + 1L
        current_start <- starts[[i]]
        current_end <- ends[[i]]
      }
    }
  }
  merged_start[[merged_count]] <- current_start
  merged_end[[merged_count]] <- current_end

  data.frame(
    start = merged_start[seq_len(merged_count)],
    end = merged_end[seq_len(merged_count)],
    length = merged_end[seq_len(merged_count)] -
      merged_start[seq_len(merged_count)] + 1L,
    stringsAsFactors = FALSE
  )
}

complement_intervals <- function(merged, sequence_length) {
  if (sequence_length < 1L) {
    return(empty_interval_table())
  }
  if (nrow(merged) == 0L) {
    return(data.frame(
      start = 1L,
      end = as.integer(sequence_length),
      length = as.integer(sequence_length),
      stringsAsFactors = FALSE
    ))
  }

  starts <- integer(nrow(merged) + 1L)
  ends <- integer(nrow(merged) + 1L)
  interval_count <- 0L
  cursor <- 1L
  for (i in seq_len(nrow(merged))) {
    if (merged$start[[i]] > cursor) {
      interval_count <- interval_count + 1L
      starts[[interval_count]] <- cursor
      ends[[interval_count]] <- merged$start[[i]] - 1L
    }
    cursor <- max(cursor, merged$end[[i]] + 1L)
  }
  if (cursor <= sequence_length) {
    interval_count <- interval_count + 1L
    starts[[interval_count]] <- cursor
    ends[[interval_count]] <- as.integer(sequence_length)
  }

  if (interval_count == 0L) {
    return(empty_interval_table())
  }
  starts <- starts[seq_len(interval_count)]
  ends <- ends[seq_len(interval_count)]
  data.frame(
    start = starts,
    end = ends,
    length = ends - starts + 1L,
    stringsAsFactors = FALSE
  )
}

prepare_peptide_matcher <- function(normalized_peptides,
                                    il_equivalent = TRUE) {
  if (nrow(normalized_peptides) == 0L) {
    return(NULL)
  }
  valid <- normalized_peptides$normalization_status == "ok" &
    !is.na(normalized_peptides$matching_peptide)
  if (!any(valid)) {
    return(NULL)
  }

  keys <- unique(normalized_peptides$matching_peptide[valid])
  valid_rows <- which(valid)
  key_ids <- match(normalized_peptides$matching_peptide[valid_rows], keys)
  rows_by_key <- split(
    valid_rows,
    factor(key_ids, levels = seq_along(keys))
  )
  key_widths <- nchar(keys)
  width_groups <- split(
    seq_along(keys),
    factor(key_widths, levels = sort(unique(key_widths)))
  )
  dictionary_groups <- unname(width_groups[
    vapply(width_groups, length, integer(1)) >= 2L
  ])
  direct_groups <- unname(width_groups[
    vapply(width_groups, length, integer(1)) < 2L
  ])
  direct_patterns <- lapply(
    direct_groups,
    function(group) Biostrings::AAString(keys[[group[[1L]]]])
  )
  dictionaries <- vector("list", length(dictionary_groups))
  if (length(dictionary_groups) > 0L) {
    for (dictionary_index in seq_along(dictionary_groups)) {
      group <- dictionary_groups[[dictionary_index]]
      encoded_patterns <- Biostrings::DNAStringSet(vapply(
        keys[group],
        encode_protein_for_dictionary,
        character(1)
      ))
      # ACtree2 allocates a large fixed-size automaton even for tiny
      # dictionaries.  Twobit is sufficient for exact matching here and
      # keeps the compiled dictionary small.  A short trusted band is enough
      # to index candidates; matchPDict still verifies the complete pattern.
      trusted_band_width <- min(
        6L,
        3L * nchar(keys[group[[1L]]])
      )
      dictionaries[[dictionary_index]] <- Biostrings::PDict(
        encoded_patterns,
        algorithm = "Twobit",
        tb.start = 1L,
        tb.width = trusted_band_width
      )
    }
  }

  list(
    keys = keys,
    rows_by_key = rows_by_key,
    dictionary_groups = dictionary_groups,
    dictionaries = dictionaries,
    direct_groups = direct_groups,
    direct_patterns = direct_patterns,
    il_equivalent = isTRUE(il_equivalent)
  )
}

find_peptide_hits <- function(sequence_text,
                              sequence_id,
                              normalized_peptides,
                              il_equivalent = TRUE,
                              matcher = NULL) {
  if (is.null(matcher)) {
    matcher <- prepare_peptide_matcher(
      normalized_peptides,
      il_equivalent = il_equivalent
    )
  } else if (!identical(
    isTRUE(matcher$il_equivalent),
    isTRUE(il_equivalent)
  )) {
    stopf("The peptide matcher was prepared with a different I/L mode.")
  }
  if (is.null(matcher)) {
    return(empty_hit_table())
  }

  raw_sequence <- toupper(as.character(sequence_text))
  matching_sequence <- if (isTRUE(il_equivalent)) {
    # raw_sequence is already upper case, so avoid a second toupper() pass.
    chartr("I", "L", raw_sequence)
  } else {
    raw_sequence
  }
  keys <- matcher$keys
  rows_by_key <- matcher$rows_by_key
  starts_by_key <- vector("list", length(keys))
  ends_by_key <- vector("list", length(keys))

  # matchPDict searches a large dictionary in one pass, but only supports
  # equal-width patterns.  Grouping by peptide width lets us retain exact
  # variable-length peptide semantics while avoiding one full subject scan
  # per peptide when many peptides share a width.
  if (length(matcher$dictionary_groups) > 0L) {
    encoded_subject <- Biostrings::DNAString(
      encode_protein_for_dictionary(matching_sequence)
    )
    for (dictionary_index in seq_along(matcher$dictionary_groups)) {
      group <- matcher$dictionary_groups[[dictionary_index]]
      dictionary_matches <- Biostrings::matchPDict(
        matcher$dictionaries[[dictionary_index]],
        encoded_subject
      )
      for (local_index in seq_along(group)) {
        key_index <- group[[local_index]]
        encoded_views <- dictionary_matches[[local_index]]
        if (length(encoded_views) == 0L) {
          next
        }
        encoded_starts <- as.integer(IRanges::start(encoded_views))
        encoded_ends <- as.integer(IRanges::end(encoded_views))
        # The nucleotide dictionary can report matches beginning inside an
        # encoded codon.  Only codon-aligned words represent amino-acid
        # peptides.
        codon_aligned <- (encoded_starts - 1L) %% 3L == 0L
        encoded_starts <- encoded_starts[codon_aligned]
        encoded_ends <- encoded_ends[codon_aligned]
        if (length(encoded_starts) == 0L) {
          next
        }
        starts_by_key[[key_index]] <-
          (encoded_starts - 1L) %/% 3L + 1L
        ends_by_key[[key_index]] <- encoded_ends %/% 3L
      }
      rm(dictionary_matches)
    }
  }

  # Use direct AA matching for groups with fewer than two peptides. This
  # avoids constructing a dictionary for a single peptide.
  if (length(matcher$direct_groups) > 0L) {
    matching_subject <- Biostrings::AAString(matching_sequence)
    for (direct_index in seq_along(matcher$direct_groups)) {
      group <- matcher$direct_groups[[direct_index]]
      key_index <- group[[1L]]
      views <- Biostrings::matchPattern(
        matcher$direct_patterns[[direct_index]],
        matching_subject,
        fixed = TRUE
      )
      if (length(views) > 0L) {
        starts_by_key[[key_index]] <- as.integer(IRanges::start(views))
        ends_by_key[[key_index]] <- as.integer(IRanges::end(views))
      }
    }
  }

  hits <- vector("list", length(keys))
  hit_index <- 0L

  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    rows <- rows_by_key[[key_index]]
    peptide_input_values <- unique(normalized_peptides$peptide_input[rows])
    canonical_values <- unique(normalized_peptides$normalized_peptide[rows])
    starts <- starts_by_key[[key_index]]
    ends <- ends_by_key[[key_index]]
    if (is.null(starts) || length(starts) == 0L) {
      next
    }
    n_matches <- length(starts)
    hit_index <- hit_index + 1L
    hits[[hit_index]] <- data.frame(
      sequence_id = rep.int(sequence_id, n_matches),
      matching_peptide = rep.int(key, n_matches),
      normalized_peptide = rep.int(
        paste(canonical_values, collapse = " | "),
        n_matches
      ),
      peptide_input_values = rep.int(
        paste(peptide_input_values, collapse = " | "),
        n_matches
      ),
      input_row_ids = rep.int(
        paste(normalized_peptides$input_row_id[rows], collapse = ","),
        n_matches
      ),
      input_count = rep.int(length(rows), n_matches),
      start = starts,
      end = ends,
      length = ends - starts + 1L,
      matched_subsequence = substring(raw_sequence, starts, ends),
      match_basis = rep.int(
        if (isTRUE(il_equivalent)) "I/L-equivalent" else "exact",
        n_matches
      ),
      stringsAsFactors = FALSE
    )
  }

  if (hit_index == 0L) {
    empty_hit_table()
  } else {
    do.call(rbind, hits[seq_len(hit_index)])
  }
}

evaluate_peptide_junctions <- function(hits, junctions) {
  require_data_frame_columns(
    hits,
    c("sequence_id", "matching_peptide", "start", "end"),
    "hits"
  )
  if (!is.null(junctions) && !is.data.frame(junctions)) {
    stopf("junctions must be a data frame or NULL.")
  }
  if (!is.null(junctions) && nrow(junctions) > 0L) {
    require_data_frame_columns(
      junctions,
      c(
        "fusion_id",
        "junction_id",
        "fusion_left_position",
        "fusion_right_position",
        "min_flank_aa"
      ),
      "junctions"
    )
  }
  if (nrow(hits) == 0L || is.null(junctions) || nrow(junctions) == 0L) {
    return(empty_junction_evaluations())
  }

  evaluations <- vector("list", nrow(junctions))
  evaluation_count <- 0L
  relevant_junctions <- junctions[
    junctions$fusion_id %in% unique(hits$sequence_id),
    ,
    drop = FALSE
  ]
  for (junction_index in seq_len(nrow(relevant_junctions))) {
    junction <- relevant_junctions[junction_index, , drop = FALSE]
    candidate <- hits$sequence_id == junction$fusion_id[[1L]] &
      junction$fusion_left_position[[1L]] <= hits$end &
      junction$fusion_right_position[[1L]] >= hits$start &
      hits$start <= junction$fusion_left_position[[1L]] &
      hits$end >= junction$fusion_right_position[[1L]]
    hit_indices <- which(candidate)
    if (length(hit_indices) == 0L) {
      next
    }
    left_flank <- junction$fusion_left_position[[1L]] -
      hits$start[hit_indices] + 1L
    right_flank <- hits$end[hit_indices] -
      junction$fusion_right_position[[1L]] + 1L
    evaluation_count <- evaluation_count + 1L
    evaluations[[evaluation_count]] <- data.frame(
      hit_row_id = hit_indices,
      sequence_id = hits$sequence_id[hit_indices],
      matching_peptide = hits$matching_peptide[hit_indices],
      junction_id = rep.int(
        junction$junction_id[[1L]],
        length(hit_indices)
      ),
      junction_left_flank = left_flank,
      junction_right_flank = right_flank,
      min_flank_aa = rep.int(
        junction$min_flank_aa[[1L]],
        length(hit_indices)
      ),
      junction_flank_pass = left_flank >= junction$min_flank_aa[[1L]] &
        right_flank >= junction$min_flank_aa[[1L]],
      stringsAsFactors = FALSE
    )
  }
  if (evaluation_count == 0L) {
    empty_junction_evaluations()
  } else {
    do.call(rbind, evaluations[seq_len(evaluation_count)])
  }
}

apply_junction_evaluations <- function(hits,
                                       evaluations,
                                       junction_assessed = TRUE) {
  validate_logical_option(junction_assessed, "junction_assessed")
  require_data_frame_columns(
    hits,
    c(
      "sequence_id",
      "matching_peptide",
      "start",
      "end"
    ),
    "hits"
  )
  if (!is.data.frame(evaluations)) {
    stopf("evaluations must be a data frame.")
  }
  if (nrow(evaluations) > 0L) {
    require_data_frame_columns(
      evaluations,
      c(
        "hit_row_id",
        "sequence_id",
        "matching_peptide",
        "junction_id",
        "junction_left_flank",
        "junction_right_flank",
        "min_flank_aa",
        "junction_flank_pass"
      ),
      "evaluations"
    )
  }
  if (nrow(hits) == 0L) {
    return(hits)
  }
  hits$crosses_junction <- FALSE
  hits$junction_assessed <- rep.int(junction_assessed, nrow(hits))
  hits$junction_flank_pass <- FALSE
  hits$junction_id <- ""
  hits$junction_left_flank <- NA_integer_
  hits$junction_right_flank <- NA_integer_
  hits$junction_min_flank_aa <- NA_integer_
  if (nrow(evaluations) == 0L) {
    return(hits)
  }

  by_hit <- split(seq_len(nrow(evaluations)), evaluations$hit_row_id)
  for (hit_id in names(by_hit)) {
    hit_row <- as.integer(hit_id)
    rows <- by_hit[[hit_id]]
    current <- evaluations[rows, , drop = FALSE]
    hits$crosses_junction[[hit_row]] <- TRUE
    hits$junction_flank_pass[[hit_row]] <- any(current$junction_flank_pass)
    scores <- pmin(
      current$junction_left_flank,
      current$junction_right_flank
    )
    best <- order(
      !current$junction_flank_pass,
      -scores,
      current$junction_id
    )[[1L]]
    hits$junction_id[[hit_row]] <- current$junction_id[[best]]
    hits$junction_left_flank[[hit_row]] <-
      current$junction_left_flank[[best]]
    hits$junction_right_flank[[hit_row]] <-
      current$junction_right_flank[[best]]
    hits$junction_min_flank_aa[[hit_row]] <-
      current$min_flank_aa[[best]]
  }
  hits
}

annotate_peptide_hits_with_junctions <- function(hits, junctions) {
  apply_junction_evaluations(
    hits,
    evaluate_peptide_junctions(hits, junctions),
    junction_assessed = !is.null(junctions) && nrow(junctions) > 0L
  )
}

classification_for_presence <- function(in_fusion, in_parent_a, in_parent_b) {
  key <- paste(as.integer(c(in_fusion, in_parent_a, in_parent_b)), collapse = "")
  switch(
    key,
    `100` = "fusion_only",
    `110` = "fusion_and_parentA",
    `101` = "fusion_and_parentB",
    `111` = "shared_all_three",
    `010` = "parentA_only",
    `001` = "parentB_only",
    `011` = "parental_only_both",
    `000` = "not_found",
    "unclassified"
  )
}

classify_peptide_hit_display_classes <- function(hits) {
  if (nrow(hits) == 0L) {
    return(character())
  }
  require_data_frame_columns(
    hits,
    c(
      "sequence_id",
      "matching_peptide",
      "crosses_junction",
      "junction_flank_pass"
    ),
    "hits"
  )
  parent_matching_peptides <- unique(
    hits$matching_peptide[hits$sequence_id != "Fusion"]
  )
  junction_assessed <- if ("junction_assessed" %in% names(hits)) {
    !is.na(hits$junction_assessed) & hits$junction_assessed
  } else {
    rep(FALSE, nrow(hits))
  }
  output <- rep("parent_reference_mapping", nrow(hits))
  fusion_rows <- which(hits$sequence_id == "Fusion")
  for (hit_row in fusion_rows) {
    shared_with_parent <- hits$matching_peptide[[hit_row]] %in%
      parent_matching_peptides
    output[[hit_row]] <- if (!junction_assessed[[hit_row]]) {
      "fusion_mapping_junction_unassessed"
    } else if (isTRUE(hits$junction_flank_pass[[hit_row]]) &&
               shared_with_parent) {
      "junction_spanning_shared_parent"
    } else if (isTRUE(hits$junction_flank_pass[[hit_row]])) {
      "junction_spanning_candidate"
    } else if (isTRUE(hits$crosses_junction[[hit_row]])) {
      "junction_crossing_below_flank_threshold"
    } else {
      "fusion_mapping_not_junction_spanning"
    }
  }
  output
}

select_best_fusion_evidence_class <- function(classes) {
  priority <- c(
    junction_spanning_candidate = 1L,
    junction_spanning_shared_parent = 2L,
    junction_crossing_below_flank_threshold = 3L,
    fusion_mapping_junction_unassessed = 4L,
    fusion_mapping_not_junction_spanning = 5L
  )
  classes <- classes[!is.na(classes) & nzchar(classes)]
  classes <- classes[classes %in% names(priority)]
  if (length(classes) == 0L) {
    return("not_in_fusion")
  }
  classes[[which.min(unname(priority[classes]))]]
}

build_peptide_summary <- function(normalized_peptides,
                                  hits,
                                  junction_assessed = TRUE) {
  if (!is.logical(junction_assessed) ||
      length(junction_assessed) != 1L ||
      is.na(junction_assessed)) {
    stopf("junction_assessed must be one non-missing TRUE/FALSE value.")
  }
  valid <- normalized_peptides$normalization_status == "ok" &
    !is.na(normalized_peptides$matching_peptide)
  keys <- unique(normalized_peptides$matching_peptide[valid])
  if (length(keys) == 0L) {
    return(data.frame(
      matching_peptide = character(),
      normalized_peptide = character(),
      peptide_input_values = character(),
      input_row_ids = character(),
      input_count = integer(),
      fusion_hits = integer(),
      parentA_hits = integer(),
      parentB_hits = integer(),
      found_in_fusion = logical(),
      found_in_parentA = logical(),
      found_in_parentB = logical(),
      presence_classification = character(),
      crosses_junction = logical(),
      junction_assessed = logical(),
      junction_flank_pass = logical(),
      junction_id = character(),
      junction_left_flank = integer(),
      junction_right_flank = integer(),
      junction_min_flank_aa = integer(),
      junction_absent_from_supplied_parents = logical(),
      fusion_evidence_class = character(),
      stringsAsFactors = FALSE
    ))
  }

  valid_rows <- which(valid)
  key_ids <- match(normalized_peptides$matching_peptide[valid_rows], keys)
  rows_by_key <- split(
    valid_rows,
    factor(key_ids, levels = seq_along(keys))
  )
  hit_rows_by_key <- if (nrow(hits) == 0L) {
    list()
  } else {
    split(seq_len(nrow(hits)), hits$matching_peptide)
  }

  n_keys <- length(keys)
  output <- list(
    matching_peptide = keys,
    normalized_peptide = character(n_keys),
    peptide_input_values = character(n_keys),
    input_row_ids = character(n_keys),
    input_count = integer(n_keys),
    fusion_hits = integer(n_keys),
    parentA_hits = integer(n_keys),
    parentB_hits = integer(n_keys),
    found_in_fusion = logical(n_keys),
    found_in_parentA = logical(n_keys),
    found_in_parentB = logical(n_keys),
    presence_classification = character(n_keys),
    crosses_junction = logical(n_keys),
    junction_assessed = logical(n_keys),
    junction_flank_pass = logical(n_keys),
    junction_id = character(n_keys),
    junction_left_flank = rep(NA_integer_, n_keys),
    junction_right_flank = rep(NA_integer_, n_keys),
    junction_min_flank_aa = rep(NA_integer_, n_keys),
    junction_absent_from_supplied_parents = logical(n_keys),
    fusion_evidence_class = character(n_keys)
  )

  has_junction_columns <- all(c(
    "crosses_junction",
    "junction_flank_pass",
    "junction_id",
    "junction_left_flank",
    "junction_right_flank"
  ) %in% names(hits))
  has_junction_threshold <- "junction_min_flank_aa" %in% names(hits)

  for (key_index in seq_along(keys)) {
    key <- keys[[key_index]]
    rows <- rows_by_key[[key_index]]
    hit_rows <- if (nrow(hits) == 0L) {
      integer()
    } else {
      hit_rows_by_key[[key]] %||% integer()
    }
    hit_sequence_ids <- if (length(hit_rows) == 0L) {
      character()
    } else {
      hits$sequence_id[hit_rows]
    }
    fusion_hits <- sum(hit_sequence_ids == "Fusion")
    parent_a_hits <- sum(hit_sequence_ids == "ParentA")
    parent_b_hits <- sum(hit_sequence_ids == "ParentB")
    fusion_hit_rows <- if (length(hit_rows) == 0L) {
      integer()
    } else {
      hit_rows[hit_sequence_ids == "Fusion"]
    }

    if (length(fusion_hit_rows) > 0L && has_junction_columns) {
      crosses_junction <- any(hits$crosses_junction[fusion_hit_rows])
      junction_flank_pass <- any(hits$junction_flank_pass[fusion_hit_rows])
      flank_rows <- fusion_hit_rows[
        !is.na(hits$junction_left_flank[fusion_hit_rows]) &
          !is.na(hits$junction_right_flank[fusion_hit_rows])
      ]
      if (length(flank_rows) > 0L) {
        flank_scores <- pmin(
          hits$junction_left_flank[flank_rows],
          hits$junction_right_flank[flank_rows]
        )
        best_flank_order <- order(
          !hits$junction_flank_pass[flank_rows],
          -flank_scores,
          hits$junction_id[flank_rows]
        )
        best_flank_row <- flank_rows[[best_flank_order[[1L]]]]
        junction_id <- hits$junction_id[[best_flank_row]]
        best_left_flank <- hits$junction_left_flank[[best_flank_row]]
        best_right_flank <- hits$junction_right_flank[[best_flank_row]]
        best_min_flank <- if (has_junction_threshold) {
          hits$junction_min_flank_aa[[best_flank_row]]
        } else {
          NA_integer_
        }
      } else {
        junction_id <- ""
        best_left_flank <- NA_integer_
        best_right_flank <- NA_integer_
        best_min_flank <- NA_integer_
      }
    } else {
      crosses_junction <- FALSE
      junction_flank_pass <- FALSE
      junction_id <- ""
      best_left_flank <- NA_integer_
      best_right_flank <- NA_integer_
      best_min_flank <- NA_integer_
    }

    junction_unique <- junction_flank_pass &&
      parent_a_hits == 0L &&
      parent_b_hits == 0L
    fusion_evidence <- if ("display_class" %in% names(hits)) {
      if (!fusion_hits) {
        "not_in_fusion"
      } else {
        selected_class <- select_best_fusion_evidence_class(
          hits$display_class[fusion_hit_rows]
        )
        if (identical(selected_class, "junction_spanning_shared_parent")) {
          "junction_spanning_but_shared_with_supplied_parent"
        } else {
          selected_class
        }
      }
    } else {
      if (!fusion_hits) {
        "not_in_fusion"
      } else if (!junction_assessed) {
        "fusion_mapping_junction_unassessed"
      } else if (junction_unique) {
        "junction_spanning_candidate"
      } else if (junction_flank_pass) {
        "junction_spanning_but_shared_with_supplied_parent"
      } else if (crosses_junction) {
        "junction_crossing_below_flank_threshold"
      } else if (parent_a_hits > 0L || parent_b_hits > 0L) {
        "parent_derived_or_shared"
      } else {
        "fusion_mapping_not_junction_spanning"
      }
    }

    output$normalized_peptide[[key_index]] <- paste(
      unique(normalized_peptides$normalized_peptide[rows]),
      collapse = " | "
    )
    output$peptide_input_values[[key_index]] <- paste(
      unique(normalized_peptides$peptide_input[rows]),
      collapse = " | "
    )
    output$input_row_ids[[key_index]] <- paste(
      normalized_peptides$input_row_id[rows],
      collapse = ","
    )
    output$input_count[[key_index]] <- length(rows)
    output$fusion_hits[[key_index]] <- fusion_hits
    output$parentA_hits[[key_index]] <- parent_a_hits
    output$parentB_hits[[key_index]] <- parent_b_hits
    output$found_in_fusion[[key_index]] <- fusion_hits > 0L
    output$found_in_parentA[[key_index]] <- parent_a_hits > 0L
    output$found_in_parentB[[key_index]] <- parent_b_hits > 0L
    output$presence_classification[[key_index]] <- classification_for_presence(
      fusion_hits > 0L,
      parent_a_hits > 0L,
      parent_b_hits > 0L
    )
    output$crosses_junction[[key_index]] <- crosses_junction
    output$junction_assessed[[key_index]] <- junction_assessed
    output$junction_flank_pass[[key_index]] <- junction_flank_pass
    output$junction_id[[key_index]] <- junction_id
    output$junction_left_flank[[key_index]] <- best_left_flank
    output$junction_right_flank[[key_index]] <- best_right_flank
    output$junction_min_flank_aa[[key_index]] <- best_min_flank
    output$junction_absent_from_supplied_parents[[key_index]] <-
      if (junction_assessed) junction_unique else NA
    output$fusion_evidence_class[[key_index]] <- fusion_evidence
  }

  as.data.frame(output, stringsAsFactors = FALSE)
}

calculate_coverage <- function(sequence_text, sequence_ids, hits) {
  coverage_rows <- vector("list", length(sequence_ids))
  uncovered_rows <- vector("list", length(sequence_ids))
  hit_rows_by_sequence <- if (nrow(hits) == 0L) {
    list()
  } else {
    split(seq_len(nrow(hits)), hits$sequence_id)
  }

  for (i in seq_along(sequence_ids)) {
    sequence_id <- sequence_ids[[i]]
    sequence_length <- nchar(sequence_text[[sequence_id]])
    sequence_hit_rows <- if (nrow(hits) == 0L) {
      integer()
    } else {
      hit_rows_by_sequence[[sequence_id]] %||% integer()
    }
    merged <- if (length(sequence_hit_rows) == 0L) {
      merge_intervals(integer(), integer())
    } else {
      merge_intervals(
        hits$start[sequence_hit_rows],
        hits$end[sequence_hit_rows]
      )
    }
    uncovered <- complement_intervals(merged, sequence_length)

    if (nrow(uncovered) > 0L) {
      uncovered_rows[[i]] <- data.frame(
        sequence_id = rep.int(sequence_id, nrow(uncovered)),
        start = uncovered$start,
        end = uncovered$end,
        length = uncovered$length,
        stringsAsFactors = FALSE
      )
    }

    covered <- if (nrow(merged) == 0L) 0L else sum(merged$length)
    coverage_rows[[i]] <- data.frame(
      sequence_id = sequence_id,
      sequence_length = sequence_length,
      n_unique_peptides = if (length(sequence_hit_rows) == 0L) {
        0L
      } else {
        length(unique(hits$matching_peptide[sequence_hit_rows]))
      },
      n_hit_occurrences = length(sequence_hit_rows),
      covered_residues = covered,
      uncovered_residues = sequence_length - covered,
      coverage_percent = round(100 * covered / sequence_length, 2),
      stringsAsFactors = FALSE
    )
  }

  list(
    summary = do.call(rbind, coverage_rows),
    uncovered_regions = if (all(vapply(uncovered_rows, is.null, logical(1)))) {
      empty_uncovered_regions()
    } else {
      do.call(
        rbind,
        uncovered_rows[!vapply(uncovered_rows, is.null, logical(1))]
      )
    }
  )
}

# Gap-aware fusion/parent alignment

make_position_vector_from_chars <- function(chars, start_position = 1L) {
  if (length(chars) == 0L) {
    return(integer())
  }
  positions <- as.integer(cumsum(chars != "-") + as.integer(start_position) - 1L)
  positions[chars == "-"] <- NA_integer_
  positions
}

alignment_start <- function(alignment, side, default = 1L) {
  aligned_object <- if (identical(side, "pattern")) {
    pwalign::pattern(alignment)
  } else {
    pwalign::subject(alignment)
  }
  value <- tryCatch(
    as.integer(IRanges::start(aligned_object)[[1]]),
    error = function(error) NA_integer_
  )
  if (is.na(value)) as.integer(default) else value
}

alignment_score <- function(alignment) {
  tryCatch(
    as.numeric(pwalign::score(alignment)[[1]]),
    error = function(error) NA_real_
  )
}

make_fusion_parent_columns <- function(alignment,
                                       parent_id,
                                       alignment_type,
                                       pattern_role,
                                       pattern_start = 1L,
                                       subject_start = 1L) {
  aligned_pattern <- as.character(pwalign::alignedPattern(alignment))[[1]]
  aligned_subject <- as.character(pwalign::alignedSubject(alignment))[[1]]
  if (nchar(aligned_pattern) != nchar(aligned_subject)) {
    stopf("Alignment returned unequal pattern and subject widths for %s.", parent_id)
  }
  if (!nzchar(aligned_pattern)) {
    return(empty_alignment_columns())
  }

  pattern_chars <- strsplit(aligned_pattern, "", fixed = TRUE)[[1]]
  subject_chars <- strsplit(aligned_subject, "", fixed = TRUE)[[1]]
  pattern_positions <- make_position_vector_from_chars(
    pattern_chars,
    pattern_start
  )
  subject_positions <- make_position_vector_from_chars(
    subject_chars,
    subject_start
  )

  if (identical(pattern_role, "fusion")) {
    fusion_chars <- pattern_chars
    parent_chars <- subject_chars
    fusion_positions <- pattern_positions
    parent_positions <- subject_positions
  } else {
    fusion_chars <- subject_chars
    parent_chars <- pattern_chars
    fusion_positions <- subject_positions
    parent_positions <- pattern_positions
  }

  raw_match <- fusion_chars != "-" &
    parent_chars != "-" &
    fusion_chars == parent_chars
  il_match <- fusion_chars != "-" &
    parent_chars != "-" &
    chartr("I", "L", fusion_chars) == chartr("I", "L", parent_chars)
  status <- rep("double_gap", length(fusion_chars))
  parent_only <- fusion_chars == "-" & parent_chars != "-"
  fusion_only <- fusion_chars != "-" & parent_chars == "-"
  mismatch <- fusion_chars != "-" & parent_chars != "-" & !raw_match
  status[parent_only] <- "parent_only"
  status[fusion_only] <- "fusion_only"
  status[raw_match] <- "match"
  status[mismatch] <- "mismatch"

  data.frame(
    parent = parent_id,
    alignment_type = alignment_type,
    alignment_column = seq_along(fusion_chars),
    fusion_pos = as.integer(fusion_positions),
    parent_pos = as.integer(parent_positions),
    fusion_aa = fusion_chars,
    parent_aa = parent_chars,
    raw_match = raw_match,
    il_equivalent_match = il_match,
    status = status,
    stringsAsFactors = FALSE
  )
}

alignment_regions_from_columns <- function(columns) {
  if (nrow(columns) == 0L) {
    return(empty_alignment_regions())
  }
  starts_new <- c(TRUE, columns$status[-1L] != columns$status[-nrow(columns)])
  region_starts <- which(starts_new)
  region_ends <- c(region_starts[-1L] - 1L, nrow(columns))
  n_regions <- length(region_starts)
  fusion_start <- rep(NA_integer_, n_regions)
  fusion_end <- rep(NA_integer_, n_regions)
  parent_start <- rep(NA_integer_, n_regions)
  parent_end <- rep(NA_integer_, n_regions)

  for (region_index in seq_len(n_regions)) {
    index <- region_starts[[region_index]]:region_ends[[region_index]]
    fusion_positions <- columns$fusion_pos[index]
    fusion_positions <- fusion_positions[!is.na(fusion_positions)]
    if (length(fusion_positions) > 0L) {
      fusion_start[[region_index]] <- min(fusion_positions)
      fusion_end[[region_index]] <- max(fusion_positions)
    }
    parent_positions <- columns$parent_pos[index]
    parent_positions <- parent_positions[!is.na(parent_positions)]
    if (length(parent_positions) > 0L) {
      parent_start[[region_index]] <- min(parent_positions)
      parent_end[[region_index]] <- max(parent_positions)
    }
  }

  data.frame(
    parent = columns$parent[region_starts],
    alignment_type = columns$alignment_type[region_starts],
    region_id = seq_len(n_regions),
    status = columns$status[region_starts],
    fusion_start = fusion_start,
    fusion_end = fusion_end,
    parent_start = parent_start,
    parent_end = parent_end,
    alignment_columns = region_ends - region_starts + 1L,
    stringsAsFactors = FALSE
  )
}

alignment_summary_row <- function(columns,
                                  parent_id,
                                  alignment_type,
                                  fusion_length,
                                  parent_length,
                                  score,
                                  substitution_matrix) {
  if (nrow(columns) == 0L) {
    return(data.frame(
      parent = parent_id,
      alignment_type = alignment_type,
      substitution_matrix = substitution_matrix,
      score = score,
      alignment_columns = 0L,
      aligned_residues = 0L,
      raw_matches = 0L,
      il_equivalent_matches = 0L,
      raw_identity_percent = 0,
      il_equivalent_identity_percent = 0,
      fusion_aligned_fraction = 0,
      parent_aligned_fraction = 0,
      stringsAsFactors = FALSE
    ))
  }
  aligned <- !is.na(columns$fusion_pos) & !is.na(columns$parent_pos)
  n_aligned <- sum(aligned)
  raw_matches <- sum(columns$raw_match, na.rm = TRUE)
  il_matches <- sum(columns$il_equivalent_match, na.rm = TRUE)
  data.frame(
    parent = parent_id,
    alignment_type = alignment_type,
    substitution_matrix = substitution_matrix,
    score = score,
    alignment_columns = nrow(columns),
    aligned_residues = n_aligned,
    raw_matches = raw_matches,
    il_equivalent_matches = il_matches,
    raw_identity_percent = if (n_aligned == 0L) 0 else round(100 * raw_matches / n_aligned, 2),
    il_equivalent_identity_percent = if (n_aligned == 0L) {
      0
    } else {
      round(100 * il_matches / n_aligned, 2)
    },
    fusion_aligned_fraction = round(n_aligned / fusion_length, 4),
    parent_aligned_fraction = round(n_aligned / parent_length, 4),
    stringsAsFactors = FALSE
  )
}

simple_amino_acid_matrix <- function() {
  alphabet <- LETTERS
  matrix_value <- matrix(
    -1,
    nrow = length(alphabet),
    ncol = length(alphabet),
    dimnames = list(alphabet, alphabet)
  )
  diag(matrix_value) <- 2
  matrix_value
}

resolve_substitution_matrix <- function(fusion_sequence,
                                        parent_sequence,
                                        requested) {
  if (is.matrix(requested)) {
    return(list(value = requested, label = "custom"))
  }
  if (!is.character(requested) || length(requested) != 1L) {
    stopf("substitution_matrix must be a matrix or one matrix name.")
  }
  if (!identical(toupper(requested), "BLOSUM62")) {
    return(list(value = requested, label = requested))
  }

  blosum_alphabet <- strsplit("ACDEFGHIKLMNPQRSTVWYBXJZ", "", fixed = TRUE)[[1]]
  observed <- unique(strsplit(
    paste0(toupper(fusion_sequence), toupper(parent_sequence)),
    "",
    fixed = TRUE
  )[[1]])
  if (all(observed %in% blosum_alphabet)) {
    list(value = "BLOSUM62", label = "BLOSUM62")
  } else {
    list(
      value = simple_amino_acid_matrix(),
      label = "simple_aa_compatible"
    )
  }
}

run_pairwise_alignment <- function(fusion_sequence,
                                   parent_sequence,
                                   parent_id,
                                   gap_opening = 10,
                                   gap_extension = 0.5,
                                   substitution_matrix = "BLOSUM62",
                                   store_alignment_objects = TRUE) {
  if (length(fusion_sequence) != 1L || length(parent_sequence) != 1L) {
    stopf("Pairwise alignment inputs must each contain exactly one sequence.")
  }
  validate_nonnegative_number(gap_opening, "gap_opening")
  validate_nonnegative_number(gap_extension, "gap_extension")
  validate_logical_option(store_alignment_objects, "store_alignment_objects")
  fusion_sequence <- toupper(as.character(fusion_sequence))
  parent_sequence <- toupper(as.character(parent_sequence))
  fusion_length <- nchar(fusion_sequence)
  parent_length <- nchar(parent_sequence)
  matrix_spec <- resolve_substitution_matrix(
    fusion_sequence,
    parent_sequence,
    substitution_matrix
  )

  global_alignment <- tryCatch(
    pwalign::pairwiseAlignment(
      pattern = Biostrings::AAString(fusion_sequence),
      subject = Biostrings::AAString(parent_sequence),
      type = "global",
      substitutionMatrix = matrix_spec$value,
      gapOpening = gap_opening,
      gapExtension = gap_extension
    ),
    error = function(error) {
      stopf("Global alignment failed for Fusion vs %s: %s", parent_id, conditionMessage(error))
    }
  )
  global_columns <- make_fusion_parent_columns(
    alignment = global_alignment,
    parent_id = parent_id,
    alignment_type = "global",
    pattern_role = "fusion"
  )
  global_summary <- alignment_summary_row(
    global_columns,
    parent_id,
    "global",
    fusion_length,
    parent_length,
    alignment_score(global_alignment),
    matrix_spec$label
  )
  if (!isTRUE(store_alignment_objects)) {
    global_alignment <- NULL
  }

  local_alignment <- tryCatch(
    pwalign::pairwiseAlignment(
      pattern = Biostrings::AAString(parent_sequence),
      subject = Biostrings::AAString(fusion_sequence),
      type = "local",
      substitutionMatrix = matrix_spec$value,
      gapOpening = gap_opening,
      gapExtension = gap_extension
    ),
    error = function(error) {
      stopf("Local alignment failed for Fusion vs %s: %s", parent_id, conditionMessage(error))
    }
  )
  local_columns <- make_fusion_parent_columns(
    alignment = local_alignment,
    parent_id = parent_id,
    alignment_type = "local",
    pattern_role = "parent",
    pattern_start = alignment_start(local_alignment, "pattern"),
    subject_start = alignment_start(local_alignment, "subject")
  )
  local_summary <- alignment_summary_row(
    local_columns,
    parent_id,
    "local",
    fusion_length,
    parent_length,
    alignment_score(local_alignment),
    matrix_spec$label
  )

  columns <- rbind(global_columns, local_columns)
  regions <- rbind(
    alignment_regions_from_columns(global_columns),
    alignment_regions_from_columns(local_columns)
  )
  summaries <- rbind(global_summary, local_summary)
  if (!isTRUE(store_alignment_objects)) {
    local_alignment <- NULL
  }

  list(
    parent = parent_id,
    global = if (isTRUE(store_alignment_objects)) global_alignment else NULL,
    local = if (isTRUE(store_alignment_objects)) local_alignment else NULL,
    columns = columns,
    regions = regions,
    summaries = summaries,
    substitution_matrix_used = matrix_spec$label
  )
}

alignment_overview_bins <- function(alignment_columns, bin_width = 25L) {
  if (!is.numeric(bin_width) || length(bin_width) != 1L ||
      is.na(bin_width) || !is.finite(bin_width) ||
      bin_width < 1L || bin_width != floor(bin_width)) {
    stopf("bin_width must be one finite positive integer.")
  }
  global_index <- alignment_columns$alignment_type == "global" &
    !is.na(alignment_columns$fusion_pos)
  local_index <- alignment_columns$alignment_type == "local" &
    !is.na(alignment_columns$fusion_pos)
  global_fusion_positions <- alignment_columns$fusion_pos[global_index]
  local <- alignment_columns[
    local_index,
    c(
      "parent",
      "fusion_pos",
      "parent_pos",
      "raw_match",
      "il_equivalent_match"
    ),
    drop = FALSE
  ]
  fusion_positions <- c(global_fusion_positions, local$fusion_pos)
  fusion_positions <- fusion_positions[!is.na(fusion_positions)]
  if (length(fusion_positions) == 0L) {
    return(data.frame(
      parent = character(),
      bin_start = integer(),
      bin_end = integer(),
      bin_mid = numeric(),
      n_aligned = integer(),
      fusion_coverage_percent = numeric(),
      raw_identity_percent = numeric(),
      il_equivalent_identity_percent = numeric(),
      status = character(),
      stringsAsFactors = FALSE
    ))
  }
  fusion_length <- max(fusion_positions)
  bin_width <- as.integer(bin_width)
  bins <- seq.int(1L, fusion_length, by = bin_width)
  bin_ends <- pmin(fusion_length, bins + bin_width - 1L)
  bin_lengths <- bin_ends - bins + 1L
  parents <- unique(alignment_columns$parent)
  output <- list()
  output_index <- 0L
  for (parent in parents) {
    parent_local <- local[
      local$parent == parent &
        !is.na(local$fusion_pos) &
        !is.na(local$parent_pos),
      ,
      drop = FALSE
    ]
    if (nrow(parent_local) == 0L) {
      n_aligned <- integer(length(bins))
      raw_matches <- numeric(length(bins))
      il_matches <- numeric(length(bins))
    } else {
      bin_ids <- pmin(
        length(bins),
        ((parent_local$fusion_pos - 1L) %/% bin_width) + 1L
      )
      n_aligned <- tabulate(bin_ids, nbins = length(bins))
      raw_matches <- numeric(length(bins))
      raw_sums <- rowsum(
        as.numeric(parent_local$raw_match),
        bin_ids,
        reorder = FALSE
      )
      raw_matches[as.integer(rownames(raw_sums))] <- raw_sums[, 1L]
      il_matches <- numeric(length(bins))
      il_sums <- rowsum(
        as.numeric(parent_local$il_equivalent_match),
        bin_ids,
        reorder = FALSE
      )
      il_matches[as.integer(rownames(il_sums))] <- il_sums[, 1L]
    }
    coverage <- 100 * n_aligned / bin_lengths
    raw_identity <- ifelse(
      n_aligned == 0L,
      NA_real_,
      100 * raw_matches / n_aligned
    )
    il_identity <- ifelse(
      n_aligned == 0L,
      NA_real_,
      100 * il_matches / n_aligned
    )
    status <- rep("aligned_with_differences", length(bins))
    status[n_aligned == 0L] <- "no_local_alignment"
    partial <- n_aligned > 0L & coverage < 50
    status[partial] <- "partial_local_alignment"
    high_exact <- n_aligned > 0L &
      !partial &
      raw_identity >= 95
    status[high_exact] <- "high_exact_identity"
    high_il <- n_aligned > 0L &
      !partial &
      !high_exact &
      il_identity >= 95
    status[high_il] <- "high_il_identity"

    output_index <- output_index + 1L
    output[[output_index]] <- data.frame(
      parent = rep.int(parent, length(bins)),
      bin_start = bins,
      bin_end = bin_ends,
      bin_mid = (bins + bin_ends) / 2,
      n_aligned = as.integer(n_aligned),
      fusion_coverage_percent = round(coverage, 1),
      raw_identity_percent = round(raw_identity, 1),
      il_equivalent_identity_percent = round(il_identity, 1),
      status = status,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

# Provenance and reusable pipeline layers

format_interval_list <- function(intervals) {
  if (nrow(intervals) == 0L) return("")
  paste(sprintf("%d-%d", intervals$start, intervals$end), collapse = "; ")
}

package_version_text <- function(package) {
  version <- tryCatch(
    utils::packageDescription(package, fields = "Version"),
    error = function(error) NA_character_
  )
  if (length(version) == 0L || is.na(version) || !nzchar(version)) {
    paste0(package, "=unavailable")
  } else {
    paste0(package, "=", version)
  }
}

# The manifest identifies inputs by file name and MD5 hash. It omits the
# directories so shared outputs do not expose local paths and runs of the same
# inputs on different machines produce the same manifest.
build_run_manifest <- function(config,
                               sequence_input,
                               peptide_input,
                               junction_input = NULL,
                               warnings = character()) {
  keys <- c(
    "generated_at_utc",
    "r_version",
    "package_versions",
    "fasta_file",
    "fasta_md5",
    "peptide_csv_file",
    "peptide_csv_md5",
    "peptide_column",
    "junction_csv_file",
    "junction_csv_md5",
    "junction_count",
    "il_equivalent_matching",
    "gap_opening",
    "gap_extension",
    "substitution_matrix",
    "alignment_matrix_used",
    "store_alignment_objects",
    "sequence_ids",
    "warnings"
  )
  values <- list(
    format(Sys.time(), tz = "UTC", usetz = TRUE),
    as.character(getRversion()),
    paste(
      vapply(
        c("Biostrings", "pwalign", "ggplot2", "htmltools", "renv"),
        package_version_text,
        character(1)
      ),
      collapse = ", "
    ),
    basename(sequence_input$source_path),
    sequence_input$source_md5,
    basename(peptide_input$source_path),
    peptide_input$source_md5,
    peptide_input$peptide_column,
    if (is.null(junction_input)) NA_character_ else basename(junction_input$source_path),
    if (is.null(junction_input)) NA_character_ else junction_input$source_md5,
    if (is.null(junction_input)) 0L else nrow(junction_input$table),
    isTRUE(config$il_equivalent),
    config$gap_opening,
    config$gap_extension,
    config$substitution_matrix,
    config$alignment_matrix_used,
    isTRUE(config$store_alignment_objects),
    paste(names(config$sequence_ids), config$sequence_ids, sep = "="),
    if (length(warnings) == 0L) "none" else warnings
  )
  data.frame(
    key = keys,
    value = vapply(values, function(value) {
      paste(value, collapse = ", ")
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

map_peptides_to_references <- function(sequence_text,
                                       normalized_peptides,
                                       junctions = empty_junction_table(),
                                       il_equivalent = TRUE) {
  validate_logical_option(il_equivalent, "il_equivalent")
  if (!is.character(sequence_text) ||
      is.null(names(sequence_text)) ||
      anyNA(names(sequence_text)) ||
      any(!nzchar(names(sequence_text)))) {
    stopf("sequence_text must be a named character vector.")
  }
  require_data_frame_columns(
    normalized_peptides,
    c(
      "input_row_id",
      "peptide_input",
      "normalized_peptide",
      "matching_peptide",
      "normalization_status"
    ),
    "normalized_peptides"
  )
  if (is.null(junctions)) {
    junctions <- empty_junction_table()
  }
  if (!is.data.frame(junctions)) {
    stopf("junctions must be a data frame or NULL.")
  }

  peptide_matcher <- prepare_peptide_matcher(
    normalized_peptides,
    il_equivalent = il_equivalent
  )
  hit_tables <- lapply(names(sequence_text), function(sequence_id) {
    find_peptide_hits(
      sequence_text[[sequence_id]],
      sequence_id,
      normalized_peptides,
      il_equivalent = il_equivalent,
      matcher = peptide_matcher
    )
  })
  peptide_hits <- if (length(hit_tables) == 0L) {
    empty_hit_table()
  } else {
    do.call(rbind, hit_tables)
  }
  junction_evaluations <- evaluate_peptide_junctions(
    peptide_hits,
    junctions
  )
  peptide_hits <- apply_junction_evaluations(
    peptide_hits,
    junction_evaluations,
    junction_assessed = nrow(junctions) > 0L
  )
  peptide_hits$display_class <- classify_peptide_hit_display_classes(
    peptide_hits
  )
  rm(peptide_matcher, hit_tables)

  coverage <- calculate_coverage(
    sequence_text,
    names(sequence_text),
    peptide_hits
  )
  list(
    peptide_hits = peptide_hits,
    junction_evaluations = junction_evaluations,
    peptide_summary = build_peptide_summary(
      normalized_peptides,
      peptide_hits,
      junction_assessed = nrow(junctions) > 0L
    ),
    coverage = coverage$summary,
    uncovered_regions = coverage$uncovered_regions
  )
}

align_fusion_to_parents <- function(sequence_text,
                                    gap_opening = 10,
                                    gap_extension = 0.5,
                                    substitution_matrix = "BLOSUM62",
                                    store_alignment_objects = FALSE) {
  required_roles <- c("Fusion", "ParentA", "ParentB")
  if (!is.character(sequence_text) ||
      is.null(names(sequence_text)) ||
      anyNA(names(sequence_text)) ||
      !all(required_roles %in% names(sequence_text))) {
    stopf(
      "sequence_text must contain named Fusion, ParentA, and ParentB sequences."
    )
  }
  for (role in required_roles) {
    validate_sequence_text(sequence_text[[role]], role)
  }
  parent_ids <- c("ParentA", "ParentB")
  alignment_results <- vector("list", length(parent_ids))
  names(alignment_results) <- parent_ids
  columns_by_parent <- vector("list", length(parent_ids))
  regions_by_parent <- vector("list", length(parent_ids))
  summaries_by_parent <- vector("list", length(parent_ids))
  for (parent_index in seq_along(parent_ids)) {
    parent_id <- parent_ids[[parent_index]]
    alignment <- run_pairwise_alignment(
      sequence_text[["Fusion"]],
      sequence_text[[parent_id]],
      parent_id,
      gap_opening = gap_opening,
      gap_extension = gap_extension,
      substitution_matrix = substitution_matrix,
      store_alignment_objects = store_alignment_objects
    )
    # Retain only the optional raw objects for the parent being returned.
    # Flattened tables are collected separately and become the canonical
    # machine-readable outputs.
    alignment_results[[parent_index]] <- alignment[
      c("parent", "global", "local", "substitution_matrix_used")
    ]
    columns_by_parent[[parent_index]] <- alignment$columns
    regions_by_parent[[parent_index]] <- alignment$regions
    summaries_by_parent[[parent_index]] <- alignment$summaries
    rm(alignment)
  }

  alignment_columns <- do.call(rbind, columns_by_parent)
  alignment_regions <- do.call(rbind, regions_by_parent)
  alignment_summaries <- do.call(rbind, summaries_by_parent)
  list(
    alignments = alignment_results,
    alignment_columns = alignment_columns,
    alignment_regions = alignment_regions,
    alignment_summaries = alignment_summaries
  )
}

run_fusion_analysis <- function(sequence_file,
                                peptide_file,
                                sequence_ids = c(
                                  Fusion = "Fusion",
                                  ParentA = "ParentA",
                                  ParentB = "ParentB"
                                ),
                                peptide_column = "peptide",
                                il_equivalent = TRUE,
                                gap_opening = 10,
                                gap_extension = 0.5,
                                substitution_matrix = "BLOSUM62",
                                junction_file = NULL,
                                store_alignment_objects = FALSE) {
  sequence_input <- read_fusion_sequences(sequence_file, sequence_ids)
  peptide_input <- read_peptide_table(
    peptide_file,
    peptide_column = peptide_column,
    il_equivalent = il_equivalent
  )
  junction_input <- read_fusion_junctions(
    junction_file,
    sequence_text = sequence_input$sequence_text,
    sequence_ids = names(sequence_input$sequence_text)
  )

  normalized <- peptide_input$normalized
  peptide_input$normalized <- NULL
  mapping <- map_peptides_to_references(
    sequence_input$sequence_text,
    normalized,
    junctions = junction_input$table,
    il_equivalent = il_equivalent
  )

  alignment <- align_fusion_to_parents(
    sequence_input$sequence_text,
    gap_opening = gap_opening,
    gap_extension = gap_extension,
    substitution_matrix = substitution_matrix,
    store_alignment_objects = store_alignment_objects
  )

  warnings <- normalized$normalization_notes[
    normalized$normalization_status != "ok"
  ]
  warnings <- unique(warnings[nzchar(warnings)])
  config <- list(
    sequence_ids = sequence_ids,
    peptide_column = peptide_input$peptide_column,
    il_equivalent = il_equivalent,
    gap_opening = gap_opening,
    gap_extension = gap_extension,
    substitution_matrix = substitution_matrix,
    store_alignment_objects = store_alignment_objects,
    junction_file = junction_input$source_path,
    alignment_matrix_used = paste(
      unique(alignment$alignment_summaries$substitution_matrix),
      collapse = ", "
    )
  )

  result <- list(
    config = config,
    sequence_input = sequence_input,
    peptide_input = peptide_input,
    normalized_peptides = normalized,
    junctions = junction_input$table,
    peptide_hits = mapping$peptide_hits,
    junction_evaluations = mapping$junction_evaluations,
    peptide_summary = mapping$peptide_summary,
    coverage = mapping$coverage,
    uncovered_regions = mapping$uncovered_regions,
    alignments = alignment$alignments,
    alignment_columns = alignment$alignment_columns,
    alignment_regions = alignment$alignment_regions,
    alignment_summaries = alignment$alignment_summaries,
    warnings = warnings,
    manifest = build_run_manifest(
      config,
      sequence_input,
      peptide_input,
      junction_input,
      warnings
    )
  )
  class(result) <- c("fusion_peptide_mapping_result", "list")
  result
}
