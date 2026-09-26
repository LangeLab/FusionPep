# Static HTML report writer for FusionPep.

html_text <- function(value, attribute = FALSE) {
  if (length(value) == 0L || is.na(value)) {
    value <- ""
  }
  htmltools::htmlEscape(as.character(value), attribute = attribute)
}

friendly_table_headers <- function(values) {
  labels <- gsub("([a-z])([A-Z])", "\\1 \\2", values, perl = TRUE)
  labels <- gsub("_", " ", labels, fixed = TRUE)
  labels <- tools::toTitleCase(labels)
  labels <- gsub("\\bId\\b", "ID", labels, perl = TRUE)
  labels <- gsub("\\bAa\\b", "AA", labels, perl = TRUE)
  labels <- gsub("\\bIl\\b", "I/L", labels, perl = TRUE)
  labels <- gsub("\\bRds\\b", "RDS", labels, perl = TRUE)
  labels <- gsub("Parent a", "Parent A", labels, fixed = TRUE)
  labels <- gsub("Parent b", "Parent B", labels, fixed = TRUE)
  explicit <- c(
    sequence_id = "Reference sequence",
    sequence_length = "Sequence length (AA)",
    length = "Length (AA)",
    fusion_hits = "Fusion hits",
    parentA_hits = "Parent A hits",
    parentB_hits = "Parent B hits",
    coverage_percent = "Coverage (%)",
    minimum_flank_aa = "Minimum flank (AA)",
    junction_flanks = "Junction flanks (AA)",
    high_similarity_fusion_residues = "High-identity fusion bins (AA)",
    high_similarity_parent_residues = "Paired parent residues (AA)",
    overall_local_exact_identity = "Overall local exact identity (%)",
    fusion_residues_locally_aligned = "Fusion residues locally aligned (%)"
  )
  mapped <- match(values, names(explicit))
  labels[!is.na(mapped)] <- unname(explicit[mapped[!is.na(mapped)]])
  labels
}

table_to_html <- function(data, max_rows = 200L) {
  if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (!is.numeric(max_rows) || length(max_rows) != 1L ||
      is.na(max_rows) || !is.finite(max_rows) ||
      max_rows < 1L || max_rows != floor(max_rows) ||
      max_rows > .Machine$integer.max) {
    stop("max_rows must be one positive integer.", call. = FALSE)
  }
  if (nrow(data) == 0L) {
    return('<p class="empty">No rows.</p>')
  }
  display_rows <- seq_len(min(nrow(data), as.integer(max_rows)))
  header <- paste(
    sprintf(
      '<th scope="col">%s</th>',
      vapply(friendly_table_headers(names(data)), html_text, character(1))
    ),
    collapse = ""
  )
  display_data <- data[display_rows, , drop = FALSE]
  formatted_columns <- lapply(display_data, function(column) {
    vapply(seq_along(column), function(index) {
      value <- column[[index]]
      html_text(if (length(value) == 0L || is.na(value) ||
                    (is.character(value) && !nzchar(trimws(value)))) {
        "N/A"
      } else {
        as.character(value)
      })
    }, character(1))
  })
  cell_classes <- ifelse(
    vapply(display_data, is.numeric, logical(1)),
    "numeric",
    ifelse(
      names(display_data) %in% c(
        "peptide", "peptide_input", "matching_peptide",
        "normalized_peptide", "peptide_input_values"
      ),
      "peptide",
      ifelse(
        names(display_data) %in% c("parent", "sequence_id", "fusion_id"),
        "reference", "text"
      )
    )
  )
  rows <- vapply(seq_along(display_rows), function(row_index) {
    values <- vapply(
      formatted_columns,
      function(column) column[[row_index]],
      character(1)
    )
    paste0(
      "<tr>",
      paste(sprintf('<td class="%s">%s</td>', cell_classes, values), collapse = ""),
      "</tr>"
    )
  }, character(1))
  table_html <- paste0(
    '<div class="table-wrap"><table><thead><tr>',
    header,
    "</tr></thead><tbody>",
    paste(rows, collapse = ""),
    "</tbody></table></div>"
  )
  if (nrow(data) > length(display_rows)) {
    paste0(
      table_html,
      '<p class="table-truncated">Showing the first ',
      length(display_rows),
      " of ",
      nrow(data),
      ' rows in this report. Complete tables are described in the <a href="#files">output manifest</a>.</p>'
    )
  } else {
    table_html
  }
}

sequence_annotation_html <- function(sequence_text,
                                      hits,
                                      sequence_id,
                                      line_width = 80L) {
  sequence_text <- toupper(as.character(sequence_text))
  chars <- strsplit(sequence_text, "", fixed = TRUE)[[1]]
  class_priority <- c(
    "junction-candidate",
    "junction-shared",
    "junction-flank-fail",
    "parent-mapping",
    "fusion-mapping"
  )
  labels <- vector("list", length(chars))
  class_rank <- integer(length(chars))
  label_values <- character()
  display_to_rank <- c(
    junction_spanning_candidate = 1L,
    junction_spanning_shared_parent = 2L,
    junction_crossing_below_flank_threshold = 3L,
    parent_reference_mapping = 4L,
    fusion_mapping_not_junction_spanning = 5L,
    fusion_mapping_junction_unassessed = 5L
  )
  if (nrow(hits) > 0L) {
    display_classes <- if ("display_class" %in% names(hits)) {
      hits$display_class
    } else {
      classify_peptide_hit_display_classes(hits)
    }
    sequence_hit_rows <- which(hits$sequence_id == sequence_id)
    if (length(sequence_hit_rows) > 0L) {
      # Store small integer label IDs per residue rather than repeatedly
      # copying long peptide strings.  The rendered tooltip remains
      # unchanged.
      label_values <- unique(hits$peptide_input_values[sequence_hit_rows])
      for (hit_row in sequence_hit_rows) {
        interval <- max(1L, hits$start[[hit_row]]):min(
          length(chars),
          hits$end[[hit_row]]
        )
        label <- hits$peptide_input_values[[hit_row]]
        label_id <- match(label, label_values)
        hit_class_rank <- unname(
          display_to_rank[display_classes[[hit_row]]]
        )
        labels[interval] <- lapply(labels[interval], function(existing) {
          unique(c(existing, label_id))
        })
        previous_rank <- class_rank[interval]
        class_rank[interval] <- ifelse(
          previous_rank == 0L,
          hit_class_rank,
          pmin(previous_rank, hit_class_rank)
        )
      }
    }
  }

  lines <- lapply(seq(1L, length(chars), by = line_width), function(start) {
    end <- min(length(chars), start + line_width - 1L)
    positions <- start:end
    rendered <- chars[positions]
    mapped <- which(lengths(labels[positions]) > 0L)
    if (length(mapped) > 0L) {
      for (line_index in mapped) {
        position <- positions[[line_index]]
        peptide_labels <- label_values[labels[[position]]]
        if (length(peptide_labels) > 0L) {
          title <- html_text(
            paste("Mapped peptide(s):", paste(peptide_labels, collapse = ", ")),
            attribute = TRUE
          )
          residue_class <- class_priority[[class_rank[[position]]]]
          rendered[[line_index]] <- sprintf(
            '<strong class="mapped-residue %s" title="%s">%s</strong>',
            html_text(residue_class, attribute = TRUE),
            title,
            chars[[position]]
          )
        }
      }
    }
    sprintf(
      '<span class="sequence-line"><span class="position">%6d</span> %s <span class="position">%6d</span></span>',
      start,
      paste(
        vapply(
          split(rendered, ceiling(seq_along(rendered) / 10)),
          paste, character(1), collapse = ""
        ),
        collapse = " "
      ),
      end
    )
  })
  paste(lines, collapse = "<br>\n")
}

alignment_block_html <- function(columns,
                                  parent_id,
                                  alignment_type = "local",
                                  chunk_width = 80L) {
  selected <- columns[
    columns$parent == parent_id & columns$alignment_type == alignment_type,
    ,
    drop = FALSE
  ]
  if (nrow(selected) == 0L) {
    return(
      sprintf(
        '<p class="empty">No %s alignment columns.</p>',
        html_text(alignment_type)
      )
    )
  }
  label_width <- max(nchar(c("Fusion", parent_id)))
  position_width <- max(nchar(as.character(
    c(selected$fusion_pos, selected$parent_pos)
  )), na.rm = TRUE)
  chunks <- split(seq_len(nrow(selected)), ceiling(seq_len(nrow(selected)) / chunk_width))
  blocks <- lapply(chunks, function(index) {
    fusion <- paste(selected$fusion_aa[index], collapse = "")
    parent <- paste(selected$parent_aa[index], collapse = "")
    markers <- ifelse(
      selected$raw_match[index],
      "|",
      ifelse(selected$il_equivalent_match[index], ":", " ")
    )
    markers[selected$status[index] %in% c("parent_only", "fusion_only", "double_gap")] <- " "
    marker_line <- paste(markers, collapse = "")
    fusion_positions <- selected$fusion_pos[index]
    parent_positions <- selected$parent_pos[index]
    position_range <- function(values) {
      values <- values[!is.na(values)]
      if (length(values) == 0L) c("-", "-") else as.character(range(values))
    }
    fusion_range <- position_range(fusion_positions)
    parent_range <- position_range(parent_positions)
    paste(
      sprintf(
        "%-*s %*s %s %*s", label_width, "Fusion",
        position_width, fusion_range[[1L]], fusion,
        position_width, fusion_range[[2L]]
      ),
      paste0(strrep(" ", label_width + position_width + 2L), marker_line),
      sprintf(
        "%-*s %*s %s %*s", label_width, parent_id,
        position_width, parent_range[[1L]], parent,
        position_width, parent_range[[2L]]
      ),
      sep = "\n"
    )
  })
  paste0(
    '<p class="alignment-key">Coordinates are 1-based. | = exact match; : = I/L-equivalent match; - = gap.</p>',
    '<div class="alignment-blocks">',
    paste0(
      '<pre class="alignment">',
      vapply(blocks, html_text, character(1)),
      "</pre>",
      collapse = "\n"
    ),
    "</div>"
  )
}

embedded_image <- function(path, alt_text, caption = alt_text, pdf_path = NULL) {
  if (length(path) == 0L || is.na(path) || !file.exists(path)) {
    return(sprintf('<p class="warning">Missing figure: %s</p>', html_text(path)))
  }
  uri <- base64enc::dataURI(file = path, mime = "image/png")
  downloads <- paste0(
    '<a href="figures/', html_text(basename(path), attribute = TRUE),
    '">600 dpi PNG</a>'
  )
  if (!is.null(pdf_path) && file.exists(pdf_path)) {
    downloads <- paste0(
      '<a href="figures/', html_text(basename(pdf_path), attribute = TRUE),
      '">Vector PDF</a> <span aria-hidden="true">|</span> ', downloads
    )
  }
  sprintf(
    '<figure><img src="%s" alt="%s"><figcaption>%s</figcaption><p class="figure-downloads">%s</p></figure>',
    html_text(uri, attribute = TRUE),
    html_text(alt_text, attribute = TRUE),
    html_text(caption),
    downloads
  )
}

# Describe the saved data files using counts from the same result as the report.
output_manifest_html <- function(result, artifact_paths) {
  # nolint start: line_length_linter. One readable sentence per output file.
  descriptions <- c(
    input_peptides = "Original peptide rows, including all supplied evidence-role and other metadata columns.",
    normalized_peptides = "One row per input peptide: normalized sequence, matching sequence, acceptance status and normalization notes.",
    sequence_metadata = "Reference roles, original FASTA headers and sequence lengths.",
    fusion_junctions = "Supplied junctions after validation: fusion and parent coordinates, inserted residues and minimum flanks.",
    peptide_hits = "Every mapped occurrence: reference, inclusive coordinates, original input row IDs, matching mode and display class.",
    peptide_summary = "One row per distinct matching sequence: reference presence, hit counts and the selected junction decision.",
    junction_evaluations = "Every crossing occurrence evaluated against each applicable junction, with both flanks and its threshold decision.",
    coverage_summary = "Per-reference peptide counts, covered and uncovered residue counts, and coverage percentage.",
    uncovered_regions = "Inclusive start and end positions of each uncovered interval in each reference.",
    alignment_columns = "Every global and local alignment column: paired residues, coordinates, gaps, exact matches and I/L-equivalent matches.",
    alignment_regions = "Contiguous runs of alignment matches, mismatches and gaps, including all regions omitted from the report preview.",
    alignment_summaries = "Per-parent global and local alignment scores, identities, aligned fractions and the substitution matrix used.",
    run_manifest = "Generation time, R and package versions, input file names and MD5 hashes, matching mode and analysis settings.",
    warnings = "Normalization warnings for excluded or unsupported input values, or an explicit statement that there were none.",
    result_rds = "Complete analysis object, including inputs, settings and all result tables. Open with readRDS() in R."
  )
  # nolint end
  row_counts <- c(
    input_peptides = nrow(result$peptide_input$input),
    normalized_peptides = nrow(result$normalized_peptides),
    sequence_metadata = length(result$sequence_input$sequence_text),
    fusion_junctions = nrow(result$junctions),
    peptide_hits = nrow(result$peptide_hits),
    peptide_summary = nrow(result$peptide_summary),
    junction_evaluations = nrow(result$junction_evaluations),
    coverage_summary = nrow(result$coverage),
    uncovered_regions = nrow(result$uncovered_regions),
    alignment_columns = nrow(result$alignment_columns),
    alignment_regions = nrow(result$alignment_regions),
    alignment_summaries = nrow(result$alignment_summaries),
    run_manifest = nrow(result$manifest)
  )
  names_to_show <- names(descriptions)[names(descriptions) %in% names(artifact_paths)]
  if (length(names_to_show) == 0L) {
    return('<p class="empty">No saved data files were supplied to this report.</p>')
  }
  rows <- vapply(names_to_show, function(name) {
    path <- artifact_paths[[name]]
    link <- if (file.exists(path)) {
      sprintf(
        '<a href="%s"><code>%s</code></a>',
        html_text(basename(path), attribute = TRUE), html_text(basename(path))
      )
    } else {
      paste0("<code>", html_text(basename(path)), "</code><br>File not available")
    }
    count <- if (name %in% names(row_counts)) {
      as.character(row_counts[[name]])
    } else {
      "Not a table"
    }
    paste0(
      '<tr><td class="manifest-file">', link,
      '</td><td class="numeric">', html_text(count),
      "</td><td>", html_text(descriptions[[name]]), "</td></tr>"
    )
  }, character(1))
  paste0(
    "<p>These files contain the complete data behind this report. Row counts exclude CSV headers. ",
    "Figure PDF and PNG downloads remain beside their figures. ",
    "Keep the report and its output folder together so file links work.</p>",
    '<div class="table-wrap"><table class="output-manifest"><thead><tr>',
    '<th scope="col">File</th><th scope="col">Rows</th>',
    '<th scope="col">Contents and use</th></tr></thead><tbody>',
    paste(rows, collapse = ""), "</tbody></table></div>"
  )
}

input_roles_for_summary <- function(result, summary) {
  input <- result$peptide_input$input
  if (nrow(summary) == 0L || nrow(input) == 0L ||
      !"evidence_role" %in% names(input)) {
    return(rep("", nrow(summary)))
  }
  vapply(summary$input_row_ids, function(row_ids) {
    if (is.na(row_ids) || !nzchar(row_ids)) return("")
    indices <- suppressWarnings(as.integer(strsplit(row_ids, ",", fixed = TRUE)[[1]]))
    indices <- indices[!is.na(indices) & indices >= 1L & indices <= nrow(input)]
    if (length(indices) == 0L) return("")
    roles <- unique(trimws(as.character(input$evidence_role[indices])))
    roles <- roles[!is.na(roles) & nzchar(roles)]
    paste(roles, collapse = " | ")
  }, character(1))
}

friendly_presence_class <- function(values) {
  labels <- c(
    fusion_only = "Fusion only",
    fusion_and_parentA = "Fusion + ParentA",
    fusion_and_parentB = "Fusion + ParentB",
    shared_all_three = "Shared by all three",
    parentA_only = "ParentA only",
    parentB_only = "ParentB only",
    parental_only_both = "Both parents only",
    not_found = "Not found"
  )
  output <- unname(labels[as.character(values)])
  output[is.na(output)] <- as.character(values)[is.na(output)]
  output
}

friendly_evidence_class <- function(values) {
  labels <- c(
    junction_spanning_candidate = "Junction candidate",
    junction_crossing_below_flank_threshold = "Crosses junction; flank fail",
    junction_spanning_but_shared_with_supplied_parent = "Junction sequence also in parent",
    parent_derived_or_shared = "Parent-derived/shared",
    fusion_mapping_junction_unassessed = "Fusion mapping; junction not assessed",
    fusion_mapping_not_junction_spanning = "Fusion mapping; not junction-spanning",
    not_in_fusion = "Not in fusion"
  )
  output <- unname(labels[as.character(values)])
  output[is.na(output)] <- as.character(values)[is.na(output)]
  output
}

friendly_input_role <- function(values) {
  output <- gsub("_", " ", as.character(values), fixed = TRUE)
  output <- tools::toTitleCase(output)
  output[is.na(values) | !nzchar(trimws(as.character(values)))] <- "Unspecified"
  output
}

peptide_interpretation_table <- function(result) {
  summary <- result$peptide_summary
  if (nrow(summary) == 0L) {
    return(data.frame(
      peptide = character(),
      input_role = character(),
      presence = character(),
      fusion_hits = integer(),
      parentA_hits = integer(),
      parentB_hits = integer(),
      crosses_junction = character(),
      junction_flanks = character(),
      fusion_evidence = character(),
      stringsAsFactors = FALSE
    ))
  }
  threshold_for <- summary$junction_min_flank_aa
  flank_text <- ifelse(
    summary$crosses_junction,
    paste0(
      summary$junction_left_flank,
      " / ",
      summary$junction_right_flank,
      " (",
      ifelse(summary$junction_flank_pass, "pass", "fail"),
      ifelse(is.na(threshold_for), "", paste0("; min ", threshold_for)),
      ")"
    ),
    "N/A"
  )
  data.frame(
    peptide = summary$peptide_input_values,
    input_role = friendly_input_role(input_roles_for_summary(result, summary)),
    presence = friendly_presence_class(summary$presence_classification),
    fusion_hits = summary$fusion_hits,
    parentA_hits = summary$parentA_hits,
    parentB_hits = summary$parentB_hits,
    crosses_junction = ifelse(
      !summary$junction_assessed, "not assessed",
      ifelse(summary$crosses_junction, "yes", "no")
    ),
    junction_flanks = flank_text,
    fusion_evidence = friendly_evidence_class(summary$fusion_evidence_class),
    stringsAsFactors = FALSE
  )
}

junction_display_table <- function(junctions) {
  if (is.null(junctions) || nrow(junctions) == 0L) {
    return(data.frame(
      junction = character(),
      parent_join = character(),
      fusion_boundary = character(),
      inserted_residues = character(),
      minimum_flank_aa = integer(),
      source = character(),
      notes = character(),
      stringsAsFactors = FALSE
    ))
  }
  upstream_labels <- ifelse(
    is.na(junctions$upstream_parent_position),
    paste(junctions$upstream_parent, "(position not supplied)"),
    paste(junctions$upstream_parent, junctions$upstream_parent_position)
  )
  downstream_labels <- ifelse(
    is.na(junctions$downstream_parent_position),
    paste(junctions$downstream_parent, "(position not supplied)"),
    paste(junctions$downstream_parent, junctions$downstream_parent_position)
  )
  data.frame(
    junction = junctions$junction_id,
    parent_join = paste(upstream_labels, "to", downstream_labels),
    fusion_boundary = paste0(
      "residue ",
      junctions$fusion_left_position,
      ifelse(
        nzchar(junctions$inserted_sequence),
        paste0(" | ", junctions$inserted_sequence, " | "),
        " | "
      ),
      "residue ",
      junctions$fusion_right_position
    ),
    inserted_residues = ifelse(
      nzchar(junctions$inserted_sequence),
      junctions$inserted_sequence,
      "None"
    ),
    minimum_flank_aa = junctions$min_flank_aa,
    source = junctions$source,
    notes = junctions$notes,
    stringsAsFactors = FALSE
  )
}

summary_cards_html <- function(result) {
  summary <- result$peptide_summary
  valid <- sum(result$normalized_peptides$normalization_status == "ok")
  candidates <- sum(
    summary$fusion_evidence_class == "junction_spanning_candidate"
  )
  below_threshold <- sum(
    summary$fusion_evidence_class == "junction_crossing_below_flank_threshold"
  )
  parent_mapped <- sum(
    summary$found_in_fusion &
      (summary$found_in_parentA | summary$found_in_parentB)
  )
  junction_assessed <- nrow(result$junctions) > 0L
  fusion_coverage <- result$coverage$coverage_percent[
    result$coverage$sequence_id == "Fusion"
  ]
  candidate_value <- if (junction_assessed) candidates else "N/A"
  candidate_detail <- if (junction_assessed) {
    "distinct matching sequences; sequence-level only"
  } else {
    "not assessed; no fusion-junction metadata supplied"
  }
  flank_failure_value <- if (junction_assessed) below_threshold else "N/A"
  flank_failure_detail <- if (junction_assessed) {
    "distinct matching sequences; the flank rule failed"
  } else {
    "not assessed; no fusion-junction metadata supplied"
  }
  card <- function(value, label, detail) {
    paste0(
      '<div class="card"><div class="card-value">',
      html_text(value),
      '</div><div class="card-label">',
      html_text(label),
      '</div><div class="card-detail">',
      html_text(detail),
      "</div></div>"
    )
  }
  paste0(
    '<div class="cards">',
    card(valid, "valid input rows", "after normalization"),
    card(
      candidate_value,
      "junction candidates",
      candidate_detail
    ),
    card(
      flank_failure_value,
      "flank failures",
      flank_failure_detail
    ),
    card(
      parent_mapped,
      "fusion + parent matches",
      "distinct matching sequences shared with a supplied parent"
    ),
    card(
      if (length(fusion_coverage) == 0L) "N/A" else paste0(fusion_coverage, "%"),
      "fusion coverage",
      "union of mapped occurrences"
    ),
    "</div>"
  )
}

alignment_overview_table <- function(result) {
  bins <- alignment_overview_bins(result$alignment_columns)
  summaries <- result$alignment_summaries[
    result$alignment_summaries$alignment_type == "local",
    ,
    drop = FALSE
  ]
  if (nrow(summaries) == 0L) {
    return(data.frame(
      parent = character(),
      high_similarity_fusion_residues = character(),
      high_similarity_parent_residues = character(),
      overall_local_exact_identity = numeric(),
      fusion_residues_locally_aligned = numeric(),
      interpretation = character(),
      stringsAsFactors = FALSE
    ))
  }
  output <- lapply(seq_len(nrow(summaries)), function(i) {
    parent <- summaries$parent[[i]]
    parent_bins <- bins[bins$parent == parent, , drop = FALSE]
    high_bins <- parent_bins[
      parent_bins$status %in% c("high_exact_identity", "high_il_identity"),
      ,
      drop = FALSE
    ]
    high_ranges <- if (nrow(high_bins) == 0L) {
      ""
    } else {
      format_interval_list(
        merge_intervals(high_bins$bin_start, high_bins$bin_end)
      )
    }
    high_parent_positions <- unlist(lapply(
      seq_len(nrow(high_bins)),
      function(bin_index) {
        bin_positions <- result$alignment_columns$parent_pos[
          result$alignment_columns$parent == parent &
            result$alignment_columns$alignment_type == "local" &
            !is.na(result$alignment_columns$parent_pos) &
            !is.na(result$alignment_columns$fusion_pos) &
            result$alignment_columns$fusion_pos >=
              high_bins$bin_start[[bin_index]] &
            result$alignment_columns$fusion_pos <=
              high_bins$bin_end[[bin_index]]
        ]
        bin_positions
      }
    ), use.names = FALSE)
    high_parent_range <- format_interval_list(
      merge_intervals(high_parent_positions, high_parent_positions)
    )
    interpretation <- if (nzchar(high_ranges)) {
      paste0(
        "High-similarity block(s): Fusion ",
        high_ranges,
        "; inspect remaining bins separately"
      )
    } else {
      "No high-similarity fusion bin under the overview thresholds"
    }
    data.frame(
      parent = parent,
      high_similarity_fusion_residues = high_ranges,
      high_similarity_parent_residues = high_parent_range,
      overall_local_exact_identity = summaries$raw_identity_percent[[i]],
      fusion_residues_locally_aligned =
        round(100 * summaries$fusion_aligned_fraction[[i]], 1),
      interpretation = interpretation,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

write_fusion_report <- function(result, output_dir, artifact_paths) {
  if (!inherits(result, "fusion_peptide_mapping_result")) {
    stopf("write_fusion_report expects a fusion_peptide_mapping_result.")
  }
  output_dir <- normalize_output_directory(output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  report_paths <- fusion_output_paths(output_dir)
  report_path <- report_paths[["report"]]
  assert_output_paths_do_not_overwrite_inputs(result, report_paths["report"])

  references <- data.frame(
    sequence_id = names(result$sequence_input$sequence_text),
    length = nchar(result$sequence_input$sequence_text),
    fasta_header = unname(result$sequence_input$source_headers),
    stringsAsFactors = FALSE
  )
  interpretation_table <- peptide_interpretation_table(result)
  junction_table <- junction_display_table(result$junctions)
  junction_evaluation_table <- result$junction_evaluations %||%
    empty_junction_evaluations()
  alignment_table <- alignment_overview_table(result)
  junction_assessed <- nrow(result$junctions) > 0L
  junction_assessment_html <- if (!junction_assessed) {
    paste0(
      '<div class="warning"><b>Junction status is not assessed.</b> ',
      "No fusion-junction metadata was supplied, so fusion matches cannot ",
      "be classified as junction-spanning or non-junction-spanning.</div>"
    )
  } else {
    ""
  }
  junction_evaluation_html <- if (nrow(junction_evaluation_table) == 0L) {
    ""
  } else {
    paste0(
      "<details><summary>Junction evaluation audit</summary>",
      "<p>Each row evaluates one mapped occurrence that contains both boundary ",
      "residues against one supplied junction. Non-crossing occurrences are ",
      "not included. The selected junction shown in the peptide summary ",
      "is the passing evaluation with the strongest complete flank pair, or ",
      "the strongest failing pair when none passes. Flanks always come from ",
      "one occurrence and one junction.</p>",
      table_to_html(junction_evaluation_table),
      "</details>"
    )
  }
  candidate_regions <- result$alignment_regions[
    result$alignment_regions$status %in% c("parent_only", "fusion_only", "mismatch"),
    , drop = FALSE
  ]
  warning_html <- if (length(result$warnings) == 0L) {
    ""
  } else {
    paste0(
      '<div class="warning"><b>Input normalization warnings</b><ul>',
      paste(
        sprintf("<li>%s</li>", vapply(result$warnings, html_text, character(1))),
        collapse = ""
      ),
      "</ul></div>"
    )
  }
  sequence_sections <- vapply(names(result$sequence_input$sequence_text), function(id) {
    paste0(
      '<details class="sequence-detail">',
      "<summary>", html_text(id), " (",
      nchar(result$sequence_input$sequence_text[[id]]), " aa)</summary>",
      '<div class="sequence">',
      sequence_annotation_html(
        result$sequence_input$sequence_text[[id]], result$peptide_hits, id
      ),
      "</div></details>"
    )
  }, character(1))
  alignment_sections <- vapply(c("ParentA", "ParentB"), function(parent_id) {
    paste0(
      '<details class="alignment-detail"><summary>Detailed local alignment: Fusion vs ',
      html_text(parent_id), "</summary>",
      alignment_block_html(result$alignment_columns, parent_id, alignment_type = "local"),
      "</details>"
    )
  }, character(1))
  sequence_annotation_note <- if (!junction_assessed) {
    paste(
      "Blue marks parent-reference mappings and grey marks fusion mappings.",
      "Junction status was not assessed because no junction metadata was supplied."
    )
  } else {
    paste(
      "Blue marks parent-reference mappings, grey marks ordinary fusion mappings,",
      "orange marks junction crossings that fail the flank rule, pink marks",
      "junction-spanning matches also found in a supplied parent, and green",
      "marks junction candidates absent from both supplied parents."
    )
  }

  junction_figure_html <- embedded_image(
    artifact_paths[["junction_evidence_plot"]],
    "Fusion-junction peptide evidence map",
    paste(
      "Figure 1. Peptide matches across each supplied junction.",
      "The flank rule is evaluated independently in each panel; inserted residues",
      "do not contribute to either flank. Labels give the matched reference sequence,",
      "full inclusive coordinates and both flank lengths.",
      "Each panel states how many crossing occurrences are shown.",
      "Passing candidates are absent from the two supplied parents, not necessarily",
      "from other proteins. These are sequence decisions, not observed peptide detections."
    ),
    artifact_paths[["junction_evidence_pdf"]]
  )
  coverage_figure_html <- embedded_image(
    artifact_paths[["peptide_coverage_plot"]],
    "Peptide mapping tracks",
    paste(
      "Figure 2. Peptide coverage of the supplied references.",
      "Coverage is the union of mapped intervals; overlapping peptides and duplicate",
      "input rows do not increase the covered-residue count.",
      "All valid input sequences contribute, including controls.",
      "Each row uses that reference's coordinates.",
      "Where classes overlap, the strongest junction interpretation is drawn on top;",
      "every occurrence remains in peptide_hits.csv."
    ),
    artifact_paths[["peptide_coverage_pdf"]]
  )
  alignment_figure_html <- embedded_image(
    artifact_paths[["alignment_status_plot"]],
    "Local sequence-similarity overview",
    paste(
      "Figure 3. Local fusion / parent sequence similarity in 25-residue bins.",
      "High-identity bins have at least 50% of their residues aligned and at least",
      "95% identity among paired residues. I/L-equivalent identity is shown separately",
      "when exact identity is below that threshold. Gaps are excluded from identity.",
      "The displayed similarity does not establish parent ancestry."
    ),
    artifact_paths[["alignment_status_pdf"]]
  )
  generated <- html_text(
    result$manifest$value[match("generated_at_utc", result$manifest$key)]
  )
  matching_mode <- if (isTRUE(result$config$il_equivalent)) {
    "I/L-equivalent"
  } else {
    "exact amino-acid"
  }
  title <- "FusionPep: peptide mapping decision report"
  css <- paste0(
    ":root{color-scheme:light;--ink:#202b36;--muted:#52616e;--accent:#155e75;",
    "--line:#d8e0e5;--soft:#f3f6f8}",
    "*{box-sizing:border-box}html{scroll-behavior:smooth}",
    "body{font-family:Arial,Helvetica,sans-serif;line-height:1.55;color:var(--ink);",
    "background:#edf1f3;margin:0;font-size:15px}",
    ".report-layout{display:grid;grid-template-columns:190px minmax(0,1fr);gap:26px;",
    "max-width:1480px;margin:24px auto;padding:0 24px;align-items:start}",
    ".report{min-width:0;background:white;padding:36px 40px;box-shadow:0 2px 18px #182c3b12}",
    ".report-sidebar{position:sticky;top:24px;max-height:calc(100vh - 48px);overflow-y:auto}",
    ".report-header{border-top:5px solid var(--accent);padding-top:22px}",
    ".eyebrow{text-transform:uppercase;letter-spacing:.12em;font-size:11px;",
    "font-weight:700;color:var(--accent);margin:0 0 10px}",
    "h1{font-family:Georgia,serif;font-size:34px;line-height:1.16;",
    "font-weight:normal;margin:0 0 16px;max-width:820px}",
    ".lead{font-size:16px;max-width:900px;margin:0 0 14px}",
    ".run-meta{font-size:12px;color:var(--muted);margin:0 0 20px}",
    "a{color:var(--accent);text-underline-offset:3px}",
    "a:focus-visible,summary:focus-visible{outline:2px solid var(--accent);outline-offset:3px}",
    ".skip-link{position:fixed;left:16px;top:-100px;padding:10px;background:white;z-index:30}",
    ".skip-link:focus{top:10px}",
    ".report-sidebar .contents{border:0;margin:0;padding:12px 0;border-radius:0}",
    ".contents>summary{padding:0 12px 8px;color:var(--muted);font-size:12px}",
    ".contents[open]>summary{margin-bottom:6px}",
    ".contents nav{display:flex;flex-direction:column;gap:3px}",
    ".contents a{display:block;text-decoration:none;font-size:13px;line-height:1.4;",
    "padding:9px 12px;border-left:3px solid transparent;color:var(--muted)}",
    ".contents a:hover{background:#e1e9ed;color:var(--ink)}",
    ".contents a[aria-current]{border-left-color:var(--accent);background:#dfeaf0;",
    "color:var(--accent);font-weight:700}",
    "section{margin-top:32px}section,[id]{scroll-margin-top:24px}",
    "h2{font-size:22px;line-height:1.3;border-bottom:1px solid var(--line);",
    "padding-bottom:10px;margin:0 0 14px;font-weight:600}",
    "h3{font-size:16px;margin:22px 0 10px}p{margin:10px 0}",
    ".cards{display:grid;grid-template-columns:repeat(5,minmax(0,1fr));",
    "gap:0;border:1px solid var(--line);margin:22px 0}",
    ".card{padding:16px 14px;border-right:1px solid var(--line)}",
    ".card:last-child{border-right:0}.card-value{font-size:28px;line-height:1.15;",
    "font-weight:600;font-variant-numeric:tabular-nums;color:var(--accent)}",
    ".card-label{font-size:12px;font-weight:700;margin-top:8px}",
    ".card-detail{font-size:11px;line-height:1.4;color:var(--muted);margin-top:4px}",
    ".warning{background:#fff7e6;border-left:3px solid #b77800;padding:12px 16px}",
    ".validation-note,.empty,.table-truncated{font-size:12px;color:var(--muted)}",
    ".validation-note{margin-bottom:0}.table-truncated{padding:8px 10px;",
    "background:var(--soft);border-left:3px solid var(--accent)}",
    ".table-wrap{overflow-x:auto;margin:14px 0}",
    "table{border-collapse:collapse;width:100%;font-size:12px;line-height:1.45}",
    "th{background:var(--soft);font-weight:600;text-align:left;color:#253d4b;",
    "border-top:1px solid var(--line);border-bottom:1px solid var(--line)}",
    "td,th{padding:9px 10px;vertical-align:top}",
    "td{border-bottom:1px solid #e5eaee;overflow-wrap:anywhere}",
    "td.numeric{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}",
    "td.reference{white-space:nowrap}",
    "td.peptide{font-family:monospace;min-width:130px;max-width:260px;word-break:break-word}",
    "tbody tr:nth-child(even){background:#fafcfd}",
    "figure{margin:24px 0 28px;break-inside:avoid;page-break-inside:avoid}",
    "figure img{display:block;width:100%;height:auto}",
    "figcaption{font-size:12px;line-height:1.55;color:#384957;margin-top:8px}",
    ".figure-downloads{font-size:11px;margin:7px 0;color:var(--muted)}",
    ".figure-downloads span{padding:0 8px}",
    "details{border:1px solid var(--line);padding:12px 16px;margin:14px 0;border-radius:3px}",
    "summary{cursor:pointer;font-size:13px;font-weight:600;color:#253d4b}",
    "details[open] summary{margin-bottom:12px}",
    ".sequence,.alignment{font-family:monospace;font-size:11px;line-height:1.85;",
    "white-space:pre;overflow-x:auto;padding:12px;background:#f8fafb}",
    ".alignment{margin:0 0 6px;break-inside:avoid;page-break-inside:avoid}",
    ".position{color:#687681}.mapped-residue{font-weight:700;color:#14202b}",
    ".parent-mapping{background:#56B4E955}.fusion-mapping{background:#7F8C8D40}",
    ".junction-flank-fail{background:#E69F0055}.junction-shared{background:#CC79A766}",
    ".junction-candidate{background:#009E7360}",
    ".alignment-key{font-size:11px;color:var(--muted)}",
    ".manifest-file{width:30%;min-width:175px}.manifest-file code{font-size:11px}",
    ".output-manifest td:last-child{min-width:230px}.output-manifest .numeric{width:75px}",
    ".report-footer{border-top:1px solid var(--line);padding-top:14px;margin-top:32px;",
    "font-size:11px;color:var(--muted)}",
    "@media(max-width:900px){body{background:white}.report-layout{display:block;margin:0;padding:0}",
    ".report-sidebar{top:0;z-index:20;max-height:70vh;background:white;",
    "border-bottom:1px solid var(--line);box-shadow:0 2px 8px #182c3b0a}",
    ".report-sidebar .contents{padding:0}.contents>summary{padding:12px 20px;font-size:13px}",
    ".contents[open]>summary{margin-bottom:0}.contents nav{padding:0 12px 12px}",
    ".report{padding:28px 24px;box-shadow:none}section,[id]{scroll-margin-top:68px}}",
    "@media(max-width:600px){.report{padding:22px 16px}h1{font-size:27px}",
    ".cards{grid-template-columns:repeat(2,minmax(0,1fr))}",
    ".card{border-bottom:1px solid var(--line)}",
    ".card:nth-child(even){border-right:0}.card:last-child{grid-column:1 / -1;border-right:0}",
    "td,th{padding:7px 8px}table{font-size:11px;min-width:720px}}",
    "@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}}",
    "@page{size:A4;margin:13mm 15mm}",
    "@media print{html{scroll-behavior:auto}body{background:white;font-size:9pt}",
    ".report-layout{display:block;margin:0;padding:0}.report-sidebar,.skip-link{display:none}",
    ".report{width:auto;max-width:none;margin:0;padding:0;box-shadow:none}",
    "h1{font-size:24pt}h2{font-size:14pt;break-after:avoid}h3{break-after:avoid}",
    ".lead{font-size:10pt}.contents,.figure-downloads{display:none}",
    ".cards{grid-template-columns:repeat(5,minmax(0,1fr));break-inside:avoid}",
    ".card{padding:8px}.card-value{font-size:18pt}.card-label{font-size:7pt}",
    ".card-detail{font-size:6.5pt}section{margin-top:20px}",
    "table{font-size:7pt;min-width:0}td,th{padding:5px 6px}.table-wrap{overflow:visible}",
    "thead{display:table-header-group}tr{break-inside:avoid}",
    "figcaption{font-size:8pt}figure{margin:14px 0}",
    "details{border:0;padding:0}details>summary{margin:12px 0;font-size:9pt}",
    "details:not([open]){display:none}#sequences:not(:has(details[open])){display:none}",
    "details>summary::marker{content:''}",
    ".sequence,.alignment{font-size:6.5pt;overflow:visible;padding:4px;line-height:1.7}",
    ".manifest-file code{font-size:6.5pt}a{color:inherit;text-decoration:none}",
    "*{-webkit-print-color-adjust:exact;print-color-adjust:exact}}"
  )
  navigation_script <- paste(
    "(() => {",
    "  const contents = document.querySelector('.contents');",
    "  const compact = window.matchMedia('(max-width: 900px)');",
    "  const resize = () => { contents.open = !compact.matches; };",
    "  compact.addEventListener('change', resize);",
    "  resize();",
    "  const links = Array.from(contents.querySelectorAll('a'));",
    "  const sections = links.map(link => document.getElementById(link.hash.slice(1)));",
    "  const update = () => {",
    "    let current = 0;",
    "    sections.forEach((section, index) => {",
    "      if (section.getBoundingClientRect().top <= 120) current = index;",
    "    });",
    "    if (innerHeight + scrollY >= document.documentElement.scrollHeight - 2) {",
    "      current = links.length - 1;",
    "    }",
    "    links.forEach((link, index) => {",
    "      if (index === current) link.setAttribute('aria-current', 'location');",
    "      else link.removeAttribute('aria-current');",
    "    });",
    "  };",
    "  const reveal = hash => {",
    "    const target = document.getElementById(hash.slice(1));",
    "    if (!target) return;",
    "    let detail = target.closest('details');",
    "    while (detail) {",
    "      detail.open = true;",
    "      detail = detail.parentElement.closest('details');",
    "    }",
    "  };",
    "  document.addEventListener('click', event => {",
    "    const link = event.target.closest('a[href^=\"#\"]');",
    "    if (!link) return;",
    "    reveal(link.hash);",
    "    if (compact.matches) contents.open = false;",
    "  });",
    "  window.addEventListener('hashchange', () => reveal(location.hash));",
    "  window.addEventListener('scroll', update, {passive: true});",
    "  window.addEventListener('resize', update);",
    "  reveal(location.hash);",
    "  update();",
    "})();",
    sep = "\n"
  )
  body <- paste0(
    '<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    "<title>", html_text(title), "</title><style>", css, "</style></head>",
    '<body><a class="skip-link" href="#overview">Skip to report</a>',
    '<div class="report-layout"><aside class="report-sidebar">',
    '<details class="contents" open><summary>Contents</summary>',
    '<nav aria-label="Report sections">',
    '<a href="#overview" aria-current="location">Overview</a>',
    '<a href="#junction">Junction decision</a>',
    '<a href="#peptides">Peptides</a><a href="#coverage">Coverage</a>',
    '<a href="#sequences">Reference sequences</a>',
    '<a href="#alignment">Parent comparison</a>',
    '<a href="#audit">Inputs and manifests</a>',
    '</nav></details></aside><main class="report">',
    '<header id="overview" class="report-header" tabindex="-1">',
    '<p class="eyebrow">FusionPep / sequence audit</p>',
    "<h1>", html_text(title), "</h1>",
    '<p class="lead">Peptide matches in the supplied fusion and parent references, ',
    "with junction and flank assessment. These sequence matches do not establish ",
    "peptide detection or biological fusion expression.</p>",
    '<p class="run-meta">Analysis generated ', generated,
    " | Peptide matching: ", html_text(matching_mode),
    " | Coordinates: 1-based, inclusive</p></header>",
    summary_cards_html(result),
    warning_html,
    '<section id="junction"><h2>Junction decision</h2>',
    junction_assessment_html,
    "<p>A <b>junction candidate</b> crosses a supplied boundary, meets the minimum ",
    "flank rule on both sides, and is absent from both supplied parent records. ",
    "Inserted residues do not count toward either flank. This remains a ",
    "sequence-level screening decision.</p>",
    junction_figure_html,
    "<details><summary>Supplied junction metadata</summary>",
    "<p>The table contains the supplied sequence-junction metadata. Parent labels ",
    "and coordinates come from that input; the mapper checks internal sequence ",
    "consistency rather than independently establishing biological ancestry.</p>",
    table_to_html(junction_table), "</details>",
    junction_evaluation_html, "</section>",
    '<section id="peptides"><h2>All peptide interpretations</h2>',
    "<p>Presence describes occurrence in the three supplied records. Input roles ",
    "are retained annotations; they do not establish detection and do not filter ",
    "the analysis. Counts refer to distinct matching sequences after the selected ",
    "I/L normalization, unless explicitly labeled as input rows or occurrences.</p>",
    table_to_html(interpretation_table), "</section>",
    '<section id="coverage"><h2>Reference coverage</h2>',
    coverage_figure_html, table_to_html(result$coverage), "</section>",
    '<section id="sequences"><h2>Annotated reference sequences</h2>',
    "<p>Open a reference to inspect residue-level mappings. ",
    html_text(sequence_annotation_note),
    " Hover over a highlighted residue for the input peptide. ",
    "Spaces group ten residues; numbers give the first and last positions on each line.</p>",
    paste(sequence_sections, collapse = ""), "</section>",
    '<section id="alignment"><h2>Parent/fusion comparison</h2>',
    alignment_figure_html,
    "<p>The intervals below identify high-identity fusion bins and the parent ",
    "residues paired within those bins. Local identity is calculated among paired ",
    "residues. Unaligned positions and gaps remain explicit in the alignment CSVs.</p>",
    table_to_html(alignment_table),
    "<details><summary>Candidate non-match regions</summary>",
    "<p>Parent-only/fusion-only alignment positions and mismatches are candidates ",
    "for review. A parent-only peptide match does not establish biological domain ",
    "loss. This table includes global and local alignments; ",
    "alignment_regions.csv retains every region.</p>",
    table_to_html(candidate_regions, max_rows = 100L), "</details>",
    paste(alignment_sections, collapse = ""), "</section>",
    '<section id="audit"><h2>Inputs and manifests</h2>',
    "<h3>Reference records</h3>", table_to_html(references),
    "<details><summary>Input normalization audit</summary>",
    if (length(result$warnings) == 0L) {
      '<p class="validation-note">No input normalization warnings.</p>'
    } else {
      warning_html
    },
    table_to_html(result$normalized_peptides), "</details>",
    "<details><summary>Run manifest: inputs, versions and settings</summary>",
    "<p>This records the analysis environment and input hashes for this run. ",
    "The file descriptions below identify the saved tables that contain each result.</p>",
    table_to_html(result$manifest), "</details>",
    '<details id="files"><summary>Output manifest: file contents and row counts</summary>',
    output_manifest_html(result, artifact_paths), "</details></section>",
    '<footer class="report-footer">Interpretation applies to the supplied records ',
    "and analysis settings. Absence from these two parents is not proteome-wide ",
    "uniqueness. PSM validation, FDR control and biological confirmation require ",
    "experimental evidence outside this sequence audit.</footer>",
    "</main></div><script>", navigation_script, "</script></body></html>"
  )
  temporary_report <- tempfile(pattern = ".fusion-report-", tmpdir = output_dir)
  on.exit(unlink(temporary_report), add = TRUE)
  writeLines(body, temporary_report, useBytes = TRUE)
  if (!file.rename(temporary_report, report_path)) {
    stopf("Could not finalize HTML report: %s", report_path)
  }
  report_path
}
