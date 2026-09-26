# Diagnostic figures for FusionPep.
#
# Plots draw the analysis result and reuse the core junction and alignment
# functions; they do not make separate scientific decisions.

peptide_display_colors <- c(
  parent_reference_mapping = "#56B4E9",
  fusion_mapping_not_junction_spanning = "#7F8C8D",
  fusion_mapping_junction_unassessed = "#7F8C8D",
  junction_crossing_below_flank_threshold = "#E69F00",
  junction_spanning_shared_parent = "#CC79A7",
  junction_spanning_candidate = "#009E73"
)

peptide_display_labels <- c(
  parent_reference_mapping = "Parent reference",
  fusion_mapping_not_junction_spanning = "Fusion; no junction crossing",
  fusion_mapping_junction_unassessed = "Fusion; junction not assessed",
  junction_crossing_below_flank_threshold = "Crosses; flank rule fails",
  junction_spanning_shared_parent = "Crosses; also in a parent",
  junction_spanning_candidate = "Junction candidate"
)

publication_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 9, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.subtitle = ggplot2::element_text(
        size = 8.5, colour = "#374151", lineheight = 1.15,
        margin = ggplot2::margin(b = 9)
      ),
      plot.caption = ggplot2::element_text(
        size = 7.5, colour = "#374151", hjust = 0, lineheight = 1.15,
        margin = ggplot2::margin(t = 9)
      ),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(10, 12, 10, 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#E5E7EB", linewidth = 0.25
      ),
      axis.text = ggplot2::element_text(size = 8, colour = "#20252B"),
      axis.title = ggplot2::element_text(size = 8.5),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 7)),
      legend.position = "bottom",
      legend.justification = "left",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 8),
      legend.key.size = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.margin = ggplot2::margin(t = 5),
      strip.text = ggplot2::element_text(
        size = 8, face = "bold", hjust = 0, lineheight = 1.2
      ),
      panel.spacing = grid::unit(7, "mm")
    )
}

peptide_matching_caption <- function(hits) {
  if (nrow(hits) == 0L) return("No mapped peptide occurrences.")
  if (!"match_basis" %in% names(hits)) {
    return("Peptide matching mode is not recorded in the supplied hits.")
  }
  paste0(
    "Peptide matching: ",
    paste(unique(hits$match_basis), collapse = ", "),
    "."
  )
}

aggregate_coverage_segments <- function(peptide_plot) {
  if (nrow(peptide_plot) == 0L) {
    return(peptide_plot)
  }
  group_key <- paste(
    as.character(peptide_plot$sequence_id),
    peptide_plot$coverage_class,
    sep = "\r"
  )
  grouped_rows <- split(seq_len(nrow(peptide_plot)), group_key)
  segments <- lapply(grouped_rows, function(rows) {
    intervals <- merge_intervals(
      peptide_plot$start[rows],
      peptide_plot$end[rows]
    )
    data.frame(
      sequence_id = rep(
        as.character(peptide_plot$sequence_id[[rows[[1L]]]]),
        nrow(intervals)
      ),
      start = intervals$start,
      end = intervals$end,
      coverage_class = rep(
        peptide_plot$coverage_class[[rows[[1L]]]],
        nrow(intervals)
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, segments)
}

plot_peptide_coverage <- function(peptide_hits,
                                  sequence_text,
                                  junctions = NULL) {
  sequence_ids <- names(sequence_text)
  if (is.null(sequence_ids)) {
    sequence_ids <- paste0("Sequence", seq_along(sequence_text))
    names(sequence_text) <- sequence_ids
  }
  coverage <- calculate_coverage(sequence_text, sequence_ids, peptide_hits)$summary
  track <- data.frame(
    sequence_id = sequence_ids,
    y = rev(seq_along(sequence_ids)),
    start = 1,
    end = nchar(sequence_text),
    stringsAsFactors = FALSE
  )
  track_labels <- sprintf(
    "%s (%d aa)\n%d covered | %.2f%%",
    sequence_ids, track$end, coverage$covered_residues, coverage$coverage_percent
  )
  junction_assessed <- !is.null(junctions) && nrow(junctions) > 0L

  plot <- ggplot2::ggplot(track) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = start - 0.5, xmax = end + 0.5,
        ymin = y - 0.13, ymax = y + 0.13
      ),
      fill = "#E5E7EB"
    ) +
    ggplot2::labs(
      title = "Peptide coverage across the supplied references",
      subtitle = paste(
        peptide_matching_caption(peptide_hits),
        "All valid input sequences contribute to coverage, including controls.",
        sep = "\n"
      ),
      caption = if (junction_assessed) {
        paste(
          "Coverage counts each residue once. A junction candidate is absent from",
          "both supplied parents; sequence matching does not establish detection.",
          sep = "\n"
        )
      } else {
        "Coverage counts each residue once. Junction status is not assessed."
      },
      x = "Residue position within each reference (1-based)",
      y = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = track$y, labels = track_labels,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, max(track$end) + 1),
      expand = ggplot2::expansion(mult = c(0, 0.015))
    ) +
    publication_plot_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(lineheight = 1.25))

  if (nrow(peptide_hits) > 0L) {
    require_data_frame_columns(peptide_hits, "display_class", "peptide_hits")
    peptide_plot <- peptide_hits[, c("sequence_id", "start", "end"), drop = FALSE]
    peptide_plot$coverage_class <- peptide_hits$display_class
    peptide_plot <- aggregate_coverage_segments(peptide_plot)
    # Draw stronger junction evidence last where classes overlap.
    peptide_plot <- peptide_plot[
      order(match(peptide_plot$coverage_class, names(peptide_display_colors))),
      ,
      drop = FALSE
    ]
    peptide_plot$y <- track$y[match(peptide_plot$sequence_id, track$sequence_id)]
    present_classes <- names(peptide_display_colors)[
      names(peptide_display_colors) %in% peptide_plot$coverage_class
    ]
    plot <- plot +
      ggplot2::geom_rect(
        data = peptide_plot,
        ggplot2::aes(
          xmin = start - 0.5, xmax = end + 0.5,
          ymin = y - 0.17, ymax = y + 0.17, fill = coverage_class
        )
      ) +
      ggplot2::scale_fill_manual(
        values = peptide_display_colors,
        breaks = present_classes,
        labels = peptide_display_labels[present_classes]
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
  }

  if (junction_assessed && "Fusion" %in% sequence_ids) {
    marker <- data.frame(
      x = (junctions$fusion_left_position + junctions$fusion_right_position) / 2,
      y = track$y[match("Fusion", track$sequence_id)],
      stringsAsFactors = FALSE
    )
    plot <- plot +
      ggplot2::geom_point(
        data = marker, ggplot2::aes(x = x, y = y + 0.29),
        inherit.aes = FALSE, shape = 25, size = 1.8,
        fill = "#20252B", colour = "#20252B", stroke = 0.2
      )
    plot$labels$caption <- paste(
      plot$labels$caption,
      "Black triangles mark supplied junctions.",
      sep = "\n"
    )
  }
  plot
}

plot_junction_evidence <- function(peptide_hits,
                                   sequence_text,
                                   junctions,
                                   flank_window = 20L,
                                   max_candidates = 100L) {
  for (name in c("flank_window", "max_candidates")) {
    value <- get(name)
    minimum <- if (identical(name, "flank_window")) 0L else 1L
    if (!is.numeric(value) || length(value) != 1L ||
        is.na(value) || !is.finite(value) ||
        value < minimum || value != floor(value) ||
        value > .Machine$integer.max) {
      stopf("%s must be one integer of at least %d.", name, minimum)
    }
  }
  if (is.null(junctions) || nrow(junctions) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text", x = 0, y = 0,
          label = "No junction metadata supplied.\nJunction status is not assessed.",
          size = 3.2, colour = "#374151", lineheight = 1.4
        ) +
        ggplot2::labs(title = "Fusion-junction peptide matches") +
        ggplot2::theme_void(base_size = 9) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 11, face = "bold"),
          plot.margin = ggplot2::margin(10, 12, 10, 10),
          plot.background = ggplot2::element_rect(fill = "white", colour = NA)
        )
    )
  }
  require_data_frame_columns(
    peptide_hits,
    c("sequence_id", "matching_peptide", "matched_subsequence", "start", "end"),
    "peptide_hits"
  )
  evaluations <- evaluate_peptide_junctions(peptide_hits, junctions)
  candidates_by_junction <- sequence_rows <- windows <- boundaries <-
    insertions <- parent_labels <- vector("list", nrow(junctions))
  panel_labels <- character(nrow(junctions))
  priority <- c(
    "junction_spanning_candidate",
    "junction_spanning_shared_parent",
    "junction_crossing_below_flank_threshold"
  )

  for (index in seq_len(nrow(junctions))) {
    junction <- junctions[index, , drop = FALSE]
    current <- evaluations[
      evaluations$junction_id == junction$junction_id, , drop = FALSE
    ]
    # Reuse the core decision for this junction, not the hit's best junction.
    annotated <- apply_junction_evaluations(peptide_hits, current)
    classes <- classify_peptide_hit_display_classes(annotated)
    selected <- current$hit_row_id
    if (length(selected) > 0L) {
      ordering <- order(
        match(classes[selected], priority),
        -pmin(current$junction_left_flank, current$junction_right_flank),
        peptide_hits$start[selected],
        peptide_hits$matching_peptide[selected]
      )
      selected <- selected[utils::head(ordering, max_candidates)]
    }
    panel <- paste0(
      junction$junction_id, " | minimum ", junction$min_flank_aa,
      " aa on each side\n", length(selected), " of ", nrow(current),
      " crossing occurrences shown"
    )
    panel_labels[[index]] <- panel
    fusion_sequence <- sequence_text[[junction$fusion_id]]
    window_start <- max(1L, junction$fusion_left_position - flank_window)
    window_end <- min(nchar(fusion_sequence), junction$fusion_right_position + flank_window)
    positions <- seq.int(window_start, window_end)
    sequence_rows[[index]] <- data.frame(
      panel = panel, position = positions, row_label = "Fusion sequence",
      amino_acid = substring(fusion_sequence, positions, positions),
      stringsAsFactors = FALSE
    )
    windows[[index]] <- data.frame(
      panel = panel, x = c(window_start - 0.5, window_end + 0.5),
      row_label = "Fusion sequence", stringsAsFactors = FALSE
    )
    boundary_positions <- unique(c(
      junction$fusion_left_position + 0.5,
      junction$fusion_right_position - 0.5
    ))
    boundaries[[index]] <- data.frame(
      panel = panel, x = boundary_positions, stringsAsFactors = FALSE
    )
    if (junction$fusion_right_position > junction$fusion_left_position + 1L) {
      insertions[[index]] <- data.frame(
        panel = panel, xmin = junction$fusion_left_position + 0.5,
        xmax = junction$fusion_right_position - 0.5,
        stringsAsFactors = FALSE
      )
    }
    parent_labels[[index]] <- data.frame(
      panel = panel, x = c(window_start, window_end),
      label = c(junction$upstream_parent, junction$downstream_parent),
      hjust = c(0, 1), stringsAsFactors = FALSE
    )

    if (length(selected) == 0L) next
    candidate <- peptide_hits[selected, , drop = FALSE]
    candidate$panel <- panel
    candidate$plot_class <- classes[selected]
    candidate$left_flank <- annotated$junction_left_flank[selected]
    candidate$right_flank <- annotated$junction_right_flank[selected]
    candidate$decision <- ifelse(
      annotated$junction_flank_pass[selected],
      ifelse(
        candidate$plot_class == "junction_spanning_shared_parent",
        "pass; parent match", "pass"
      ),
      "fail"
    )
    peptide_labels <- vapply(candidate$matched_subsequence, function(value) {
      starts <- seq.int(1L, nchar(value), by = 28L)
      paste(substring(value, starts, starts + 27L), collapse = "\n")
    }, character(1))
    candidate$row_label <- paste0(
      peptide_labels, "\n", candidate$start, "-", candidate$end,
      "; flanks ", candidate$left_flank, " / ", candidate$right_flank,
      " aa: ", candidate$decision
    )
    candidate$xmin <- pmax(candidate$start - 0.5, window_start - 0.5)
    candidate$xmax <- pmin(candidate$end + 0.5, window_end + 0.5)
    candidate$clipped_left <- candidate$start < window_start
    candidate$clipped_right <- candidate$end > window_end
    candidates_by_junction[[index]] <- candidate
  }

  sequence_rows <- do.call(rbind, sequence_rows)
  windows <- do.call(rbind, windows)
  boundaries <- do.call(rbind, boundaries)
  insertions <- do.call(rbind, insertions)
  parent_labels <- do.call(rbind, parent_labels)
  candidates <- do.call(rbind, candidates_by_junction)
  row_levels <- c(
    "Fusion sequence",
    if (!is.null(candidates)) rev(unique(candidates$row_label)) else character()
  )
  sequence_rows$row_label <- factor(sequence_rows$row_label, levels = row_levels)
  sequence_rows$panel <- factor(sequence_rows$panel, levels = panel_labels)

  plot <- ggplot2::ggplot(sequence_rows) +
    ggplot2::geom_blank(
      data = windows, ggplot2::aes(x = x, y = row_label), inherit.aes = FALSE
    )
  if (!is.null(insertions)) {
    plot <- plot +
      ggplot2::geom_rect(
        data = insertions,
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = "#F0E442", alpha = 0.22
      )
  }
  plot <- plot +
    ggplot2::geom_vline(
      data = boundaries, ggplot2::aes(xintercept = x),
      colour = "#59616B", linetype = "dashed", linewidth = 0.35
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = position, y = row_label, label = amino_acid),
      family = "mono", size = 2.3, colour = "#20252B"
    ) +
    ggplot2::geom_text(
      data = parent_labels,
      ggplot2::aes(x = x, y = Inf, label = label, hjust = hjust),
      inherit.aes = FALSE, vjust = 1.5, fontface = "bold",
      size = 2.6, colour = "#374151"
    ) +
    ggplot2::facet_wrap(~panel, ncol = 1, scales = "free") +
    ggplot2::scale_y_discrete(
      limits = function(values) row_levels[row_levels %in% values],
      expand = ggplot2::expansion(add = c(0.55, 0.8))
    ) +
    ggplot2::scale_x_continuous(
      breaks = function(limits) {
        step <- max(1, 5 * ceiling(diff(limits) / 60))
        seq(ceiling(limits[[1L]] / step) * step, limits[[2L]], by = step)
      },
      expand = ggplot2::expansion(mult = c(0.015, 0.015))
    ) +
    ggplot2::labs(
      title = "Peptide matches across the supplied junction",
      subtitle = paste(
        peptide_matching_caption(peptide_hits),
        "Flanks and decisions are evaluated separately for each junction.",
        sep = "\n"
      ),
      caption = paste(
        "Dashed lines mark the supplied parent boundaries; yellow marks inserted residues.",
        "Labels give reference subsequences, full coordinates and flanks. Arrows mark clipped ends.",
        "Candidates are absent from both supplied parents. Full evaluations: junction_evaluations.csv.",
        sep = "\n"
      ),
      x = "Fusion residue position (1-based)", y = NULL
    ) +
    publication_plot_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.2, lineheight = 1.2),
      panel.grid.major.x = ggplot2::element_blank()
    )

  if (!is.null(candidates) && nrow(candidates) > 0L) {
    candidates$row_label <- factor(candidates$row_label, levels = row_levels)
    present_classes <- names(peptide_display_colors)[
      names(peptide_display_colors) %in% candidates$plot_class
    ]
    plot <- plot +
      ggplot2::geom_tile(
        data = candidates,
        ggplot2::aes(
          x = (xmin + xmax) / 2, y = row_label,
          width = xmax - xmin, fill = plot_class
        ),
        height = 0.28, inherit.aes = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = peptide_display_colors,
        breaks = present_classes,
        labels = peptide_display_labels[present_classes]
      ) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE))
    for (side in c("left", "right")) {
      clipped <- candidates[
        candidates[[paste0("clipped_", side)]], , drop = FALSE
      ]
      if (nrow(clipped) == 0L) next
      clipped$tip <- if (identical(side, "left")) clipped$xmin else clipped$xmax
      clipped$tail <- clipped$tip + if (identical(side, "left")) 1 else -1
      plot <- plot +
        ggplot2::geom_segment(
          data = clipped,
          ggplot2::aes(x = tail, xend = tip, y = row_label, yend = row_label),
          inherit.aes = FALSE, linewidth = 0.35, colour = "#20252B",
          arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")
        )
    }
  }
  plot
}

plot_alignment_status <- function(alignment_columns, junctions = NULL) {
  bins <- alignment_overview_bins(alignment_columns)
  if (nrow(bins) == 0L) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }
  bins$parent <- factor(bins$parent, levels = rev(unique(bins$parent)))
  bins$status <- factor(
    bins$status,
    levels = c(
      "high_exact_identity",
      "high_il_identity",
      "aligned_with_differences",
      "partial_local_alignment",
      "no_local_alignment"
    )
  )
  plot <- ggplot2::ggplot(
    bins,
    ggplot2::aes(
      x = bin_mid,
      y = parent,
      width = bin_end - bin_start + 1L,
      height = 0.82,
      fill = status
    )
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.2, show.legend = TRUE) +
    ggplot2::scale_fill_manual(
      values = c(
        high_exact_identity = "#009E73",
        high_il_identity = "#0072B2",
        aligned_with_differences = "#E69F00",
        partial_local_alignment = "#F4A582",
        no_local_alignment = "#D9D9D9"
      ),
      labels = c(
        high_exact_identity = "Exact identity >=95%",
        high_il_identity = "I/L identity >=95%",
        aligned_with_differences = "Lower identity",
        partial_local_alignment = "<50% of bin aligned",
        no_local_alignment = "No local alignment"
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Fusion / parent sequence similarity",
      subtitle = paste(
        "Local alignments summarized in 25-residue fusion bins.",
        "Sequence similarity does not establish parent ancestry.",
        sep = "\n"
      ),
      caption = paste(
        "High identity requires at least 50% of the bin to align. Identity excludes gaps.",
        if (!is.null(junctions) && nrow(junctions) > 0L) {
          "Dashed lines mark supplied junction boundaries."
        } else {
          "No junction metadata supplied."
        },
        sep = "\n"
      ),
      x = "Fusion residue position (1-based)",
      y = NULL,
      fill = "Local alignment summary"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.005, 0.02))
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 2, byrow = TRUE)) +
    publication_plot_theme()
  if (!is.null(junctions) && nrow(junctions) > 0L) {
    boundary_data <- data.frame(
      x = c(
        junctions$fusion_left_position + 0.5,
        junctions$fusion_right_position - 0.5
      )
    )
    plot <- plot +
      ggplot2::geom_vline(
        data = boundary_data,
        ggplot2::aes(xintercept = x),
        inherit.aes = FALSE,
        colour = "#20252B",
        linetype = "dashed",
        linewidth = 0.35
      )
  }
  plot
}

# Save the same drawing at manuscript width in raster and vector formats.
write_publication_figure <- function(plot, png_path, pdf_path, height) {
  width <- 180 / 25.4
  ggplot2::ggsave(
    png_path, plot, width = width, height = height,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    pdf_path, plot, device = grDevices::cairo_pdf,
    width = width, height = height, units = "in",
    bg = "white", limitsize = FALSE
  )
}
