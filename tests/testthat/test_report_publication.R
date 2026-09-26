publication_result <- run_fusion_analysis(
  file.path(project_root, "input", "sequences.fasta"),
  file.path(project_root, "input", "peptides.csv"),
  junction_file = file.path(project_root, "input", "fusion_junctions.csv")
)

testthat::test_that("the supplied example retains its peptide decisions and coverage", {
  hits <- publication_result$peptide_hits
  fusion <- hits[hits$sequence_id == "Fusion", , drop = FALSE]
  testthat::expect_identical(nrow(hits), 9L)
  testthat::expect_equal(fusion$start, c(390L, 386L, 390L, 398L))
  testthat::expect_equal(fusion$end, c(407L, 397L, 394L, 407L))
  testthat::expect_identical(
    fusion$display_class,
    c(
      "junction_spanning_candidate",
      "junction_crossing_below_flank_threshold",
      "fusion_mapping_not_junction_spanning",
      "fusion_mapping_not_junction_spanning"
    )
  )
  testthat::expect_equal(
    publication_result$coverage$covered_residues, c(22L, 52L, 68L)
  )
  testthat::expect_equal(
    publication_result$coverage$coverage_percent, c(2.76, 5.9, 14.72)
  )
  testthat::expect_equal(
    publication_result$junction_evaluations$junction_left_flank, c(5L, 9L)
  )
  testthat::expect_equal(
    publication_result$junction_evaluations$junction_right_flank, c(12L, 2L)
  )
})

testthat::test_that("junction legend labels stay attached to every category subset", {
  sequences <- c(Fusion = "AAAXBBB", ParentA = "YAAXBBY", ParentB = "BBBB")
  junctions <- data.frame(
    fusion_id = "Fusion", junction_id = "synthetic_join",
    upstream_parent = "ParentA", downstream_parent = "ParentB",
    fusion_left_position = 3L, fusion_right_position = 5L, min_flank_aa = 2L,
    stringsAsFactors = FALSE
  )
  peptides <- c("AAAXBBB", "AAXBB", "AXBB")
  mapping <- map_peptides_to_references(
    sequences, normalize_peptides(peptides, il_equivalent = FALSE),
    junctions, il_equivalent = FALSE
  )
  expected_labels <- c(
    junction_spanning_candidate = "Junction candidate",
    junction_spanning_shared_parent = "Crosses; also in a parent",
    junction_crossing_below_flank_threshold = "Crosses; flank rule fails"
  )
  for (mask in seq_len(7L)) {
    selected <- peptides[as.logical(intToBits(mask)[seq_len(3L)])]
    hits <- mapping$peptide_hits[
      mapping$peptide_hits$matching_peptide %in% selected, , drop = FALSE
    ]
    plot <- plot_junction_evidence(hits, sequences, junctions)
    scale <- ggplot2::ggplot_build(plot)$plot$scales$get_scales("fill")
    breaks <- scale$get_breaks()
    testthat::expect_equal(
      as.character(scale$get_labels()),
      unname(expected_labels[breaks]),
      info = paste("Category subset", mask)
    )
    testthat::expect_length(breaks, length(selected))
  }
})

testthat::test_that("each plotted junction retains its own flank decision", {
  sequences <- c(Fusion = "AAAAAAAAAAAA", ParentA = "AAA", ParentB = "AAAA")
  junctions <- data.frame(
    fusion_id = c("Fusion", "Fusion"),
    junction_id = c("first", "second"),
    upstream_parent = c("ParentA", "ParentA"),
    downstream_parent = c("ParentB", "ParentB"),
    fusion_left_position = c(3L, 10L),
    fusion_right_position = c(4L, 11L),
    min_flank_aa = c(3L, 3L),
    stringsAsFactors = FALSE
  )
  mapping <- map_peptides_to_references(
    sequences,
    normalize_peptides("AAAAAAAAAAAA", il_equivalent = FALSE),
    junctions, il_equivalent = FALSE
  )
  plot <- plot_junction_evidence(
    mapping$peptide_hits, sequences, junctions, flank_window = 1L
  )
  tile <- Filter(
    function(layer) inherits(layer$geom, "GeomTile"), plot$layers
  )[[1L]]$data
  testthat::expect_identical(
    tile$plot_class,
    c("junction_spanning_candidate", "junction_crossing_below_flank_threshold")
  )
  testthat::expect_equal(tile$left_flank, c(3L, 10L))
  testthat::expect_equal(tile$right_flank, c(9L, 2L))
  testthat::expect_true(all(tile$clipped_left))
  testthat::expect_identical(tile$clipped_right, c(TRUE, FALSE))
  testthat::expect_silent(ggplot2::ggplot_build(plot))
})

testthat::test_that("single residue peptide coverage has a visible inclusive width", {
  sequences <- c(Fusion = "AXB", ParentA = "AA", ParentB = "BB")
  mapping <- map_peptides_to_references(
    sequences, normalize_peptides("X", il_equivalent = FALSE),
    il_equivalent = FALSE
  )
  plot <- plot_peptide_coverage(mapping$peptide_hits, sequences)
  built <- ggplot2::ggplot_build(plot)
  rectangles <- built$data[[2L]]
  testthat::expect_equal(rectangles$xmin, 1.5)
  testthat::expect_equal(rectangles$xmax, 2.5)
  testthat::expect_equal(rectangles$xmax - rectangles$xmin, 1)
})

testthat::test_that("alignment text puts both sequences in the same columns", {
  columns <- publication_result$alignment_columns
  selected <- columns[
    columns$parent == "ParentB" & columns$alignment_type == "local", , drop = FALSE
  ]
  html <- alignment_block_html(columns, "ParentB", chunk_width = 40L)
  text <- strsplit(html, '<pre class="alignment">', fixed = TRUE)[[1L]][[2L]]
  text <- strsplit(text, "</pre>", fixed = TRUE)[[1L]][[1L]]
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  fusion <- paste(selected$fusion_aa[seq_len(40L)], collapse = "")
  parent <- paste(selected$parent_aa[seq_len(40L)], collapse = "")
  fusion_column <- regexpr(fusion, lines[[1L]], fixed = TRUE)[[1L]]
  parent_column <- regexpr(parent, lines[[3L]], fixed = TRUE)[[1L]]
  testthat::expect_gt(fusion_column, 1L)
  testthat::expect_equal(fusion_column, parent_column)
  testthat::expect_match(lines[[1L]], "^Fusion[[:space:]]+[0-9]+")
  testthat::expect_match(lines[[3L]], "^ParentB[[:space:]]+[0-9]+")
  testthat::expect_match(lines[[1L]], "[0-9]+$")
  testthat::expect_length(
    gregexpr('<pre class="alignment">', html, fixed = TRUE)[[1L]],
    ceiling(nrow(selected) / 40L)
  )
})

testthat::test_that("missing junction metadata remains unassessed in report tables", {
  result <- publication_result
  mapping <- map_peptides_to_references(
    result$sequence_input$sequence_text, result$normalized_peptides
  )
  result[names(mapping)] <- mapping
  result$junctions <- empty_junction_table()
  interpretations <- peptide_interpretation_table(result)
  testthat::expect_true(all(interpretations$crosses_junction == "not assessed"))
  testthat::expect_false(any(interpretations$crosses_junction == "no"))
})

testthat::test_that("report region truncation states the complete row count", {
  report <- write_fusion_report(
    publication_result, tempfile("publication-report-"), list()
  )
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(html, "Showing the first 100 of 203 rows")
  testthat::expect_match(html, "alignment_regions.csv", fixed = TRUE)
  testthat::expect_match(html, "@media print", fixed = TRUE)
  testthat::expect_match(html, 'aria-label="Report sections"', fixed = TRUE)
})

testthat::test_that("empty junction plots retain explicit assessment and count text", {
  result <- publication_result
  no_metadata <- plot_junction_evidence(
    empty_hit_table(), result$sequence_input$sequence_text, empty_junction_table()
  )
  testthat::expect_silent(ggplot2::ggplot_build(no_metadata))
  no_hits <- plot_junction_evidence(
    empty_hit_table(), result$sequence_input$sequence_text, result$junctions
  )
  testthat::expect_true(any(grepl(
    "0 of 0 crossing occurrences shown", no_hits$data$panel, fixed = TRUE
  )))
  testthat::expect_silent(ggplot2::ggplot_build(no_hits))
})

testthat::test_that("direct joins and missing parent coordinates stay distinct", {
  junction <- publication_result$junctions
  junction$inserted_sequence <- ""
  junction$fusion_right_position <- junction$fusion_left_position + 1L
  junction$upstream_parent_position <- NA_integer_
  junction$downstream_parent_position <- NA_integer_
  display <- junction_display_table(junction)
  testthat::expect_identical(display$inserted_residues, "None")
  testthat::expect_identical(display$fusion_boundary, "residue 394 | residue 395")
  testthat::expect_match(display$parent_join, "position not supplied", fixed = TRUE)
})

testthat::test_that("adjacent parent intervals are merged in the report overview", {
  overview <- alignment_overview_table(publication_result)
  parent_a <- overview[overview$parent == "ParentA", , drop = FALSE]
  testthat::expect_identical(parent_a$high_similarity_parent_residues, "1-375")
  testthat::expect_type(parent_a$overall_local_exact_identity, "double")
})

testthat::test_that("limited junction panels state the complete occurrence count", {
  plot <- plot_junction_evidence(
    publication_result$peptide_hits,
    publication_result$sequence_input$sequence_text,
    publication_result$junctions,
    max_candidates = 1L
  )
  testthat::expect_true(all(grepl(
    "1 of 2 crossing occurrences shown", plot$data$panel, fixed = TRUE
  )))
  candidates <- Filter(
    function(layer) inherits(layer$geom, "GeomTile"), plot$layers
  )[[1L]]$data
  testthat::expect_equal(nrow(candidates), 1L)
  testthat::expect_identical(candidates$plot_class, "junction_spanning_candidate")
  testthat::expect_match(plot$labels$subtitle, "I/L-equivalent", fixed = TRUE)
})

testthat::test_that("the output manifest explains files and reports their complete row counts", {
  output <- tempfile("manifest-files-")
  dir.create(output)
  paths <- list(
    peptide_hits = file.path(output, "peptide_hits.csv"),
    junction_evaluations = file.path(output, "junction_evaluations.csv"),
    result_rds = file.path(output, "fusion_peptide_mapper_result.rds")
  )
  write_csv_output(publication_result$peptide_hits, paths$peptide_hits)
  write_csv_output(publication_result$junction_evaluations, paths$junction_evaluations)
  saveRDS(publication_result, paths$result_rds)
  html <- output_manifest_html(publication_result, paths)
  testthat::expect_match(
    html, 'peptide_hits.csv</code></a></td><td class="numeric">9</td>', fixed = TRUE
  )
  testthat::expect_match(
    html, 'junction_evaluations.csv</code></a></td><td class="numeric">2</td>', fixed = TRUE
  )
  testthat::expect_match(html, "original input row IDs", fixed = TRUE)
  testthat::expect_match(html, "each applicable junction", fixed = TRUE)
  testthat::expect_match(html, "Open with readRDS() in R", fixed = TRUE)
  testthat::expect_match(html, "Row counts exclude CSV headers", fixed = TRUE)
})

testthat::test_that("unavailable output files are explained without broken download links", {
  missing_path <- tempfile("missing-coverage-", fileext = ".csv")
  html <- output_manifest_html(
    publication_result, list(coverage_summary = missing_path)
  )
  testthat::expect_match(html, "File not available", fixed = TRUE)
  testthat::expect_false(grepl("<a href=", html, fixed = TRUE))
  testthat::expect_match(
    output_manifest_html(publication_result, list()),
    "No saved data files were supplied", fixed = TRUE
  )
})

testthat::test_that("the report groups output files in a closed manifest disclosure", {
  report <- write_fusion_report(
    publication_result, tempfile("manifest-report-"), list()
  )
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(
    html,
    '<details id="files"><summary>Output manifest: file contents and row counts</summary>',
    fixed = TRUE
  )
  testthat::expect_false(grepl("Machine-readable artifacts", html, fixed = TRUE))
  audit_start <- regexpr('<section id="audit">', html, fixed = TRUE)[[1L]]
  manifest_start <- regexpr('<details id="files">', html, fixed = TRUE)[[1L]]
  testthat::expect_gt(manifest_start, audit_start)
  testthat::expect_match(html, 'href="#overview"', fixed = TRUE)
})
