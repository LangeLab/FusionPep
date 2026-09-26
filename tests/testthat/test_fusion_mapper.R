testthat::test_that("normalization handles supported modifications and flanks", {
  normalized <- normalize_one_peptide(
    "K.M[Oxidation (M)]PEPTIDE.R",
    il_equivalent = TRUE
  )
  testthat::expect_identical(normalized$normalization_status, "ok")
  testthat::expect_identical(normalized$normalized_peptide, "MPEPTIDE")
  testthat::expect_identical(normalized$matching_peptide, "MPEPTLDE")
  testthat::expect_match(normalized$normalization_notes, "flanking residues")

  unsupported <- normalize_one_peptide("K.PEP[UnknownModification].R")
  testthat::expect_identical(unsupported$normalization_status, "unsupported_annotation")
  starred <- normalize_one_peptide("PEP*")
  testthat::expect_identical(starred$normalization_status, "unsupported_annotation")
  unimod <- normalize_one_peptide("PEP[UniMod:35]")
  testthat::expect_identical(unimod$normalization_status, "ok")
  unimod_name <- normalize_one_peptide("PEP[UniMod:Oxidation]")
  testthat::expect_identical(unimod_name$normalization_status, "ok")
  unknown_unimod <- normalize_one_peptide("PEP[UniMod:999999]")
  testthat::expect_identical(
    unknown_unimod$normalization_status,
    "unsupported_annotation"
  )
  unknown_label <- normalize_one_peptide("PEP[unknown label]")
  testthat::expect_identical(
    unknown_label$normalization_status,
    "unsupported_annotation"
  )
})

testthat::test_that("normalization caches duplicate raw peptide values", {
  normalized <- normalize_peptides(c("PEPTIDE", "PEPTIDE", "PEP[+15.99]"))
  testthat::expect_identical(normalized$input_row_id, 1:3)
  testthat::expect_identical(
    normalized$normalized_peptide[1:2],
    c("PEPTIDE", "PEPTIDE")
  )
  testthat::expect_identical(
    normalized$normalization_status,
    c("ok", "ok", "ok")
  )
  testthat::expect_identical(
    normalized$normalization_notes[1],
    normalized$normalization_notes[2]
  )
})

testthat::test_that("literal NA remains a valid peptide value", {
  peptide_file <- tempfile(fileext = ".csv")
  writeLines(c("peptide", "NA"), peptide_file)
  input <- read_peptide_table(peptide_file, il_equivalent = FALSE)
  testthat::expect_identical(input$normalized$peptide_input, "NA")
  testthat::expect_identical(input$normalized$normalization_status, "ok")
})

testthat::test_that("overlapping peptide occurrences are all mapped", {
  peptides <- normalize_peptides("AAA", il_equivalent = FALSE)
  hits <- find_peptide_hits("AAAAA", "Fusion", peptides, il_equivalent = FALSE)
  testthat::expect_equal(hits$start, c(1L, 2L, 3L))
  testthat::expect_equal(hits$end, c(3L, 4L, 5L))
  testthat::expect_false(any(
    c("matched_subsequence_normalized", "junction_spanning", "junction_ids") %in%
      names(hits)
  ))
})

testthat::test_that("same-width dictionary matching preserves peptide coordinates", {
  peptides <- normalize_peptides(
    c("AAB", "ABB", "BBB", "BBC", "CAA"),
    il_equivalent = FALSE
  )
  hits <- find_peptide_hits(
    "AAABBBCAAABBB",
    "Fusion",
    peptides,
    il_equivalent = FALSE
  )
  starts <- split(hits$start, hits$matching_peptide)
  testthat::expect_equal(starts[["AAB"]], c(2L, 9L))
  testthat::expect_equal(starts[["ABB"]], c(3L, 10L))
  testthat::expect_equal(starts[["BBB"]], c(4L, 11L))
  testthat::expect_equal(starts[["BBC"]], 5L)
  testthat::expect_equal(starts[["CAA"]], 7L)
})

testthat::test_that("coverage is calculated as a union of intervals", {
  hits <- data.frame(
    sequence_id = c("Fusion", "Fusion"),
    matching_peptide = c("AAAA", "BBBB"),
    start = c(1L, 3L),
    end = c(4L, 6L),
    stringsAsFactors = FALSE
  )
  coverage <- calculate_coverage(
    c(Fusion = "AAAAAAA"),
    "Fusion",
    hits
  )
  testthat::expect_equal(coverage$summary$covered_residues, 6L)
  testthat::expect_equal(coverage$summary$coverage_percent, 85.71)
  testthat::expect_equal(coverage$uncovered_regions$start, 7L)
  testthat::expect_equal(coverage$uncovered_regions$end, 7L)
})

testthat::test_that("gap coordinates remain explicit", {
  testthat::expect_equal(
    make_position_vector_from_chars(
      strsplit("ABC-DE", "", fixed = TRUE)[[1]],
      start_position = 1L
    ),
    c(1L, 2L, 3L, NA_integer_, 4L, 5L)
  )
  testthat::expect_equal(
    make_position_vector_from_chars(
      strsplit("--ABC", "", fixed = TRUE)[[1]],
      start_position = 4L
    ),
    c(NA_integer_, NA_integer_, 4L, 5L, 6L)
  )
})

testthat::test_that("empty alignment details name the requested alignment type", {
  local_message <- alignment_block_html(
    empty_alignment_columns(),
    "ParentA",
    alignment_type = "local"
  )
  testthat::expect_match(local_message, "No local alignment columns")
})

testthat::test_that("empty local alignments remain valid results", {
  alignment <- run_pairwise_alignment(
    "AAAA",
    "RRRR",
    "ParentA",
    store_alignment_objects = FALSE
  )
  local_summary <- alignment$summaries[
    alignment$summaries$alignment_type == "local",
    ,
    drop = FALSE
  ]
  testthat::expect_equal(nrow(alignment$columns[
    alignment$columns$alignment_type == "local",
    ,
    drop = FALSE
  ]), 0L)
  testthat::expect_equal(local_summary$aligned_residues, 0L)
  testthat::expect_equal(local_summary$fusion_aligned_fraction, 0)
})

testthat::test_that("alignment parameters require finite non-negative values", {
  testthat::expect_error(
    run_pairwise_alignment(
      "AAAA",
      "AAAA",
      "ParentA",
      gap_opening = -1,
      store_alignment_objects = FALSE
    ),
    "gap_opening"
  )
  testthat::expect_error(
    run_pairwise_alignment(
      "AAAA",
      "AAAA",
      "ParentA",
      gap_extension = Inf,
      store_alignment_objects = FALSE
    ),
    "gap_extension"
  )
  testthat::expect_error(
    alignment_overview_bins(empty_alignment_columns(), bin_width = 1.5),
    "positive integer"
  )
})

testthat::test_that("junction annotations require both-sided flanks", {
  junction_file <- tempfile(fileext = ".csv")
  writeLines(c(
    "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,upstream_parent_position,downstream_parent_position,min_flank_aa,source,notes",
    "Fusion,j1,ParentA,ParentB,4,5,4,2,2,test,synthetic"
  ), junction_file)
  sequence_text <- c(
    Fusion = "AAAXBBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  junctions <- read_fusion_junctions(junction_file, sequence_text)
  normalized <- normalize_peptides(c("AXBB", "XBB"), il_equivalent = FALSE)
  hits <- find_peptide_hits(
    sequence_text[["Fusion"]],
    "Fusion",
    normalized,
    il_equivalent = FALSE
  )
  annotated <- annotate_peptide_hits_with_junctions(hits, junctions$table)
  passing <- annotated[annotated$matching_peptide == "AXBB", , drop = FALSE]
  short <- annotated[annotated$matching_peptide == "XBB", , drop = FALSE]
  testthat::expect_true(passing$crosses_junction)
  testthat::expect_true(passing$junction_flank_pass)
  testthat::expect_equal(passing$junction_left_flank, 2L)
  testthat::expect_equal(passing$junction_right_flank, 2L)
  testthat::expect_true(short$crosses_junction)
  testthat::expect_false(short$junction_flank_pass)
})

testthat::test_that("multiple junctions retain complete flank evaluations", {
  hits <- data.frame(
    sequence_id = "Fusion",
    matching_peptide = "PEPTIDE",
    normalized_peptide = "PEPTIDE",
    peptide_input_values = "PEPTIDE",
    input_row_ids = "1",
    input_count = 1L,
    start = 1L,
    end = 12L,
    length = 12L,
    matched_subsequence = "PEPTIDE",
    match_basis = "exact",
    stringsAsFactors = FALSE
  )
  junctions <- data.frame(
    fusion_id = c("Fusion", "Fusion"),
    junction_id = c("j1", "j2"),
    fusion_left_position = c(4L, 8L),
    fusion_right_position = c(5L, 9L),
    min_flank_aa = c(6L, 6L),
    stringsAsFactors = FALSE
  )
  evaluations <- evaluate_peptide_junctions(hits, junctions)
  annotated <- apply_junction_evaluations(hits, evaluations)

  testthat::expect_equal(
    evaluations[, c(
      "junction_id",
      "junction_left_flank",
      "junction_right_flank"
    )],
    data.frame(
      junction_id = c("j1", "j2"),
      junction_left_flank = c(4L, 8L),
      junction_right_flank = c(8L, 4L),
      stringsAsFactors = FALSE
    )
  )
  testthat::expect_identical(annotated$junction_id, "j1")
  testthat::expect_equal(annotated$junction_left_flank, 4L)
  testthat::expect_equal(annotated$junction_right_flank, 8L)
  testthat::expect_false(annotated$junction_flank_pass)
})

testthat::test_that("summary flank values come from one peptide occurrence", {
  normalized <- normalize_peptides("PEPTIDE", il_equivalent = FALSE)
  hits <- data.frame(
    sequence_id = rep("Fusion", 3L),
    matching_peptide = rep("PEPTIDE", 3L),
    normalized_peptide = rep("PEPTIDE", 3L),
    peptide_input_values = rep("PEPTIDE", 3L),
    input_row_ids = rep("1", 3L),
    input_count = rep(1L, 3L),
    start = c(1L, 2L, 3L),
    end = c(7L, 8L, 9L),
    length = rep(7L, 3L),
    matched_subsequence = rep("PEPTIDE", 3L),
    match_basis = rep("exact", 3L),
    crosses_junction = rep(TRUE, 3L),
    junction_flank_pass = c(FALSE, TRUE, FALSE),
    junction_id = rep("j1", 3L),
    junction_left_flank = c(3L, 2L, 1L),
    junction_right_flank = c(1L, 2L, 3L),
    stringsAsFactors = FALSE
  )
  summary <- build_peptide_summary(normalized, hits)
  testthat::expect_equal(summary$junction_left_flank, 2L)
  testthat::expect_equal(summary$junction_right_flank, 2L)
})

testthat::test_that("mapping and alignment layers are independently reusable", {
  sequence_text <- c(
    Fusion = "AAAXYBBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  normalized <- normalize_peptides("AXY", il_equivalent = FALSE)
  mapping <- map_peptides_to_references(
    sequence_text,
    normalized,
    il_equivalent = FALSE
  )
  testthat::expect_true(any(mapping$peptide_hits$sequence_id == "Fusion"))
  testthat::expect_equal(nrow(mapping$coverage), 3L)

  alignment <- align_fusion_to_parents(
    sequence_text,
    store_alignment_objects = FALSE
  )
  testthat::expect_identical(names(alignment$alignments), c("ParentA", "ParentB"))
  testthat::expect_true(nrow(alignment$alignment_summaries) >= 2L)
})

testthat::test_that("alignment overview preserves disjoint parent intervals", {
  match_status <- c(rep(TRUE, 25L), rep(FALSE, 25L), rep(TRUE, 25L))
  columns <- data.frame(
    parent = rep("ParentA", 75L),
    alignment_type = rep("local", 75L),
    alignment_column = seq_len(75L),
    fusion_pos = seq_len(75L),
    parent_pos = seq_len(75L),
    fusion_aa = rep("A", 75L),
    parent_aa = rep("A", 75L),
    raw_match = match_status,
    il_equivalent_match = match_status,
    status = ifelse(match_status, "match", "mismatch"),
    stringsAsFactors = FALSE
  )
  summaries <- data.frame(
    parent = "ParentA",
    alignment_type = "local",
    substitution_matrix = "simple_aa_compatible",
    score = 1,
    alignment_columns = 75L,
    aligned_residues = 75L,
    raw_matches = 50L,
    il_equivalent_matches = 50L,
    raw_identity_percent = 66.67,
    il_equivalent_identity_percent = 66.67,
    fusion_aligned_fraction = 1,
    parent_aligned_fraction = 1,
    stringsAsFactors = FALSE
  )
  overview <- alignment_overview_table(list(
    alignment_columns = columns,
    alignment_summaries = summaries
  ))
  testthat::expect_identical(
    overview$high_similarity_fusion_residues,
    "1-25; 51-75"
  )
  testthat::expect_identical(
    overview$high_similarity_parent_residues,
    "1-25; 51-75"
  )
})

testthat::test_that("junction metadata cannot use Fusion as a parent", {
  junction_file <- tempfile(fileext = ".csv")
  writeLines(c(
    "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,min_flank_aa",
    "Fusion,j1,Fusion,ParentB,4,5,2"
  ), junction_file)
  sequence_text <- c(
    Fusion = "AAAXBBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  testthat::expect_error(
    read_fusion_junctions(junction_file, sequence_text),
    "parent values must refer to parent sequences"
  )
})

testthat::test_that("junction metadata cannot identify a parent as the fusion", {
  junction_file <- tempfile(fileext = ".csv")
  writeLines(c(
    "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,min_flank_aa",
    "ParentA,j1,ParentA,ParentB,4,5,2"
  ), junction_file)
  sequence_text <- c(
    Fusion = "AAAXBBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  testthat::expect_error(
    read_fusion_junctions(junction_file, sequence_text),
    "fusion_id values must be 'Fusion'"
  )
})

testthat::test_that("junction text values keep literal NA and leading zeros", {
  junction_file <- tempfile(fileext = ".csv")
  writeLines(c(
    "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,upstream_parent_position,downstream_parent_position,inserted_sequence,min_flank_aa",
    "Fusion,NA,ParentA,ParentB,3,6,NA,NA,NA,2",
    "Fusion,01,ParentA,ParentB,2,3,,,,1"
  ), junction_file)
  sequence_text <- c(
    Fusion = "AAANABBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  junctions <- read_fusion_junctions(junction_file, sequence_text)$table
  testthat::expect_identical(junctions$junction_id, c("NA", "01"))
  testthat::expect_identical(junctions$inserted_sequence, c("NA", ""))
  testthat::expect_identical(junctions$upstream_parent_position, c(NA_integer_, NA_integer_))
  testthat::expect_identical(junctions$fusion_left_position, c(3L, 2L))
})

testthat::test_that("junction coordinates must be finite integers in range", {
  sequence_text <- c(
    Fusion = "AAAXBBB",
    ParentA = "AAAX",
    ParentB = "YBBB"
  )
  for (bad in c("Inf", "4.5", "3e10", "NA", "")) {
    junction_file <- tempfile(fileext = ".csv")
    writeLines(c(
      "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,min_flank_aa",
      paste0("Fusion,j1,ParentA,ParentB,", bad, ",5,2")
    ), junction_file)
    testthat::expect_error(
      read_fusion_junctions(junction_file, sequence_text),
      "must contain integer values",
      info = bad
    )
  }
})

testthat::test_that("FASTA validation rejects missing required records", {
  fasta <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">Fusion",
    "MPEPTIDE",
    ">ParentA",
    "MPEPTIDE"
  ), fasta)
  testthat::expect_error(
    read_fusion_sequences(fasta),
    "Could not find FASTA record"
  )
})

testthat::test_that("FASTA invalid residues are rejected before parsing", {
  fasta <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">Fusion",
    "ACD?EF",
    ">ParentA",
    "ACDEFG",
    ">ParentB",
    "ACDEFG"
  ), fasta)
  testthat::expect_error(
    read_fusion_sequences(fasta),
    "contains unsupported character"
  )
})

testthat::test_that("end-to-end mapping classifies peptide sequence matches", {
  fasta <- tempfile(fileext = ".fasta")
  peptides <- tempfile(fileext = ".csv")
  writeLines(c(
    ">Fusion demo",
    "MPEPTIDEALPHABRAVCDEFGHIKLMNPQRSTVWY",
    ">ParentA demo",
    "MPEPTIDEALPHAFUSIONAQRKLMNQRSTVWYACDEFGHIKL",
    ">ParentB demo",
    "MGKPEPTIDEBETABRAVCDEFGHIKLMNPQRSTVWY"
  ), fasta)
  writeLines(c(
    "peptide",
    "PEPTIDE",
    "ALPHA",
    "BRAVCDE",
    "HAB",
    "FUSIONA",
    "BETAB"
  ), peptides)

  result <- run_fusion_analysis(
    sequence_file = fasta,
    peptide_file = peptides,
    il_equivalent = TRUE
  )
  classifications <- setNames(
    result$peptide_summary$presence_classification,
    result$peptide_summary$normalized_peptide
  )
  testthat::expect_false("classification" %in% names(result$peptide_summary))
  evidence_classes <- setNames(
    result$peptide_summary$fusion_evidence_class,
    result$peptide_summary$normalized_peptide
  )
  testthat::expect_identical(classifications[["HAB"]], "fusion_only")
  testthat::expect_identical(
    evidence_classes[["HAB"]],
    "fusion_mapping_junction_unassessed"
  )
  hab_hit <- result$peptide_hits[
    result$peptide_hits$matching_peptide == "HAB" &
      result$peptide_hits$sequence_id == "Fusion",
    ,
    drop = FALSE
  ]
  testthat::expect_identical(hab_hit$display_class, "fusion_mapping_junction_unassessed")
  testthat::expect_identical(classifications[["FUSIONA"]], "parentA_only")
  testthat::expect_identical(classifications[["BETAB"]], "parentB_only")
  testthat::expect_identical(classifications[["PEPTIDE"]], "shared_all_three")
  testthat::expect_true(all(result$alignment_summaries$alignment_type %in% c("global", "local")))
  testthat::expect_true(any(
    result$alignment_summaries$substitution_matrix == "simple_aa_compatible"
  ))
  testthat::expect_true(any(result$alignment_columns$status == "match"))
  alignment_bins <- alignment_overview_bins(result$alignment_columns)
  testthat::expect_true(any(
    alignment_bins$parent == "ParentA" &
      alignment_bins$status == "high_exact_identity"
  ))
  testthat::expect_true(any(
    alignment_bins$parent == "ParentB" &
      alignment_bins$status == "high_exact_identity"
  ))
  testthat::expect_null(result$alignments$ParentA$global)
  testthat::expect_null(result$alignments$ParentB$local)
  testthat::expect_match(
    summary_cards_html(result),
    "not assessed; no fusion-junction metadata supplied"
  )
  testthat::expect_silent(
    ggplot2::ggplot_build(
      plot_peptide_coverage(
        result$peptide_hits,
        result$sequence_input$sequence_text,
        result$junctions
      )
    )
  )
  no_junction_report <- write_fusion_report(
    result,
    tempfile("fusion-no-junction-report-"),
    list()
  )
  no_junction_report_text <- paste(
    readLines(no_junction_report, warn = FALSE),
    collapse = "\n"
  )
  testthat::expect_match(
    no_junction_report_text,
    "Junction status is not assessed"
  )
})

testthat::test_that("HTML report and machine-readable outputs are generated", {
  fasta <- file.path(project_root, "input", "sequences.fasta")
  peptides <- file.path(project_root, "input", "peptides.csv")
  junctions <- file.path(project_root, "input", "fusion_junctions.csv")
  result <- run_fusion_analysis(
    sequence_file = fasta,
    peptide_file = peptides,
    junction_file = junctions
  )
  output_dir <- tempfile("fusion-output-")
  artifacts <- write_fusion_outputs(result, output_dir)
  report <- write_fusion_report(result, output_dir, artifacts)

  testthat::expect_true(file.exists(report))
  testthat::expect_true(file.exists(artifacts$peptide_hits))
  testthat::expect_true(file.exists(artifacts$fusion_junctions))
  testthat::expect_true(file.exists(artifacts$junction_evaluations))
  testthat::expect_true(file.exists(artifacts$input_peptides))
  testthat::expect_true(file.exists(artifacts$peptide_coverage_plot))
  testthat::expect_true(file.exists(artifacts$alignment_status_plot))
  for (name in c(
    "peptide_coverage_pdf", "junction_evidence_pdf", "alignment_status_pdf"
  )) {
    testthat::expect_true(file.exists(artifacts[[name]]))
    testthat::expect_identical(
      rawToChar(readBin(artifacts[[name]], what = "raw", n = 4L)), "%PDF"
    )
  }
  testthat::expect_true(file.exists(artifacts$result_rds))
  reloaded <- readRDS(artifacts$result_rds)
  testthat::expect_s3_class(reloaded, "fusion_peptide_mapping_result")
  testthat::expect_identical(
    names(reloaded$junction_evaluations),
    names(result$junction_evaluations)
  )
  testthat::expect_identical(result$junctions$inserted_sequence, "A")
  testthat::expect_identical(result$junctions$fusion_left_position, 394L)
  testthat::expect_identical(result$junctions$fusion_right_position, 396L)
  junction_peptide <- result$peptide_summary[
    result$peptide_summary$peptide_input_values == "ITQGKAIETQSSSSEEIV",
    ,
    drop = FALSE
  ]
  testthat::expect_equal(nrow(junction_peptide), 1L)
  testthat::expect_true(junction_peptide$crosses_junction)
  testthat::expect_true(junction_peptide$junction_flank_pass)
  testthat::expect_identical(
    junction_peptide$fusion_evidence_class,
    "junction_spanning_candidate"
  )
  short_junction_peptide <- result$peptide_summary[
    result$peptide_summary$peptide_input_values == "LSSCITQGKAIE",
    ,
    drop = FALSE
  ]
  testthat::expect_true(short_junction_peptide$crosses_junction)
  testthat::expect_false(short_junction_peptide$junction_flank_pass)
  report_text <- paste(readLines(report, warn = FALSE), collapse = "\n")
  testthat::expect_match(report_text, "Annotated reference sequences")
  testthat::expect_match(report_text, "Junction decision")
  testthat::expect_match(report_text, "Peptide mapping tracks")
  testthat::expect_match(report_text, "valid input rows")
  testthat::expect_match(
    report_text,
    "supplied[[:space:]]+sequence-junction metadata"
  )
})

testthat::test_that("large report tables state when rows are truncated", {
  html <- table_to_html(
    data.frame(
      peptide = paste0("PEP", 1:3),
      score = 1:3,
      stringsAsFactors = FALSE
    ),
    max_rows = 2L
  )
  testthat::expect_match(html, "Showing the first 2 of 3 rows")
  testthat::expect_match(
    html,
    'Complete tables are described in the <a href="#files">output manifest</a>',
    fixed = TRUE
  )
})

testthat::test_that("missing evidence roles are rendered as unspecified", {
  result <- list(
    peptide_input = list(
      input = data.frame(
        peptide = "PEPTIDE",
        evidence_role = NA_character_,
        stringsAsFactors = FALSE
      )
    )
  )
  summary <- data.frame(input_row_ids = "1", stringsAsFactors = FALSE)
  roles <- input_roles_for_summary(result, summary)
  testthat::expect_identical(unname(roles), "")
  testthat::expect_identical(unname(friendly_input_role(roles)), "Unspecified")
})

testthat::test_that("saved outputs record input file names, not directories", {
  input_dir <- tempfile("fusionpep-private-inputs-")
  dir.create(input_dir)
  file.copy(
    file.path(project_root, "input", c("sequences.fasta", "peptides.csv", "fusion_junctions.csv")),
    input_dir
  )
  result <- run_fusion_analysis(
    sequence_file = file.path(input_dir, "sequences.fasta"),
    peptide_file = file.path(input_dir, "peptides.csv"),
    junction_file = file.path(input_dir, "fusion_junctions.csv")
  )
  output_dir <- tempfile("fusion-output-")
  artifacts <- write_fusion_outputs(result, output_dir)
  artifacts$report <- write_fusion_report(result, output_dir, artifacts)

  manifest <- setNames(result$manifest$value, result$manifest$key)
  testthat::expect_identical(
    unname(manifest[c("fasta_file", "peptide_csv_file", "junction_csv_file")]),
    c("sequences.fasta", "peptides.csv", "fusion_junctions.csv")
  )
  input_dir_name <- basename(input_dir)
  text_files <- c(
    unlist(artifacts[grepl("[.](csv|txt|html)$", unlist(artifacts))])
  )
  for (path in text_files) {
    testthat::expect_false(
      any(grepl(input_dir_name, readLines(path, warn = FALSE), fixed = TRUE)),
      info = path
    )
  }
  saved <- readRDS(artifacts$result_rds)
  saved_strings <- unlist(rapply(
    unclass(saved), function(x) x, classes = "character", how = "unlist"
  ))
  testthat::expect_false(any(grepl(input_dir_name, saved_strings, fixed = TRUE)))
  testthat::expect_identical(saved$sequence_input$source_path, "sequences.fasta")
  testthat::expect_identical(saved$config$junction_file, "fusion_junctions.csv")
  testthat::expect_identical(
    result$sequence_input$source_path,
    normalizePath(file.path(input_dir, "sequences.fasta"))
  )
})

testthat::test_that("outputs cannot overwrite input files", {
  input_dir <- tempfile("fusion-input-")
  dir.create(input_dir)
  fasta <- file.path(input_dir, "sequences.fasta")
  peptides <- file.path(input_dir, "peptides.csv")
  junctions <- file.path(input_dir, "fusion_junctions.csv")
  writeLines(c(
    ">Fusion",
    "AAAXYBBB",
    ">ParentA",
    "AAAX",
    ">ParentB",
    "YBBB"
  ), fasta)
  writeLines(c("peptide", "AXY"), peptides)
  writeLines(c(
    "fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,min_flank_aa",
    "Fusion,j1,ParentA,ParentB,4,5,1"
  ), junctions)
  result <- run_fusion_analysis(
    sequence_file = fasta,
    peptide_file = peptides,
    junction_file = junctions,
    il_equivalent = FALSE
  )
  testthat::expect_error(
    write_fusion_outputs(result, input_dir),
    "overwrite input file"
  )
})
