# Interpreting results

Read the sequence-presence result separately from the junction assessment. Presence describes where a normalized peptide occurs. Junction assessment asks whether an occurrence contains both supplied boundary residues and enough residues on each side.

All conclusions refer to the selected protein records and matching settings. A result from FusionPep is sequence evidence, not a new spectrum identification or a biological confirmation.

## Sequence presence

`peptide_summary.csv` has one row per matching peptide sequence, with the contributing input row IDs. Its `presence_classification` values describe the three supplied references:

| Value | Sequence occurrence |
| --- | --- |
| `fusion_only` | Fusion only |
| `fusion_and_parentA` | Fusion and ParentA |
| `fusion_and_parentB` | Fusion and ParentB |
| `shared_all_three` | Fusion, ParentA, and ParentB |
| `parentA_only` | ParentA only |
| `parentB_only` | ParentB only |
| `parental_only_both` | Both parents, without a fusion match |
| `not_found` | None of the three references |

`fusion_only` does not establish a junction-spanning match. A peptide can occur only in the supplied fusion without crossing its defined junction. Likewise, `parentA_only` does not establish that a biological fusion lacks a parent domain; it describes this peptide and these references.

## Junction crossing and flanks

For an occurrence from residue `start` through `end`, crossing requires `start <= fusion_left_position` and `end >= fusion_right_position`. Both boundary residues must be present.

The two flank lengths include the boundary residues. The left flank is the left boundary position minus the occurrence start plus one. The right flank is the occurrence end minus the right boundary position plus one. Inserted residues between the boundaries count toward neither flank.

The occurrence passes the flank rule only when both lengths meet that junction's `min_flank_aa`. A long flank on one side does not compensate for a short flank on the other. This configurable rule is a screening heuristic, not a probability or FDR estimate.

The [worked example](Worked-Example) compares a peptide with flanks of 5 and 12 residues against one with flanks of 9 and 2 residues, using a minimum of three on each side.

## Junction evidence classes

The current `fusion_evidence_class` in `peptide_summary.csv` uses these values:

| Value | Meaning |
| --- | --- |
| `junction_spanning_candidate` | Crosses, passes both flanks, and is absent from both supplied parents |
| `junction_spanning_but_shared_with_supplied_parent` | Crosses and passes, but also occurs in a supplied parent |
| `junction_crossing_below_flank_threshold` | Crosses, but at least one flank is too short |
| `fusion_mapping_not_junction_spanning` | Maps to the fusion without crossing a supplied junction |
| `fusion_mapping_junction_unassessed` | Maps to the fusion without a supplied junction definition |
| `not_in_fusion` | Has no mapped occurrence in the supplied fusion |

No class establishes an observed or validated peptide. In particular, absence from the two supplied parents is not proteome-wide uniqueness.

The occurrence table uses `display_class` for plotting. Its shared-parent junction value is `junction_spanning_shared_parent`; the longer value above is the peptide summary label. Consult `presence_classification` for parent sharing among other fusion matches.

## Repeated occurrences and multiple junctions

`peptide_hits.csv` retains every mapped occurrence, including overlapping matches. `junction_evaluations.csv` retains one row for each crossing occurrence and applicable supplied junction, with both flank lengths and that junction's threshold. Non-crossing occurrences have no row in the evaluation table.

Peptide summaries choose the strongest junction interpretation across that peptide's occurrences. They do not replace the occurrence and evaluation tables. Use `hit_row_id` in `junction_evaluations.csv` to identify the one-based row in `peptide_hits.csv`; then use `input_row_ids` to return to the original peptide table. CSV row numbers here exclude the header.

When a peptide crosses more than one junction, inspect each evaluation independently. Do not combine the left flank from one junction with the right flank or threshold from another.

## Missing metadata and excluded peptides

Missing junction metadata means unassessed. It must not be interpreted as evidence that a peptide fails to cross a junction. The R API can produce this state explicitly with `junction_file = NULL`.

An excluded input sequence differs from a valid peptide with no match. Read `normalized_peptides.csv` and the report's input-normalization panel to identify missing, blank, unsupported, or invalid inputs and the accompanying notes.

With I/L-equivalent matching, peptide rows differing only by I and L can share the same matching sequence. Review `match_basis` and `matched_subsequence` in the occurrence table when exact letters matter. Use `--exact-il` for a new analysis that keeps them distinct.

## Coverage and alignments

Coverage counts the union of mapped residue intervals on each selected reference. Overlaps and duplicate input rows do not inflate the covered-residue count. All valid input sequences contribute, including controls; evidence-role and confidence columns do not filter coverage.

Global and local alignments compare the supplied fusion with each parent. Gap-aware alignment columns preserve the respective residue coordinates. The alignment overview is a visual summary of sequence similarity, not evidence of parent ancestry or domain loss. See [Reports and figures](Reports-and-Figures) for the overview thresholds and display limits.

## Using a result in a scientific claim

State the reference records, matching mode, normalization choices, input selection, and flank requirement. Describe a passing result as a sequence candidate under those settings.

PSM validity, FDR-controlled peptide evidence, protein inference, and biological fusion expression require evidence from the upstream experiment or further validation. The report cannot supply that evidence from peptide strings alone.
