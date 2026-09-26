# Input formats

A FusionPep analysis selects one fusion protein and two parental proteins, reads a peptide CSV, and uses an explicit junction CSV when junction assessment is required. Coordinates refer to amino-acid residues, starting at one and including both interval endpoints.

## Protein FASTA

For the command-line defaults, use these header identifiers followed by the corresponding protein sequences:

```text
>Fusion
>ParentA
>ParentB
```

These are header examples, not a complete FASTA. Each record needs a non-empty, ungapped protein sequence. Descriptions after the first whitespace are allowed, so accessions, isoforms, and source information can remain in the header. Each selected identifier must match exactly one record, and the three roles must select different records.

FusionPep selects the configured records if a FASTA contains additional entries. It does not compare peptides against those extra records. Use the [R interface](R-Interface) to map other header identifiers onto the three roles without editing the runner.

Sequence text is converted to uppercase and whitespace is removed. Standard amino-acid symbols and `B`, `X`, `Z`, `J`, `U`, and `O` are accepted. Ambiguity symbols are matched literally, not as wildcards. Gaps, stop symbols, and unsupported characters are rejected.

Verify that the file contains protein sequences. A string made only of `A`, `C`, `G`, and `T` cannot be identified as DNA from its letters alone because all four are also amino-acid symbols.

Keep the source accession or equivalent identifier, record version or release when available, isoform identity, and the reason for choosing each reference with the input files. FusionPep records the source headers and input file hashes but does not retrieve or verify those references remotely.

## Peptide CSV

The default required column is `peptide`. This row is taken from the bundled example:

```csv
peptide,evidence_role,reference
ITQGK,parent_control,UniProt PML P29590-derived sequence
```

Use `--peptide-column=NAME` when the sequence column has a different name. Column lookup accepts an exact name or a unique case-insensitive match. Other columns remain in `input_peptides.csv` and the R result.

Retain relevant upstream metadata such as search engine, score, q-value, target/decoy status, sample identifier, and evidence role. FusionPep does not apply confidence thresholds or use those columns to include or exclude rows. If you need a filtered peptide set, select it before running the mapper and record that selection.

Controls and reference peptides also contribute to coverage when they are valid input sequences. `evidence_role` explains a row; it is not a coverage filter.

### Normalization and row identity

The normalizer handles whitespace, letter case, terminal underscores, flanking notation such as `K.PEPTIDE.R`, and supported modification annotations in brackets, parentheses, or braces. Supported annotations include recognized modification names, selected UniMod tokens, numeric mass shifts, and recognized isotope or reporter labels. Unsupported or malformed annotations are reported and excluded from matching.

The original input, normalized sequence, matching sequence, row identifier, status, and normalization notes are retained. Do not replace an excluded value with a guessed peptide just to obtain a match.

By default, I and L are equivalent for matching. `--exact-il` keeps them distinct. This setting can change parent sharing and junction classification. Figure labels use the matched reference subsequence, which can differ from an input string under I/L equivalence.

Rows that produce the same matching sequence share a mapping result with the contributing input row IDs recorded. Duplicate rows do not add new covered residues, and distinct occurrences of a peptide are retained.

## Junction CSV

Junctions are supplied as coordinates in the selected fusion sequence. Multiple rows can describe different junctions in that fusion; each `junction_id` must be unique.

The required fields are:

| Field                   | Meaning                                                         |
| ----------------------- | --------------------------------------------------------------- |
| `fusion_id`             | Canonical role `Fusion`                                         |
| `junction_id`           | Unique identifier for this junction                             |
| `upstream_parent`       | Upstream role, `ParentA` or `ParentB`                           |
| `downstream_parent`     | The other parent role                                           |
| `fusion_left_position`  | Last upstream-parent residue in fusion coordinates              |
| `fusion_right_position` | First downstream-parent residue in fusion coordinates           |
| `min_flank_aa`          | Required residue count on each side, an integer of at least one |

Additional fields describe the supplied reference and join:

| Field                        | Meaning                                            |
| ---------------------------- | -------------------------------------------------- |
| `upstream_parent_position`   | Corresponding residue in the upstream parent       |
| `downstream_parent_position` | Corresponding residue in the downstream parent     |
| `inserted_sequence`          | Residues between the two fusion boundary positions |
| `source`                     | Source record or reference                         |
| `notes`                      | Explanation of the supplied junction               |

Supply both parent coordinates or neither; leave both cells blank or write `NA` to omit them. Elsewhere in the junction CSV, `NA` is literal text: a `junction_id` of `NA` is an identifier and an `inserted_sequence` of `NA` is asparagine followed by alanine. Identifiers are kept exactly as written, including leading zeros. When supplied, the boundary residues must agree exactly with those parent residues. This checks local consistency; it does not establish ancestry or breakpoint mechanism.

For a direct join, the right position is the left position plus one and `inserted_sequence` is empty. If there are intervening residues, provide their exact sequence. Its length must equal the right position minus the left position minus one. All coordinates must fit within the relevant sequences.

The bundled reference uses the following row:

```csv
fusion_id,junction_id,upstream_parent,downstream_parent,fusion_left_position,fusion_right_position,upstream_parent_position,downstream_parent_position,inserted_sequence,min_flank_aa,source,notes
Fusion,PML_RARA_bcr3,ParentA,ParentB,394,396,394,61,A,3,AAA60126.1,"PML residue 394 is followed by one breakpoint-derived alanine and RARA residue 61"
```

The value three is the bundled example's flank requirement, not a universal threshold or an implicit default for a missing CSV field. See [Interpreting results](Interpreting-Results) for how both flank lengths are calculated.

The reusable R analysis allows missing junction metadata and reports an unassessed junction state. An empty supplied junction CSV is rejected. The CLI defaults to the example junction file, so use explicit paths for your own analysis.
