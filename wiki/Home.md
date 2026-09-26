# FusionPep

**Peptide mapping and junction review for fusion proteins.**

FusionPep maps supplied peptides onto a proposed fusion protein and its two supplied parents. It reports sequence positions, coverage, parent comparisons, and whether mapped occurrences cross explicitly defined junctions. An HTML report brings the decisions, alignments, and figures together for review.

Use it when you already have the protein references and peptide sequences you want to inspect. A peptide list or a peptide-spectrum match (PSM) export can supply the sequences; FusionPep retains the accompanying columns but does not perform spectrum searching or confidence filtering.

## Start here

1. [Get started](Getting-Started): Set up the R environment and run the bundled example.
1. [Prepare your inputs](Input-Formats): Supply protein references, a peptide table, and junction coordinates.
1. [Interpret the results](Interpreting-Results): Distinguish parent sharing, junction crossing, and the flank rule.
1. [Use the report and figures](Reports-and-Figures): Inspect decisions, export figures, and find the complete data tables.

The [worked example](Worked-Example) follows two supplied PML::RARA fusion-region peptides through their junction assessments. It includes actual results from the bundled reference inputs.

## What a run gives you

- Every mapped occurrence in the three selected reference sequences, with one-based inclusive coordinates.
- Junction assessments that retain both flank lengths and the configured threshold.
- Residue coverage and global and local fusion-to-parent alignments.
- An HTML report, PNG and vector PDF figures, CSV tables, an R result object, and input/settings records.

The report is the starting point. Its collapsible manifests explain the saved files and the settings used for that run.

## What the results establish

A sequence match establishes that a normalized peptide occurs in a supplied reference under the selected matching settings. A junction candidate additionally meets an explicit flank rule and is absent from the two supplied parents. Neither result establishes peptide detection, PSM validity, statistical confidence, or biological fusion expression.

The comparison covers the supplied references. It does not establish proteome-wide uniqueness. Read [Interpreting results](Interpreting-Results) before using a classification in a biological claim.

## Use FusionPep from R

The [R interface](R-Interface) separates analysis from file writing and lets you reuse the mapping and alignment tables in scripts. [Architecture](Architecture) explains the current code organization. [Troubleshooting](Troubleshooting) covers common input and environment checks.
