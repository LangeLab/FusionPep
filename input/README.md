# FusionPep reference example

The default FASTA is a real human PML::RARA reference example, not a synthetic protein:

- `ParentA` is canonical human PML, UniProt `P29590-1`, sequence version 3, 882 aa.
- `ParentB` is canonical human RARA, UniProt `P10276-1`, sequence version 2, 462 aa.
- `Fusion` is the primary NCBI protein record `AAA60126.1` (PML-RAR protein), 797 aa, derived by conceptual translation from GenBank cDNA `M73779.1`.
- `fusion_junctions.csv` records the sequence-level join used by this reference: PML residue 394, one breakpoint-derived alanine, then RARA residue 61. The explicit inserted residue is important because a fusion-protein junction is not always a simple parent-protein concatenation.

The fusion record and parent records are reference sequences. The peptide CSV contains reference peptides and parent controls for demonstrating mapping and junction annotation; these rows are not observed PSMs from an MS experiment conducted by this project. Replace it with the peptide/PSM export from the actual experiment for biological interpretation.

## Peptide sources

`LSSCITQGKAIE` appears as a predicted PML-RARA peptide in Table II of [Conlon et al. (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3790285/). That entry is not marked as experimentally detected in the paper. In the bundled fusion, it crosses the supplied junction with only two downstream-parent residues, so the supplied three-residue flank rule reports it as crossing below threshold.

The CSV attributes `ITQGKAIETQSSSSEEIV` to [Gambacorti-Passerini et al. (1993)](https://doi.org/10.1182/blood.v81.5.1369.bloodjournal8151369), with a separator omitted. The article identity is confirmed, but the exact peptide transcription from its full text remains unconfirmed. Treat this row as a supplied reference peptide until that attribution is resolved. The input's `published_fusion_region_peptide` label is retained metadata, not an independent verification. Its sequence match and passing flank decision can be reproduced from the bundled files.

The remaining rows are controls drawn from the two supplied parent references. All valid input sequences, including controls, contribute to coverage.

## References

- [NCBI protein AAA60126.1](https://www.ncbi.nlm.nih.gov/protein/AAA60126.1).
- [Source cDNA M73779.1](https://www.ncbi.nlm.nih.gov/nuccore/M73779.1).
- [UniProt PML P29590](https://www.uniprot.org/uniprotkb/P29590).
- [UniProt RARA P10276](https://www.uniprot.org/uniprotkb/P10276).
- [PML::RARA breakpoint biology](https://pmc.ncbi.nlm.nih.gov/articles/PMC7139833/).
- [PML/RARA fusion-region peptide study, Blood (1993)](https://doi.org/10.1182/blood.v81.5.1369.bloodjournal8151369).
- [Original PML-RAR description, Cell (1991)](https://pubmed.ncbi.nlm.nih.gov/1652368/).
- [Predicted fusion-peptide sequences, Conlon et al. (2013)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3790285/).
- [FusionPro junction-peptide search criteria (2019)](https://pmc.ncbi.nlm.nih.gov/articles/PMC6683003/) 
- [Identification of gene fusions from human lung cancer mass spectrometry data](https://pmc.ncbi.nlm.nih.gov/articles/PMC4042237/).
