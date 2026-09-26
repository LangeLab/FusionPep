# FusionPep

FusionPep maps peptide sequences supplied directly or from peptide-spectrum match (PSM) exports to a proposed fusion protein and its two selected parent proteins. It reports mapped occurrences, parent sharing, residue coverage, and pairwise fusion-to-parent alignments. When junction coordinates are supplied, it checks whether an occurrence crosses a junction and meets the required minimum flank length on both sides. FusionPep analyzes sequence matches; it does not search spectra, validate PSMs, or establish peptide detection, FDR-controlled evidence, or biological fusion expression.

Absence from the two supplied parent sequences does not show that a peptide is unique across the proteome.

## Quick start

From the project root, with R 4.6.0 or newer and `renv` available, run the reference example:

```bash
Rscript setup_renv.R
Rscript run_fusion_mapper.R --output=results/example
```

Setup installs dependencies into the project-local library and writes a lockfile snapshot. It requires network access; compiling dependencies may require a system build toolchain. PDF figure export requires Cairo graphics. See [Getting started](wiki/Getting-Started.md) for prerequisites, setup details, commands, and path behavior.

Open `results/example/fusion_peptide_mapper_report.html`. Keep the output directory together when sharing linked PDFs and data files. Reusing an output directory replaces generated files with the same names.

The bundled PML::RARA example uses reference peptide sequences and controls, not PSMs measured by this project. The [worked example](wiki/Worked-Example.md) explains its inputs, decisions, and coverage. [Input source notes](input/README.md) record the verified sources and the longer peptide's unconfirmed literature attribution.

## Documentation

Start with the [FusionPep Wiki home](wiki/Home.md). The `wiki/` directory contains the GitHub Wiki Markdown sources, including its sidebar and footer.

- [Getting started](wiki/Getting-Started.md): prerequisites, setup, commands, options, and path behavior.
- [Input formats](wiki/Input-Formats.md): protein records, peptide tables, and junction definitions.
- [Interpreting results](wiki/Interpreting-Results.md): mapping classes, junction decisions, coverage, and scientific limits.
- [Reports and figures](wiki/Reports-and-Figures.md): report navigation, figures, printing, and saved outputs.
- [R interface](wiki/R-Interface.md): reusable functions and output writing.
- [Troubleshooting](wiki/Troubleshooting.md): setup, input, and interpretation checks.
- [Architecture](wiki/Architecture.md): analysis flow and source organization.

For your own analysis, supply paths for the protein FASTA, peptide table, and junction table. The CLI retains bundled example files when an input path is omitted; see [Input formats](wiki/Input-Formats.md) before using those defaults.

## Citation

The root [CITATION.cff](CITATION.cff) records the software citation. When this file is on a GitHub repository's default branch, GitHub's ["Cite this repository" panel](https://docs.github.com/en/enterprise-cloud%40latest/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files) can export the citation in APA or BibTeX. Use the BibTeX entry in an R Markdown or Quarto bibliography. FusionPep is an R project, not an installable R package, so `citation()` without a package name returns the citation for R itself.

## Authors

Enes K. Ergin, Agustina Conrrero, and Lange Lab.

## Development

Run the tests with the existing local environment active:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

`Rscript check_project.R` also regenerates files in `results/`. Preserve any earlier run before using it. For a separate analysis, choose a new directory with `--output=PATH`. [Architecture](wiki/Architecture.md) explains the source files and development checks.

## License

FusionPep code and documentation use the [MIT license](LICENSE). The bundled reference sequences retain the source attribution in [input source notes](input/README.md).
