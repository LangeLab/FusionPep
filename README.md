<!-- markdownlint-disable MD033 MD041 -->
<h1 align="center">FusionPep</h1>

<p align="center">
  Peptide mapping and junction review for fusion proteins.
</p>

<p align="center">
  <img src="https://img.shields.io/badge/R-4.6%2B-2D7D46?style=flat-square&logo=r&logoColor=white" alt="R 4.6+">
  <img src="https://img.shields.io/badge/Bioconductor-3.23-1A81C2?style=flat-square" alt="Bioconductor 3.23">
  <img src="https://img.shields.io/badge/status-pre--release-C17D10?style=flat-square" alt="Pre-release">
  <a href="https://github.com/LangeLab/FusionPep/actions/workflows/ci.yml"><img src="https://img.shields.io/github/actions/workflow/status/LangeLab/FusionPep/ci.yml?branch=main&style=flat-square&logo=github&label=CI" alt="CI"></a>
  <a href="https://codecov.io/gh/LangeLab/FusionPep"><img src="https://img.shields.io/codecov/c/github/LangeLab/FusionPep?branch=main&style=flat-square&logo=codecov&logoColor=white" alt="Coverage"></a>
  <img src="https://img.shields.io/badge/license-MIT-4B9D6E?style=flat-square" alt="MIT">
</p>

<p align="center">
  <a href="NEWS.md"><img src="https://img.shields.io/badge/changelog-NEWS-E05D44?style=flat-square" alt="Changelog"></a>
  <a href="CITATION.cff"><img src="https://img.shields.io/badge/cite-CITATION.cff-0066CC?style=flat-square" alt="Citation"></a>
  <a href="https://github.com/LangeLab/FusionPep/wiki"><img src="https://img.shields.io/badge/docs-Wiki-0F766E?style=flat-square" alt="Docs"></a>
</p>

FusionPep maps peptide sequences, supplied directly or from a peptide-spectrum match (PSM) export, to a proposed fusion protein and its two parent proteins. It reports every mapped occurrence, parent sharing, residue coverage, and global and local fusion-to-parent alignments. When junction coordinates are supplied, it checks whether each occurrence crosses a junction with the required number of residues on both sides. The results come together in an HTML report, with CSV tables, an R result object, and publication-sized PNG and PDF figures.

FusionPep analyzes sequence matches. It does not search spectra, validate PSMs, or establish peptide detection, FDR-controlled evidence, or biological fusion expression. Absence from the two supplied parents does not show that a peptide is unique across the proteome.

## Quick start

FusionPep needs R 4.6.0 or newer with `renv` installed. Setup installs the locked dependency versions into a project-local library:

```bash
git clone https://github.com/LangeLab/FusionPep.git
cd FusionPep
Rscript setup_renv.R
Rscript run_fusion_mapper.R --output=results/example
```

Open `results/example/fusion_peptide_mapper_report.html`. The bundled PML::RARA example uses literature reference peptides and controls, not PSMs measured by this project. The [worked example](https://github.com/LangeLab/FusionPep/wiki/Worked-Example) explains its decisions, and the [input source notes](input/README.md) record the verified sources and the longer peptide's unconfirmed literature attribution.

Setup needs network access, and compiling dependencies can need a system build toolchain; on Linux, also install the zlib development headers (for example `zlib1g-dev`). PDF figures need an R build with Cairo graphics.

## Analyze your own data

Supply a protein FASTA with records named `Fusion`, `ParentA`, and `ParentB`, a peptide CSV with a `peptide` column, and a junction CSV:

```bash
Rscript run_fusion_mapper.R \
  --fasta=input/my-fusion.fasta \
  --peptides=input/my-peptides.csv \
  --junctions=input/my-junctions.csv \
  --output=results/my-fusion
```

Pass every input explicitly: an omitted input falls back to the bundled example file. Relative paths resolve from the project root. [Input formats](https://github.com/LangeLab/FusionPep/wiki/Input-Formats) describes each file, and `Rscript run_fusion_mapper.R --help` lists all options.

## Use from R

The same analysis is available as R functions, which return the results without writing files:

```r
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(path)

result <- run_fusion_analysis(
  sequence_file = "input/sequences.fasta",
  peptide_file = "input/peptides.csv",
  junction_file = "input/fusion_junctions.csv"
)
result$peptide_summary
```

[R interface](https://github.com/LangeLab/FusionPep/wiki/R-Interface) covers the arguments, the returned object, and writing outputs and reports.

## Documentation

The [wiki](https://github.com/LangeLab/FusionPep/wiki) holds the detailed documentation. Its source is the `wiki/` directory of this repository.

- [Getting started](https://github.com/LangeLab/FusionPep/wiki/Getting-Started): prerequisites, setup, commands, options, and path behavior.
- [Input formats](https://github.com/LangeLab/FusionPep/wiki/Input-Formats): protein records, peptide tables, and junction definitions.
- [Interpreting results](https://github.com/LangeLab/FusionPep/wiki/Interpreting-Results): mapping classes, junction decisions, coverage, and scientific limits.
- [Reports and figures](https://github.com/LangeLab/FusionPep/wiki/Reports-and-Figures): report navigation, figures, printing, and saved outputs.
- [Worked example](https://github.com/LangeLab/FusionPep/wiki/Worked-Example): the bundled PML::RARA inputs and their results.
- [Troubleshooting](https://github.com/LangeLab/FusionPep/wiki/Troubleshooting): setup, input, and interpretation checks.
- [Architecture](https://github.com/LangeLab/FusionPep/wiki/Architecture): analysis flow, source organization, and development checks.

## Development

Run the tests with the project library installed:

```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

`Rscript check_project.R` runs the full local gate: parsing, lint, tests, the bundled example, and output checks. It regenerates `results/`, so keep any earlier run elsewhere. Continuous integration runs the same gate with a coverage floor on Linux, and the tests on macOS and Windows. Changes are recorded in [NEWS.md](NEWS.md). Pushing a `vX.Y.Z` tag that matches the version publishes a GitHub release with those notes.

Issue reports and contributions are welcome through [GitHub issues](https://github.com/LangeLab/FusionPep/issues).

## Citation

If you use FusionPep in your work, please cite it:

```bibtex
@software{ergin_fusionpep_2026,
  author  = {Ergin, Enes K. and Conrrero, Agustina and Lange, Philipp F.},
  title   = {{FusionPep}: Peptide mapping and junction review for fusion proteins},
  year    = {2026},
  url     = {https://github.com/LangeLab/FusionPep},
  license = {MIT},
}
```

[CITATION.cff](CITATION.cff) carries the same metadata for tools that read it, including GitHub's _Cite this repository_ button.

## License

MIT License; see [LICENSE](LICENSE). The bundled reference sequences keep the source attribution in [input/README.md](input/README.md).
