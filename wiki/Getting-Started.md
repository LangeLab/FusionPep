# Getting started

FusionPep currently runs as an R project through `run_fusion_mapper.R`. Start from a local copy containing that script, `setup_renv.R`, `DESCRIPTION`, `renv.lock`, and the bundled `input/` files. Run the commands below from the project root.

## Prerequisites

- R 4.6.0 or newer. The current lockfile records R 4.6.1 and Bioconductor 3.23.
- No separate `renv` installation: `renv/activate.R` downloads the locked `renv` version when R starts in the project directory.
- An R installation with Cairo graphics for vector PDF export.
- Network access for initial dependency setup. A system build toolchain may be needed when dependencies compile from source. On Linux, Bioconductor packages compile from source and need the zlib development headers (`zlib1g-dev` on Debian and Ubuntu, `zlib-devel` on Fedora); `pak` installs the other system libraries it can identify when it has administrator rights.

You can check the R version and Cairo graphics with:

```bash
Rscript -e 'print(getRversion()); print(capabilities("cairo"))'
```

Cairo must report `TRUE` for the current PDF export function. The analysis dependencies are `Biostrings`, `IRanges`, `pwalign`, `ggplot2`, `htmltools`, and `base64enc`.

## Set up the environment

```bash
Rscript setup_renv.R
```

The setup script activates a project-local `renv` library, establishes `pak` if needed, and installs dependencies with `pak`. It installs exactly the versions recorded in `renv.lock`, including `testthat` for the development tests, and does not rewrite the lockfile. A Bioconductor package whose locked version has since been superseded is installed by `renv` from that Bioconductor release's source archive. Package libraries and caches stay under `renv/`; the script does not install into a global R library.

Routine analysis checks the project-local packages against the lockfile. It does not install missing packages. Even `--help` currently runs after these environment checks, so complete setup first.

## Run the reference example

Choose a separate output directory:

```bash
Rscript run_fusion_mapper.R --output=results/example
```

Open `results/example/fusion_peptide_mapper_report.html` in a browser. The same directory contains the tables, R result, warnings, and a `figures/` directory with PNG and PDF exports. See the [worked example](Worked-Example) for the expected peptide decisions and coverage values.

The bundled peptide rows are reference sequences and controls. They are not PSMs measured by this project. The [worked example](Worked-Example) explains their sources and the longer peptide's unconfirmed literature attribution.

Reusing an output directory replaces the generated files with those names. Use a new directory when you want to retain an earlier run. Running the script without `--output` writes directly into `results/`.

## Run your own inputs

Prepare the three files using [Input formats](Input-Formats). The following is a command template: create the named files or replace the paths with your own.

```bash
Rscript run_fusion_mapper.R \
  --fasta=input/my-fusion.fasta \
  --peptides=input/my-peptides.csv \
  --junctions=input/my-junctions.csv \
  --output=results/my-fusion
```

Pass all applicable input paths explicitly. Omitted inputs retain their example defaults, including the junction CSV. FusionPep does not infer a replacement junction from your new FASTA.

The default FASTA identifiers are `Fusion`, `ParentA`, and `ParentB`. Put accessions and other descriptions after those identifiers in the headers. The default peptide column is `peptide`; use `--peptide-column=NAME` for another column name. Custom FASTA identifiers can also be selected through the [R interface](R-Interface).

## Command-line options

The runner accepts these options:

| Option                  | Meaning                               | Default                      |
| ----------------------- | ------------------------------------- | ---------------------------- |
| `--fasta=PATH`          | Protein FASTA                         | `input/sequences.fasta`      |
| `--peptides=PATH`       | Peptide or PSM CSV                    | `input/peptides.csv`         |
| `--junctions=PATH`      | Explicit junction CSV                 | `input/fusion_junctions.csv` |
| `--output=PATH`         | Directory for generated files         | `results/`                   |
| `--peptide-column=NAME` | Sequence column in the peptide CSV    | `peptide`                    |
| `--exact-il`            | Keep I and L distinct during matching | I/L equivalence              |
| `--help`, `-h`          | Print the available options           | N/A                          |

Use `--name=value` for options with values. Unknown or repeated options are rejected. Quote an entire argument when its path contains spaces.

Relative input and output paths resolve from the project root, even when the runner is invoked from another working directory. Use absolute paths for files elsewhere. R functions follow the R session's working directory instead.

To inspect mappings without supplied junction metadata, use `junction_file = NULL` through the R interface. The CLI does not have a flag that disables its default junction file.

Continue with [Interpreting results](Interpreting-Results) or [Troubleshooting](Troubleshooting).
