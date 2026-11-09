# Metatranscriptomics of acute conjunctivitis: figure code

This repository contains the code that generates Figures 2, 3 and 4 of our manuscript on RNA
metatranscriptomics of conjunctival swabs from patients with acute conjunctivitis (pathogen spectrum,
resistance landscape and viral genomes; manuscript under review). Each script reads processed tables and
writes the corresponding figure as a vector PDF.

## Contents

| Folder | Script | Output |
|---|---|---|
| `figure2` | `figure2.R` | Figure 2: whole-community heatmap (A); Shannon index and Bray-Curtis PCoA by center (B) |
| `figure3` | `01_arg_quantification.py` | ARG filtering (coverage breadth >= 50%; DPM > 1 or read count > 10), per-library DPM matrix and the 35 high-frequency ARGs |
| | `02_taxon_arg_correlation.py` | Spearman correlations between prioritized species and ARGs |
| | `03_arg_contig_colocation.py` | Species-ARG pairs found on the same assembled contig |
| | `figure3.R` | Figure 3: ARGs per center (A), prevalence and abundance of high-frequency ARGs (B), species-ARG correlation (C) |
| `figure4` | `figure4.R` | Figure 4: HAdV whole-genome phylogeny (A); HCoV-229E phylogeny and spike amino-acid diversity (B) |

Read processing, taxonomic and functional profiling, ARG annotation, viral identification and tree inference
used published software; versions, databases and parameters are listed in
[`docs/upstream_processing.md`](docs/upstream_processing.md) and in the STAR Methods of the paper.

## Requirements

- R 4.5.2 with ggplot2, patchwork, pheatmap, gtable, readr, dplyr, vegan, ape, ggtree, treeio, gridGraphics and png
- Python 3.10 or later with pandas, numpy and scipy

Versions used for the paper are listed in [`environment/`](environment/). The published figures were rendered on
macOS with the `quartz` device (Arial embedded); on other platforms the scripts use `cairo_pdf`, which can give
small differences in font rendering.

## Input data

Place the input tables in `data/` or point `CMT_DATA_DIR` to another directory. File names and columns are
described in [`docs/input_files.md`](docs/input_files.md).

## Running

Run from the repository root:

```bash
Rscript figure2/figure2.R

python figure3/01_arg_quantification.py
python figure3/02_taxon_arg_correlation.py
python figure3/03_arg_contig_colocation.py
Rscript figure3/figure3.R

Rscript figure4/figure4.R
```

Figures and intermediate tables are written to `output/` (or to `CMT_OUTPUT_DIR`), one subfolder per figure.

## Data availability

Sequencing reads are deposited at [repository and accession to be added]. Processed result tables are provided
as Supplementary Tables of the article; reference genome accessions are listed in Table S7.

## License

MIT License; see [`LICENSE`](LICENSE).

## Citation

If you use this code, please cite the associated article (reference to be added upon publication).
