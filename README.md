# Metatranscriptomics of acute conjunctivitis: figure code

R and Python scripts that generate Figures 2, 3 and 4 of our manuscript on RNA metatranscriptomics of
conjunctival swabs from patients with acute conjunctivitis (manuscript under review). Upstream processing (read
trimming, host removal, taxonomic profiling, ARG annotation, assembly and tree inference) used published software,
as described in the STAR Methods of the paper.

## Running

Put the input tables listed below in `data/` and run from the repository root. Figures (vector PDF) and
intermediate tables are written to `output/`, one subfolder per figure.

```bash
Rscript figure2/figure2.R                   # Figure 2: community heatmap; Shannon index and PCoA by center

python figure3/01_arg_quantification.py     # ARG filtering, per-library DPM matrix, 35 high-frequency ARGs
python figure3/02_taxon_arg_correlation.py  # species-ARG Spearman correlations
python figure3/03_arg_contig_colocation.py  # species-ARG pairs on the same assembled contig
Rscript figure3/figure3.R                   # Figure 3: ARG landscape

Rscript figure4/figure4.R                   # Figure 4: HAdV and HCoV-229E phylogenies, spike diversity
```

Versions used for the paper: R 4.5.2 with ggplot2 4.0.2, patchwork 1.3.2, pheatmap 1.0.13, gtable 0.3.6,
readr 2.2.0, dplyr 1.2.1, vegan 2.7.3, ape 5.8.1, ggtree 4.0.5, treeio 1.34.0, gridGraphics 0.5.1 and png 0.1.9;
Python 3.13.12 with pandas 3.0.1, numpy 2.4.3 and scipy 1.18.0. The published figures were rendered on macOS
(`quartz` device); on other platforms the scripts use `cairo_pdf`, which can give small differences in fonts.

## Input files

| File | Figure | Content |
|---|---|---|
| `library_metadata.csv` | 2, 3 | one row per analysed library: `sample`, `center_code`, `library_type` |
| `figure2_display_matrix.csv` | 2A | displayed taxa (rows) by library: Kraken2/Bracken relative abundance |
| `figure2_display_rowanno.csv` | 2A | `kingdom` of each displayed taxon |
| `bracken_species_relative_abundance.csv` | 2B | Kraken2/Bracken species relative abundance: `name` and one column per library |
| `arg_read_level_metrics.tsv` | 3 | per library and CARD reference: `sample`, `feature_id`, `read_count`, `coverage_breadth_pct`, `depth_per_million` |
| `arg_feature_annotation.tsv` | 3 | CARD annotation per `feature_id`, including `gene_name` and `Primary_Drug_Class` |
| `species_relative_abundance.tsv` | 3C | MetaPhlAn4 relative abundance (%) of the prioritized species, one column per library |
| `species_annotation.csv` | 3C | `species`, `evidence_tier`, `max_ra_pct` |
| `arg_contig_host_assignments.tsv` | 3C | ARG-host assignments on assembled contigs: `sample`, `contig`, `best_hit_aro`, `host_taxon` |
| `hadv_tree.treefile`, `hadv_tip_annotation.csv` | 4A | HAdV IQ-TREE tree (Newick) and tip annotation |
| `library_id_map.csv` | 4 | published library identifiers (column 1) and the identifiers used in the tree files (column 2) |
| `hcov229e_tree.treefile`, `hcov229e_reference_table.csv` | 4B | HCoV-229E IQ-TREE tree and reference metadata |
| `hcov229e_spike_aa_matrix.csv`, `hcov229e_spike_domains.csv`, `hcov229e_private_sites.csv` | 4B | spike amino acids at polymorphic sites, spike domain of each site and study-specific substitutions |

## Data availability

Host-depleted sequencing reads are deposited in the Genome Sequence Archive of the National Genomics Data Center
(China National Center for Bioinformation) under BioProject PRJCA079448. Processed result tables are provided as
Supplementary Tables of the article; reference sequences used in the phylogenies are listed in Table S7.

## License and citation

MIT License (see `LICENSE`). If you use this code, please cite the associated article (reference to be added upon
publication).
