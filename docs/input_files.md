# Input files

All inputs are read from `data/` (or from the directory given by `CMT_DATA_DIR`). Library identifiers must be
consistent across files. Columns not listed below are ignored.

## Shared

**`library_metadata.csv`**: one row per analysed library (34 libraries), in the column order used for the
Figure 2 heatmap. Used by Figures 2 and 3.

| Column | Content |
|---|---|
| `sample` | library identifier |
| `center_code` | `SHGH` (Center 1), `SY` (Center 2), `HG` (Center 3) or `MIX` (cross-center pooled library) |
| `library_type` | `Individual/same-patient combined` or `Cross-patient pooled` |

## Figure 2

**`figure2_display_matrix.csv`**: taxa shown in the heatmap (rows; grouped by kingdom in display order) by
library (columns); host-excluded, microbial-renormalized Kraken2/Bracken relative abundance as a fraction.

**`figure2_display_rowanno.csv`**: one row per displayed taxon (same row names), column `kingdom`
(`Viruses`, `Bacteria`, `Fungi` or `Other eukaryota`).

**`bracken_species_relative_abundance.csv`**: Kraken2/Bracken species-level relative abundance. Column `name`
(species) and one column per library; taxonomy columns (`taxonomy_id`, `taxonomy_lvl`, `genus_name`,
`genus_id`, `family_name`, `family_id`, `order_name`, `order_id`) may be present and are ignored.

## Figure 3

**`arg_read_level_metrics.tsv`**: read-level ARG alignment metrics (DIAMOND BLASTX against CARD), one row per
library and ARG reference.

| Column | Content |
|---|---|
| `sample` | library identifier |
| `feature_id` | CARD reference label, e.g. `gb|CAE51638.1|ARO:3002641|APH(3')-Ia` |
| `read_count` | reads assigned to the reference |
| `coverage_breadth_pct` | percentage of reference protein residues covered by at least one read |
| `depth_per_million` | DPM: mean per-residue depth (aligned reference residues / reference length), divided by total library reads, x 10^6 |

**`arg_feature_annotation.tsv`**: one row per ARG reference.

| Column | Content |
|---|---|
| `feature_id` | as above |
| `gene_name`, `ARO_ID` | CARD gene name and ARO identifier |
| `Primary_Drug_Class` | primary drug class used in the paper |
| `CARD_Drug_Class`, `CARD_Resistance_Mechanism`, `CARD_AMR_Gene_Family` | CARD annotation |

**`species_relative_abundance.tsv`**: MetaPhlAn4 relative abundance (%) of the 20 prioritized species. First
column: MetaPhlAn clade name ending in `s__Genus_species`; one column per library.

**`species_annotation.csv`**: one row per prioritized species: `species`, `evidence_tier`
(`Conjunctivitis-associated`, `Indirect ocular/infection support` or `Background detection`), `max_ra_pct`
(maximum relative abundance, %). Species with `Background detection` are not shown in Figure 3C.

**`arg_contig_host_assignments.tsv`**: high-confidence ARG-host assignments on assembled contigs (RGI Perfect or
Strict hit and host sequence on the same contig, host confirmed against NCBI nr and nt at >= 90% nucleotide
identity). Columns used: `sample`, `contig`, `best_hit_aro` (RGI best-hit name), `host_taxon`.

## Figure 4

**`hadv_tree.treefile`**: IQ-TREE maximum-likelihood tree of the study and reference HAdV genomes (Newick),
including the outgroup `OUTGROUP__AC_000008_1`.

**`hadv_tip_annotation.csv`**: one row per tree tip: `name` (tip label), `label` (library identifier for study
genomes, accession for references), `source_type` (`reference` or study), `hadv_major_type`.

**`library_id_map.csv`**: column 1 = published library identifier (e.g. `Center 2-05`), column 2 = library
identifier used in the tip table and in the HCoV-229E files (both are listed in Table S2).

**`hcov229e_tree.treefile`**: IQ-TREE tree of the study HCoV-229E genome and references, including the outgroup
HCoV-NL63 `NC_005831.2`; node labels carry SH-aLRT/ultrafast bootstrap support.

**`hcov229e_reference_table.csv`**: one row per sequence: `accession`, `lineage_or_type` (`Study sequence` for
the study genome), `country_or_region`, `year`.

**`hcov229e_spike_aa_matrix.csv`**: spike amino acids at the polymorphic sites; row names = sequence names in the
tree, columns `S<position>`.

**`hcov229e_spike_domains.csv`**: `col_index` (column of the amino-acid matrix) and `domain` (spike domain).

**`hcov229e_private_sites.csv`**: substitutions specific to the study genome: `ref_pos`, `consensus_aa`,
`study_aa`, `label` (e.g. `S:N20K`).
