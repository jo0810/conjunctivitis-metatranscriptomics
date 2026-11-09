# Upstream processing

Steps upstream of the scripts in this repository used published software with the versions, databases and
parameters below (as described in the STAR Methods). They produce the input tables listed in
[`input_files.md`](input_files.md).

## Read processing

| Step | Software | Settings |
|---|---|---|
| Adapter and quality trimming | fastp v0.24.0 | adapter detection, quality filtering and low-complexity removal; Phred quality >= 20, average quality >= 20, <= 5 N bases, length >= 50 bp |
| Host and rRNA removal | Bowtie2 v2.5.4 | removal of reads aligning to the GRCh38.p14 primary assembly and to a SILVA release 138.2 LSURef/SSURef NR99 rRNA index |

## Taxonomic profiling (Figures 2 and 3)

| Step | Software | Settings |
|---|---|---|
| Taxonomic classification | Kraken2 v2.1.3, Bracken v3.0.1 | PlusPF database (k2_pluspf_20240904); Bracken abundance at species, genus and family level with 150-bp read length |
| Species-level profiles | MetaPhlAn4 v4.2.4 | Bowtie2 mapping to the mpa_vJun23_CHOCOPhlAnSGB_202403 marker database |

## Antimicrobial resistance genes (Figure 3)

| Step | Software | Settings |
|---|---|---|
| Read-level ARG detection | DIAMOND v2.1.10 BLASTX against CARD v4.0.1 | more-sensitive mode; e-value 1e-5; identity >= 80%; query coverage >= 50%; up to five targets per read; best hit per mate by bitscore, then e-value, then identity |
| Read-level metrics | (summary of alignments) | read count; coverage breadth (% of reference residues covered by >= 1 read); DPM (aligned reference residues / reference length / total library reads x 10^6) |
| Contig assembly | rnaSPAdes | per library, post-host non-rRNA reads |
| ARG open reading frames on contigs | RGI with CARD | contig mode; Perfect and Strict hits retained, Loose hits discarded |
| Contig host assignment | DIAMOND BLASTX against NCBI nr; BLAST+ v2.16.0+ BLASTN against NCBI nt | nr: e-value 1e-5, up to five targets, best hit by bitscore resolved to species or genus from the NCBI taxonomy identifier; assignments retained only when confirmed against nt at >= 90% nucleotide identity |

## Viral genomes and phylogenetics (Figure 4)

| Step | Software | Settings |
|---|---|---|
| Viral contig assembly | MEGAHIT v1.2.9 | post-host non-rRNA reads |
| Candidate viral contigs | DIAMOND v2.1.10 BLASTX against NCBI nr | contigs with top hits in the kingdom Viruses |
| Filtering | BLAST+ v2.16.0+ BLASTN against NCBI nt | removal of host-derived sequences, endogenous viral elements, vector sequences and background signals |

Human adenovirus: nine near-complete study genomes, 19 HAdV-D reference genomes and one HAdV-C outgroup.

```bash
mafft --adjustdirection --auto --thread 8 --inputorder hadv_29_unaligned.fasta > hadv_29_aligned.fasta
iqtree -s hadv_29_aligned.fasta -m MFP -B 1000 --alrt 1000 -T 4 --seed 20260717 --prefix hadv_fullgenome_ml
```

MAFFT v7.526; IQ-TREE v3.0.1 with ModelFinder (GTR+F+I+R3 selected by BIC), 1,000 ultrafast bootstrap and
1,000 SH-aLRT replicates.

HCoV-229E: maximum-likelihood tree with IQ-TREE under GTR+I+G with 1,000 ultrafast bootstrap and 1,000 SH-aLRT
replicates, rooted with HCoV-NL63 (NC_005831.2). Spike amino-acid diversity was summarized across 26 sequences at
100 polymorphic sites; the two substitutions specific to the study genome were checked against 127
quality-filtered full-genome HCoV-229E references.
