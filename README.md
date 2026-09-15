# BSP457 — Porcine Airway ALI Culture PRCV Challenge

Sequencing and analysis pipeline for a porcine bronchial/tracheal air-liquid
interface (ALI) culture challenge experiment: two parallel analyses on the
same 48-library design — (1) challenge-virus genome sequencing and variant
calling against the reference genome, and (2) host (pig) transcriptome
differential expression.

## Experimental design

6 animals (4100, 4101, 4103, 4108, 4194, 4198), each contributing bronchial
and tracheal ALI cultures, mock-treated or infected, sampled at 6 and 24
hours post-infection (hpi). Full 2 (site) x 2 (time) x 2 (infection)
factorial, 48 libraries total (one per animal per condition).

## Repository structure

```
GITHUB/
├── scripts/
│   ├── kraken_illumina.sh              # kraken2 classification
│   ├── kraken_illumina_extract.sh      # extract virus-classified reads
│   ├── guide.sh                        # per-sample mapping + variant calling
│   ├── filter.sh                       # variant filtering + caller concordance
│   ├── coverage_viz.R                  # viral genome coverage figures
│   ├── Snakefile.kallisto              # host transcript quantification (not reviewed here)
│   ├── Snakemake.qc                    # read QC/trimming (not reviewed here)
│   ├── 1.R                             # likely: host DE analysis (not confirmed)
│   └── abundanceplot.R                 # likely: PCA/TPM figures (not confirmed)
├── DATA/
│   ├── DQ811787.fasta / .fasta.fai / .gb   # PRCV reference genome
│   └── adapters.fa
├── TABLES/                             # DE result tables, sleuth tables, TPMs
└── FIGURES/                            # PCA, volcano, venn, coverage plots
```

Four scripts not reviewed for this README (`Snakefile.kallisto`,
`Snakemake.qc`, `1.R`, `abundanceplot.R`) are listed for completeness only;
their descriptions above are inferred from filenames and output files, not
from reading their contents.

## Reference genome

`DATA/DQ811787.fasta` — Porcine respiratory coronavirus (PRCV), GenBank
accession DQ811787. `.fai` is the samtools faidx index; `.gb` is the
GenBank flat-file annotation. `DATA/adapters.fa` holds adapter sequences
for read trimming.

## Pipeline

### 1. Read QC / trimming
`scripts/Snakemake.qc` — not reviewed for this README.

### 2. Taxonomic classification and viral read extraction

`scripts/kraken_illumina.sh`
- Runs kraken2 (`--db /home/share/`, 3 threads) on every sample in `FASTQS/`.
- Output: `KRAKEN_FILES/<sample>` (classification), `KRAKEN_REPORT/<sample>` (report).

`scripts/kraken_illumina_extract.sh`
- Pulls reads classified to taxon 11146 (the challenge-virus lineage; run
  with `--include-children --include-parents`) out of each sample's kraken2
  output, via KrakenTools' `extract_kraken_reads.py`.
- Output: `KRAKEN_VALID_FILES/<sample>.fastq`.
- Note: a second taxon (`9822`, commented `# pig`) is assigned but never
  used in an extraction call — looks unfinished, presumably intended for a
  parallel host-read extraction.

### 3. Viral genome mapping, variant calling and consensus

`scripts/guide.sh <sample>` (one sample per invocation; the header comment
shows this run in parallel over `KRAKEN_VALID_FILES/*.fastq`)

- Maps reads to `DQ811787.fasta` with `minimap2 -ax sr`.
- `samtools`: sort by queryname → `fixmate -m` → sort by coordinate →
  `markdup -r -s` → `BAM_FINAL/<sample>.merged.bam` (deduplicated,
  coordinate-sorted, indexed).
- QC: `samtools flagstat`, `samtools coverage` (→ `COVERAGE/`), `samtools
  depth` (→ `DEPTH/`).
- Variant calling, two callers, both haploid (`--ploidy 1` / `-p 1`):
  - freebayes (`-F 0.01 --min-alternate-count 1 --min-alternate-fraction
    0.001`) → `FB_VCF_FILES/`
  - `bcftools mpileup | bcftools call -cv --ploidy 1` → `BCF_VCF_FILES/`
  - both normalised with `bcftools norm -m-both`
- Consensus sequences from the BCFtools calls, two haplotype-selection
  modes (`-H LA`, `-H LR`) → `CONSENSUS/`.

### 4. Variant filtering and caller concordance

`scripts/filter.sh`

- FreeBayes calls: SNPs only, `QUAL>=20`, `INFO/DP>=5`, `INFO/AO>=5`,
  `INFO/AF>0.05`.
- BCFtools calls: same depth/quality/alt-read thresholds, computed manually
  from the `DP`/`DP4` INFO tags (`awk`), biallelic SNPs only.
- Both callers additionally restricted to position `500` to
  `genome_length-500` (trims the first/last 500 bp of the reference).
- Output: `FILTERED_VCF/<sample>.{freebayes,bcftools}.filtered.vcf.gz`
  (bgzipped, tabix-indexed).
- `SNP_TABLES/<sample>.merged_snps.tsv` — union of both callers' SNPs
  (CHROM/POS/REF/ALT/caller); per-sample FreeBayes/BCFtools/union counts
  are printed to stdout.

### 5. Viral genome coverage visualisation

`scripts/coverage_viz.R`

- Reads every `COVERAGE/*.coverage.txt` (`samtools coverage` output,
  %-covered column), parses tissue/timepoint/infection/animal from the
  filename ("traheal" is matched as an alternate spelling of tracheal).
- `FIGURES/coverage.png` — % genome covered per sample, by the 8
  site x time x infection groups.
- `FIGURES/coverage_6_vs_24hpi_paired.png` — per-animal 6 -> 24 hpi
  trajectories (infected samples only, complete pairs only), faceted by
  tissue.
- `FIGURES/coverage_change_24_minus_6hpi.png` — within-animal change in
  %coverage (24 hpi minus 6 hpi), by tissue.

### 6. Host transcript quantification

`scripts/Snakefile.kallisto` — not reviewed for this README. Output
consistent with kallisto pseudo-alignment of each sample against the
Ensembl *Sus scrofa* Sscrofa11.1 (release 115) transcriptome, one output
directory per sample under `results/kallisto/<site>/<infection>/<time>hpi/
<animal>/`.

### 7. Sample sheet

`SLEUTH/table.csv` — one row per library: `sample`, `site`, `time`,
`infection`, `path` (to that library's kallisto output directory). 48 rows
(6 animals x 2 sites x 2 timepoints x 2 infection states).

### 8. Host differential expression

Transcript-level (limma/voom, blocked on animal) and gene-level (sleuth
likelihood-ratio tests) differential expression across the site x time x
infection factorial design, plus STRING protein-protein interaction and
KEGG gene set enrichment follow-up on the resulting gene lists. Script not
confirmed (likely `scripts/1.R`).

- Design: `~0 + group + animal` (group = site:time:infection, 8 levels;
  animal blocks the 6 pigs, each contributing one sample per group — a
  complete randomised block design).
- limma contrasts: infected-vs-mock within each site x time cell (4);
  infection, site and time main effects (3); infection x site and
  infection x time interactions (2).
- Significance: Benjamini-Hochberg FDR < 0.05 and |log2FC| > 1.39
  (~2.6-fold).
- sleuth LRTs: main effects of site/time/infection, plus the 2-way and
  3-way interaction terms.
- STRING PPI (v12, *Sus scrofa*, taxon 9823, score >= 700) and KEGG GSEA
  (clusterProfiler `gseKEGG`) run per contrast on each contrast's DE gene
  list.

## Outputs

### TABLES/

- `limma.res_<contrast>.csv` — full transcript-level limma/voom results,
  one file per contrast (9 contrasts).
- `<contrast>.all.csv` / `<contrast>.DE.csv` — same results annotated with
  gene IDs/symbols, and the FDR+logFC-significant subset.
- `DE_counts_summary.csv` — DE transcript counts per contrast.
- `all_DE_union_primary_contrasts.csv` — union of DE transcripts across the
  4 site x time infected-vs-mock contrasts.
- `all_DE_multiple_transcripts.csv` — genes with more than one DE
  transcript in that union.
- `sleuth_DE_<term>.csv` — gene-level LRT results (main effects and
  interaction terms).
- `tpm_values.csv` — transcript TPM matrix used for PCA.

### FIGURES/

- `PC1.PC2.*`, `PC3.PC4.*`, `plot_sample_heatmap.*` — sample-level QC.
- `<contrast>.volcano.pdf` — one volcano plot per limma contrast.
- `DE_venn_site_time.*` — overlap of DE genes across the 4 site x time
  contrasts.
- `DE_gene_expression_panels.*` — per-gene expression across all 8 groups
  for the union of DE transcripts.
- `pairs_logFC_primary_contrasts.pdf` — pairwise correlation of logFC
  across the 4 site x time contrasts.
- `plotSA.pdf` — voom/limma mean-variance trend diagnostic.
- `coverage*.png` — viral genome coverage figures (see section 5).

## Dependencies

- kraken2, KrakenTools (`extract_kraken_reads.py`)
- minimap2, samtools, freebayes, bcftools, htslib (bgzip/tabix)
- kallisto, sleuth
- R: tidyverse, limma, edgeR, STRINGdb, igraph, ggraph, clusterProfiler,
  org.Ss.eg.db, rtracklayer
- R and package versions not pinned in this README — record `sessionInfo()`
  alongside results if reproducibility matters.

## Known issues

- `kraken_illumina_extract.sh` sets an unused `taxon=9822` variable (see
  section 2).
- STRING/KEGG GSEA outputs (from step 8) write to `STRING/` and `KEGG/`,
  which aren't shown in the `TABLES/`/`FIGURES/` listing this README was
  based on — confirm those steps have been run if those directories are
  missing.
- Scripts not reviewed for this README (`Snakefile.kallisto`,
  `Snakemake.qc`, `1.R`, `abundanceplot.R`): descriptions above are
  inferred from filenames and outputs, not from reading their contents —
  check them directly before relying on this README for those steps.
  

# Methods

## Sample collection

RNA-seq read libraries were generated from bronchial and tracheal air–liquid interface (ALI) epithelial cultures derived from six pigs (IDs 4100, 4101, 4103, 4108, 4194, 4198). For each animal, cultures from both airway sites were mock-treated or infected with PRCV at MOI XX, and sampled at 6 and 24 hpi. ...

## Transcript quantification and differential expression testing

Reads were pseudo-aligned and transcript abundance quantified with kallisto using 100 bootstraps (Bray et al 2016) and the pig Sscrofa11.1 reference assembly (Warr et al 2020) and the Ensembl release-115 annotation (Martin et al 2023). A transcript-to-gene lookup table was built from the Ensembl data using R package rtracklayer v1.58.0 (REF). Transcript-level counts were imported with tximport v1.26.1 (Soneson et al 2015) and edgeR v3.40.2 (Robinson et al 2010). The animal ID was a fixed blocking factor in this randomised block design experiment. After filtering the 57,764 transcripts, 32,646 remained. The read libraries were TMM-normalised and modelled with voom (Law et al 2014). Linear models were fitted in limma v3.54.2 (Ritchie et al 2015) with empirical Bayes moderation. Transcripts were tested for differential expression (DE) across the main effects of infection, time and anatomical site using 48 samples, as well as for infection within each time-site subset of 12 samples. Transcripts with DE had a Benjamini–Hochberg FDR < 0.05 and a log2 fold change (log2FC) < |1|. This corresponded to a power to detect transcripts with true DE of 44.4% for 6 samples per group, 72.9% for 12 samples per group, 95.7% for 24 samples per group given the expression variation in this dataset. Genes with DE were tested using their aggregated to Ensembl gene IDs using Sleuth v0.30.1 (Pimentel et al 2017). This tested a full model (modelling infection, time and site) against a reduced model iteratively removing one of the three terms in turn. The full and reduced models were compared by likelihood-ratio test (LRT) to detect genes with DE that had q <0.05. Transcripts and genes with DE were modelling with dplyr v1.2.1 (REF), ggplot2 v4.0.3 (REF) and tidyr v1.3.2 (REF).

## Protein–protein interaction network analysis

To retrieve protein-protein interaction (PPI) data, senes with DE based on the infection comparison were queried against STRING v12 (Szklarczyk et al 2023) using R packages STRINGdb v2.10.1 (REF), ggraph v2.2.2 and igraph v2.3.3 for Sus scrofa using a combined-score threshold of 400. PPI networks were restricted to their largest connected component. PPI hub subnetworks were defined as nodes with degree ≥ 10.

## Pathway enrichment

For each of the main limma DE comparisons above, transcript-level moderated t-statistics were averaged per gene, and mapped to their Entrez IDs using org.Ss.eg.db v3.16.0. This was used to rank genes for gene set enrichment analysis (GSEA) against KEGG pathways (Kanehisa & Goto 2000) with clusterProfiler v4.6.2 (Wu et al 2021), gene set sizes of 10 to 500, and a p < 0.05.

---

# Results

Sleuth found 10,517 genes with DE based on the anatomical site, and 11,001 associated with time (6 vs 24 hpi). No genes with DE were found based on infection status. Nine genes with DE were found for the site-time interaction effects, and none for the site-infection and time-infection interactions.

Across the transcripts with DE associated with anatomical site, 4,417 were found. For time, 3,603 transcripts with DE were found. For infection, 38 were found, corresponding to 35 unique genes. Within infection at the bronchial site there 2 DE transcripts at 6 hpi and 44 at 24 hpi. Within infection at the trachael site, there were none at 6 hpi and three at 24 hpi. There was one DE transcript associated with the infection-time interaction, and none for the infection-site and time-site interactions. 






# References

Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). Near-optimal probabilistic RNA-seq quantification. *Nature Biotechnology*, 34(5), 525–527. https://doi.org/10.1038/nbt.3519

Kanehisa, M., & Goto, S. (2000). KEGG: Kyoto Encyclopedia of Genes and Genomes. *Nucleic Acids Research*, 28(1), 27–30. https://doi.org/10.1093/nar/28.1.27

Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014). voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biology*, 15(2), R29. https://doi.org/10.1186/gb-2014-15-2-r29

Martin, F. J., Amode, M. R., Aneja, A., et al. (2023). Ensembl 2023. *Nucleic Acids Research*, 51(D1), D933–D941.

Pimentel, H., Bray, N. L., Puente, S., Melsted, P., & Pachter, L. (2017). Differential analysis of RNA-seq incorporating quantification uncertainty. *Nature Methods*, 14(7), 687–690. https://doi.org/10.1038/nmeth.4324

R Core Team. R: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria. [version/year to be added]

Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*, 43(7), e47. https://doi.org/10.1093/nar/gkv007

Robinson, M. D., McCarthy, D. J., & Smyth, G. K. (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics*, 26(1), 139–140. https://doi.org/10.1093/bioinformatics/btp616

Soneson, C., Love, M. I., & Robinson, M. D. (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Research*, 4, 1521. https://doi.org/10.12688/f1000research.7563.1

Szklarczyk, D., Kirsch, R., Koutrouli, M., Nastou, K., Mehryary, F., Hachilif, R., Gable, A. L., Fang, T., Doncheva, N. T., Pyysalo, S., Bork, P., Jensen, L. J., & von Mering, C. (2023). The STRING database in 2023: protein–protein association networks and functional enrichment analyses for any sequenced genome of interest. *Nucleic Acids Research*, 51(D1), D638–D646. https://doi.org/10.1093/nar/gkac1000

Warr, A., Affara, N., Aken, B., et al. (2020). An improved pig reference genome sequence to enable pig genetics and genomics research. *GigaScience*, 9(6), giaa051. https://doi.org/10.1093/gigascience/giaa051

Wu, T., Hu, E., Xu, S., Chen, M., Guo, P., Dai, Z., Feng, T., Zhou, L., Tang, W., Zhan, L., Fu, X., Liu, S., Bo, X., & Yu, G. (2021). clusterProfiler 4.0: A universal enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141. https://doi.org/10.1016/j.xinn.2021.100141
