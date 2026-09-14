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
