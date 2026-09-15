# =============================================================================
# Porcine airway explant RNA-seq: site x time x infection factorial design
# Continues from your existing Part 1 (setup/grouping/data import) and
# Part 2 (PCA plots) -- this picks up from Part 3 onward.
#
# Requires from Part 1: `s2c`, a data.frame with columns
#   sample, site (bronchial/tracheal), time (6/24), infection (infected/mock),
#   path
# (matches your table.csv). If s2c isn't already in your session, it's
# re-read below.
#
# *** DATA CHECK ***
# Your message said time = 6 vs 12 hpi, but table.csv actually has time =
# 6 and 24. The code below uses 6/24 to match the file -- if 12 was correct
# and the table has a typo, fix the sheet (or the `full$time` line below)
# before running.
#
# KEY DESIGN DECISIONS (change these if you want something different):
#
#  1. Model: 8-group factorial (site x time x infection) + pig as a FIXED
#     BLOCKING FACTOR. Every pig (4100, 4101, 4103, 4108, 4194, 4198)
#     contributes a sample to all 8 site/time/infection combinations, i.e.
#     this is a complete randomised block design. Blocking on pig removes
#     animal-to-animal baseline variation from every contrast below, which
#     should give real power gains over an unblocked model. sleuth (Part 5)
#     has no blocking/random-effect mechanism, so the limma results here are
#     the better-powered, primary analysis; sleuth is a complementary,
#     bootstrap-uncertainty-aware check that ignores the pig structure.
#
#  2. Primary contrasts: infected vs mock, done separately within each
#     site x time cell (4 contrasts) -- the direct analogue of your old
#     script's "D0 vs Dx" contrasts, just swapping "day" for "site x time".
#     Main effects of infection/site/time and two infection interaction
#     contrasts (does the infection effect differ by site / by time) are
#     also built in.
#
#  3. logFC sign: every contrast is written as infected - mock, so a
#     positive logFC already means "up with infection". The sign-flip hack
#     in the old make_volcano() (needed because the D0-Dx contrasts ran
#     backwards) has been removed.
#
#  4. Natural-spline-over-time model: removed. You now have 2 timepoints,
#     so a df=3 spline isn't fittable or meaningful.
#
#  5. PRRSV viral-load trajectory section: removed. It was keyed to a
#     different dataset's sample sheet (hardcoded "P22-xxxx" IDs) that
#     doesn't correspond to these explant samples. If you have matched
#     viral quantification for this experiment, that section can be rebuilt
#     against the new sample names.
#
#  6. KEGG GSEA / STRING PPI: previously combined every day-comparison into
#     one undirected sqrt(sum(t^2)) ranking. That doesn't make sense here --
#     "infected vs mock" has a clear direction, and mixing it with site/time
#     contrasts wouldn't be meaningful -- so both are now run separately per
#     contrast using the signed t-statistic.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tximport)
library(edgeR)
library(limma)
library(ggVennDiagram)
library(rtracklayer)
library(sleuth)
library(STRINGdb)
library(igraph)
library(ggraph)
library(clusterProfiler)
library(org.Ss.eg.db)

sessioninfo::package_info("dplyr")
sessioninfo::package_info("tidyr")
sessioninfo::package_info("ggplot2")
sessioninfo::package_info("tximport")
sessioninfo::package_info("edgeR")
sessioninfo::package_info("limma")
sessioninfo::package_info("ggVennDiagram")
sessioninfo::package_info("rtracklayer")
sessioninfo::package_info("sleuth")
sessioninfo::package_info("STRINGdb")
sessioninfo::package_info("igraph")
sessioninfo::package_info("ggraph")
sessioninfo::package_info("clusterProfiler")
sessioninfo::package_info("org.Ss.eg.db")

s2c <- read.csv("../SLEUTH/table.csv", header = TRUE, sep = "\t")
s2c$path
####### 3 - get biomart data ##################################################
# Unchanged from before -- builds the transcript-to-gene lookup table used
# throughout. (Kept as rtracklayer/GTF-based rather than biomaRt, since
# biomaRt was unreliable for this genome build.)

gtf <- import("Sscrofa11.1.115.gtf.gz")
gtf_df <- as.data.frame(gtf)
t2g <- gtf_df %>%
  filter(type == "transcript") %>%
  dplyr::select(    target_id          = transcript_id,
    transcript_version,
    ens_gene           = gene_id,
    ext_gene           = gene_name,
    target_biotype     = transcript_biotype  )
t2g$target_id <- paste(t2g$target_id, t2g$transcript_version, sep = ".")
t2g$transcript_version <- NULL
str(t2g)

####### 4 - limma exploration, site x time x infection ########################

## 4a. metadata & design ------------------------------------------------------

full <- s2c
full$site      <- factor(full$site)
full$infection <- factor(full$infection, levels = c("mock", "infected")) 
# mock = reference
full$time      <- factor(full$time)   # 6 / 24 hpi -- see DATA CHECK note above
full$animal    <- factor(sub("^[A-Za-z]+_([0-9]+)_.*$", "\\1", full$sample))
full$group     <- factor(paste(full$site, full$time, full$infection, sep = "_"))

table(full$group)             # sanity check -- should be 6 samples per group
table(full$animal, full$group)  # sanity check -- every pig x every group, once each

design <- model.matrix(~0 + group + animal, data = full)
colnames(design) <- make.names(colnames(design))
str(design)

## 4b. import counts, filter, normalise, fit ----------------------------------

files <- paste(s2c$path, "abundance.h5", sep = "/")
txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)

y <- DGEList(txi.kallisto$counts)
colnames(y) <- full$sample
dim(y) # 57764

keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
dim(y)   # 32646 transcripts retained after filtering

y <- calcNormFactors(y)
v <- voom(y, design)

fit <- lmFit(v, design)
fit <- eBayes(fit, trend = TRUE)
pdf("plotSA.pdf")
plotSA(fit)   # check the mean-variance trend -- keep trend=TRUE in eBayes() if there is one
dev.off()

## 4c. contrasts ---------------------------------------------------------------
# Column names below follow model.matrix()'s default naming: "group" + level,
# e.g. groupbronchial_6_infected. Edit/add contrasts here as needed.

cont <- makeContrasts(
  Bronchial_6hpi_Inf_vs_Mock  = groupbronchial_6_infected  - groupbronchial_6_mock,
  Bronchial_24hpi_Inf_vs_Mock = groupbronchial_24_infected - groupbronchial_24_mock,
  Tracheal_6hpi_Inf_vs_Mock   = grouptracheal_6_infected   - grouptracheal_6_mock,
  Tracheal_24hpi_Inf_vs_Mock  = grouptracheal_24_infected  - grouptracheal_24_mock,

  Infection_overall = (groupbronchial_6_infected + groupbronchial_24_infected +
                        grouptracheal_6_infected + grouptracheal_24_infected) / 4 -
                       (groupbronchial_6_mock + groupbronchial_24_mock +
                        grouptracheal_6_mock + grouptracheal_24_mock) / 4,

  Site_Bronchial_vs_Tracheal = (groupbronchial_6_infected + groupbronchial_6_mock +
                                 groupbronchial_24_infected + groupbronchial_24_mock) / 4 -
                                (grouptracheal_6_infected + grouptracheal_6_mock +
                                 grouptracheal_24_infected + grouptracheal_24_mock) / 4,

  Time_24_vs_6hpi = (groupbronchial_24_infected + groupbronchial_24_mock +
                      grouptracheal_24_infected + grouptracheal_24_mock) / 4 -
                     (groupbronchial_6_infected + groupbronchial_6_mock +
                      grouptracheal_6_infected + grouptracheal_6_mock) / 4,

  Infection_x_Site = ((groupbronchial_6_infected - groupbronchial_6_mock) +
                       (groupbronchial_24_infected - groupbronchial_24_mock)) / 2 -
                      ((grouptracheal_6_infected - grouptracheal_6_mock) +
                       (grouptracheal_24_infected - grouptracheal_24_mock)) / 2,

  Infection_x_Time = ((groupbronchial_24_infected - groupbronchial_24_mock) +
                       (grouptracheal_24_infected - grouptracheal_24_mock)) / 2 -
                      ((groupbronchial_6_infected - groupbronchial_6_mock) +
                       (grouptracheal_6_infected - grouptracheal_6_mock)) / 2,

  levels = design
)

primary_contrasts <- c("Bronchial_6hpi_Inf_vs_Mock", "Bronchial_24hpi_Inf_vs_Mock",
                        "Tracheal_6hpi_Inf_vs_Mock", "Tracheal_24hpi_Inf_vs_Mock")

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE)

## 4d. per-contrast results -----------------------------------------------------

res_list <- list()
for (cn in colnames(cont)) {
  res_list[[cn]] <- topTable(fit2, coef = cn, number = Inf)
  write.csv(res_list[[cn]], paste0("limma.res_", cn, ".csv")) }

de_counts <- sapply(res_list, function(r) sum(r$adj.P.Val <= 0.05 & abs(r$logFC) > 1, na.rm = TRUE))
print(de_counts)
write.csv(data.frame(contrast = names(de_counts), n_DE_transcripts = de_counts),
          "DE_counts_summary.csv", row.names = FALSE)

## 4e. volcano plots -------------------------------------------------------------

make_volcano <- function(results, comparison, t2g, fdr_cutoff =0.05, logfc_cutoff=1){
  transcript_ids <- rownames(results)
  if (grepl("\\.", transcript_ids[1])) {
    results$target_id <- rownames(results)
    results_with_genes <- left_join(
      data.frame(target_id = rownames(results), results, row.names = NULL),
      t2g %>% dplyr::select(target_id, ens_gene, ext_gene),
      by = "target_id")
  } else {
    results$ensembl_transcript_id <- rownames(results)
    results_with_genes <- left_join(
      data.frame(ensembl_transcript_id = rownames(results), results, row.names = NULL),
      t2g %>% dplyr::select(ensembl_transcript_id = target_id, ens_gene, ext_gene),
      by = "ensembl_transcript_id") }
  results_with_genes$gene_label <- ifelse(
    !is.na(results_with_genes$ext_gene) & results_with_genes$ext_gene != "",
    results_with_genes$ext_gene,
    ifelse(!is.na(results_with_genes$ens_gene) & results_with_genes$ens_gene != "",
           results_with_genes$ens_gene,
           ifelse(grepl("\\.", transcript_ids[1]),
                  results_with_genes$target_id,
                  results_with_genes$ensembl_transcript_id)) )

  # logFC is already "infected - mock" from the contrast itself -- positive
  # already means up with infection, so (unlike the old script) no sign flip.
  results_with_genes$significant <- ifelse(
    results_with_genes$adj.P.Val < fdr_cutoff & abs(results_with_genes$logFC) > logfc_cutoff,
    "Significant", "Not")

  top_genes <- results_with_genes %>% filter(significant == "Significant")
  cat(comparison, ":", nrow(top_genes), "DE transcripts\n")

  write.csv(results_with_genes, paste0( comparison, ".all.csv"), row.names = FALSE)
  write.csv(top_genes, paste0( comparison, ".DE.csv"), row.names = FALSE)

  top_genes2 <- top_genes %>% arrange(adj.P.Val) %>% slice_head(n = 10)

  pdf(paste0( comparison, ".volcano.pdf"), width = 6, height = 4)
  p <- ggplot(results_with_genes, aes(x = logFC, y = -log10(adj.P.Val), color = significant)) +
    geom_point(shape = 1, size = 1.2, alpha = 0.2) +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "darkgray") +
    scale_color_manual(values = c(Significant = "red", Not = "black")) +
    geom_text_repel(data = top_genes2, aes(label = gene_label), box.padding = 0.4,
                     min.segment.length = 0, max.overlaps = 220, force = 10, size = 6) +
    labs(x = "log2(FC)", y = "-log10(FDR)", title = comparison) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 16, face = "bold"),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12))
  print(p)
  dev.off()
  return(results_with_genes)
}

vol_list <- list()
for (cn in names(res_list)) {
  vol_list[[cn]] <- make_volcano(res_list[[cn]], cn, t2g) }

## 4f. how correlated is the infection response across site/time? --------------
# (replaces the old pairs() plot -- this uses the actual per-contrast logFCs,
# not an arbitrary slice of a single topTable)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y))
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if (missing(cex.cor)) cex.cor <- 1 + 0.4 / strwidth(txt)
  text(0.5, 0.5, txt, cex = 1 + cex.cor * Cor) }

pdf("pairs_logFC_primary_contrasts.pdf", width = 6, height = 6)
pairs(fit2$coefficients[, primary_contrasts], cex = 0.1,
      upper.panel = panel.cor, lower.panel = panel.smooth, cex.labels = 1.2)
dev.off()

## 4g. overlap of DE genes across the 4 site x time comparisons ----------------

de_gene_sets <- lapply(primary_contrasts, function(cn) {
  df <- vol_list[[cn]]
  unique(na.omit(df$ens_gene[df$significant == "Significant"])) })
names(de_gene_sets) <- primary_contrasts

p_venn <- ggVennDiagram(de_gene_sets, label = "count") +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(title = "DE genes (infected vs mock), by site x time")
ggsave("DE_venn_site_time.pdf", p_venn, width = 7, height = 7)
ggsave("DE_venn_site_time.png", p_venn, width = 7, height = 7, dpi = 300)

## 4h. union of DE transcripts across the 4 primary contrasts ------------------

all_de <- bind_rows(lapply(primary_contrasts, function(cn) {
  # mutate() (not $<-) so this doesn't error when a contrast has 0 DE rows
  # (base R's `df$contrast <- cn` throws "replacement has 1 row, data has 0"
  # on a 0-row data.frame -- this is what halted the run at Tracheal_6hpi)
  vol_list[[cn]] %>%
    filter(significant == "Significant") %>%
    mutate(contrast = cn) }))
write.csv(all_de, "all_DE_union_primary_contrasts.csv", row.names = FALSE)
cat("Unique DE transcripts (>=1 primary contrast):", length(unique(all_de$target_id)), "\n")
cat("Unique DE genes (>=1 primary contrast):", length(unique(na.omit(all_de$ens_gene))), "\n")

gene_counts <- all_de %>%
  distinct(target_id, ens_gene) %>%
  count(ens_gene, name = "n_transcripts")
multi_transcript_genes <- gene_counts %>% filter(n_transcripts > 1)
multi_transcript_de <- all_de %>%
  semi_join(multi_transcript_genes, by = "ens_gene") %>%
  arrange(ens_gene)
write.csv(multi_transcript_de, "all_DE_multiple_transcripts.csv", row.names = FALSE)
cat("Genes with >1 DE transcript:", nrow(multi_transcript_genes), "\n")

## 4i. expression of DE transcripts across all 8 groups -------------------------

colnames(v$E) <- full$sample
group_levels <- c("bronchial_6_mock", "bronchial_6_infected",
                   "bronchial_24_mock", "bronchial_24_infected",
                   "tracheal_6_mock", "tracheal_6_infected",
                   "tracheal_24_mock", "tracheal_24_infected")

gene_lookup <- t2g %>% dplyr::select(target_id, ens_gene, ext_gene) %>% distinct()
gene_lookup$symbol <- ifelse(!is.na(gene_lookup$ext_gene) & gene_lookup$ext_gene != "",
                              gene_lookup$ext_gene, gene_lookup$ens_gene)

de_transcripts <- unique(all_de$target_id)
expr_mat <- as.data.frame(v$E)
expr_mat$target_id <- rownames(expr_mat)
expr_mat <- expr_mat %>%
  filter(target_id %in% de_transcripts) %>%
  left_join(gene_lookup %>% dplyr::select(target_id, symbol), by = "target_id")
expr_mat$plot_name <- ifelse(is.na(expr_mat$symbol) | expr_mat$symbol == "",
                              expr_mat$target_id, paste0(expr_mat$symbol, " | ", expr_mat$target_id))

expr_long <- expr_mat %>%
  pivot_longer(cols = -c(target_id, symbol, plot_name), names_to = "sample", values_to = "expression") %>%
  left_join(full[, c("sample", "site", "time", "infection", "group")], by = "sample")
expr_long$group <- factor(expr_long$group, levels = group_levels)

# cluster genes by their across-group mean profile so facets aren't alphabetical
expr_cluster <- expr_long %>%
  group_by(target_id, group) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_expr)
cluster_ids <- expr_cluster$target_id
cluster_matrix <- as.matrix(expr_cluster[, -1])
rownames(cluster_matrix) <- cluster_ids
gene_cor <- cor(t(cluster_matrix), use = "pairwise.complete.obs")
hc <- hclust(as.dist(1 - gene_cor), method = "average")
gene_order <- cluster_ids[hc$order]

plot_lookup <- expr_long %>% dplyr::select(target_id, plot_name) %>% distinct()
plot_order <- plot_lookup$plot_name[match(gene_order, plot_lookup$target_id)]
expr_long$plot_name <- factor(expr_long$plot_name, levels = plot_order)

lab_df <- expr_long %>% dplyr::select(plot_name, symbol) %>% distinct()
labeller_gene <- function(x) lab_df$symbol[match(x, lab_df$plot_name)]

p <- ggplot(expr_long, aes(x = group, y = expression, colour = infection)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 1.5, alpha = 0.6) +
  stat_summary(aes(group = infection), fun = mean, geom = "line", linewidth = 0.6, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 1.3) +
  scale_colour_manual(values = c(mock = "grey50", infected = "firebrick")) +
  facet_wrap(~plot_name, ncol = 10, scales = "free_y",
             labeller = labeller(plot_name = labeller_gene)) +
  labs(x = NULL, y = "Normalised expression (voom logCPM)") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 5, face = "bold"),
        axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        panel.spacing = unit(0.02, "lines"))

n_genes <- length(unique(expr_long$plot_name))
n_rows <- ceiling(n_genes / 10)
ggsave("DE_gene_expression_panels.pdf", p, width = 22, height = max(6, n_rows * 2),
       dpi = 300, limitsize = FALSE)
ggsave("DE_gene_expression_panels.png", p, width = 22, height = max(6, n_rows * 2),
       dpi = 300, limitsize = FALSE, bg = "white")

####### 5 - Sleuth across genes (gene-level LRTs) ##############################
# NOTE: sleuth has no blocking/random-effect term, so these do not account for
# the repeated sampling of the same 6 pigs across all 8 conditions. Treat the
# pig-blocked limma results above as primary; use these as a complementary,
# bootstrap-uncertainty-aware check.

so2 <- sleuth_prep(s2c, target_mapping = t2g, aggregation_column = "ens_gene",
                    extra_bootstrap_summary = TRUE, read_bootstrap_tpm = TRUE)

so2 <- sleuth_fit(so2, ~ site + time + infection, "full")

# main effect of infection, controlling for site & time
so2 <- sleuth_fit(so2, ~ site + time, "reduced_infection")
so2 <- sleuth_lrt(so2, "reduced_infection", "full")
sr_infection <- sleuth_results(so2, "reduced_infection:full", "lrt", show_all = FALSE)
sr_infection_sig <- filter(sr_infection, qval <= 0.05)
write.csv(sr_infection_sig, "sleuth_DE_infection_main_effect.csv", row.names = FALSE)
cat("Genes with a main effect of infection:", nrow(sr_infection_sig), " - 0 \n") 

# main effect of site, controlling for time & infection
so2 <- sleuth_fit(so2, ~ time + infection, "reduced_site")
so2 <- sleuth_lrt(so2, "reduced_site", "full")
sr_site <- sleuth_results(so2, "reduced_site:full", "lrt", show_all = FALSE)
sr_site_sig <- filter(sr_site, qval <= 0.05)
write.csv(sr_site_sig, "sleuth_DE_site_main_effect.csv", row.names = FALSE)
cat("Genes with a main effect of site:", nrow(sr_site_sig), " - 10,517\n")

# main effect of time, controlling for site & infection
so2 <- sleuth_fit(so2, ~ site + infection, "reduced_time")
so2 <- sleuth_lrt(so2, "reduced_time", "full")
sr_time <- sleuth_results(so2, "reduced_time:full", "lrt", show_all = FALSE)
sr_time_sig <- filter(sr_time, qval <= 0.05)
write.csv(sr_time_sig, "sleuth_DE_time_main_effect.csv", row.names = FALSE)
cat("Genes with a main effect of time:", nrow(sr_time_sig), " - 11,101\n")

# Extension: fit the full factorial (full2 = ~site*time*infection), then test
# the 3-way interaction first, and tests each 2-way interaction
# (site×infection, time×infection, site×time) by dropping it from the no-3-way
# model rather than from full2 — testing a 2-way term in the presence of the
# 3-way interaction that contains it would violate marginality. The
# site×infection and time×infection tests are the sleuth-side counterparts of
# your Infection_x_Site/Infection_x_Time limma contrasts


# interaction terms: full factorial model
so2 <- sleuth_fit(so2, ~ site * time * infection, "full2")
 
# 3-way interaction (site x time x infection): does the site-dependence of
# the infection effect itself change between 6 and 24 hpi? Drop the 3-way
# term, keep every 2-way interaction and main effect.
so2 <- sleuth_fit(so2, ~ site + time + infection + site:time + site:infection + time:infection, "no_3way")
so2 <- sleuth_lrt(so2, "no_3way", "full2")
sr_3way <- sleuth_results(so2, "no_3way:full2", "lrt", show_all = FALSE)
sr_3way_sig <- filter(sr_3way, qval <= 0.05)
write.csv(sr_3way_sig, "sleuth_DE_site_x_time_x_infection.csv", row.names = FALSE)
cat("Genes with a site x time x infection interaction:", nrow(sr_3way_sig), "\n")
 
# 2-way interactions, each tested by dropping just that term from the
# no-3-way model (not from full2) -- testing a 2-way term in the presence
# of the 3-way interaction that contains it violates marginality.
so2 <- sleuth_fit(so2, ~ site + time + infection + site:time + time:infection, "no_site_infection")
so2 <- sleuth_lrt(so2, "no_site_infection", "no_3way")
sr_site_infection <- sleuth_results(so2, "no_site_infection:no_3way", "lrt", show_all = FALSE)
sr_site_infection_sig <- filter(sr_site_infection, qval <= 0.05)
write.csv(sr_site_infection_sig, "sleuth_DE_site_x_infection.csv", row.names = FALSE)
cat("Genes where the infection effect differs by site:", nrow(sr_site_infection_sig), "\n")
 
so2 <- sleuth_fit(so2, ~ site + time + infection + site:time + site:infection, "no_time_infection")
so2 <- sleuth_lrt(so2, "no_time_infection", "no_3way")
sr_time_infection <- sleuth_results(so2, "no_time_infection:no_3way", "lrt", show_all = FALSE)
sr_time_infection_sig <- filter(sr_time_infection, qval <= 0.05)
write.csv(sr_time_infection_sig, "sleuth_DE_time_x_infection.csv", row.names = FALSE)
cat("Genes where the infection effect differs by time:", nrow(sr_time_infection_sig), "\n")
 
so2 <- sleuth_fit(so2, ~ site + time + infection + site:infection + time:infection, "no_site_time")
so2 <- sleuth_lrt(so2, "no_site_time", "no_3way")
sr_site_time <- sleuth_results(so2, "no_site_time:no_3way", "lrt", show_all = FALSE)
sr_site_time_sig <- filter(sr_site_time, qval <= 0.05)
write.csv(sr_site_time_sig, "sleuth_DE_site_x_time.csv", row.names = FALSE)
cat("Genes where the site effect differs by time:", nrow(sr_site_time_sig), "\n")

#####################  PCA #############

# Get TPM values, filter, make PCA plots

tpm_data <- so2$obs_norm %>% # extract TPM values
  dplyr::select(target_id, sample, tpm) %>%
  tidyr::pivot_wider(names_from = sample, values_from = tpm)
write.csv(tpm_data, "tpm_values.csv", row.names = FALSE)

# Convert to matrix format for PCA
tpm_matrix <- as.matrix(tpm_data[,-1])
rownames(tpm_matrix) <- tpm_data$target_id
var_per_gene <- apply(tpm_matrix, 1, var) # filter for non-zero values
tpm_matrix_filtered <- tpm_matrix[var_per_gene > 0, ]
pca_result <- prcomp(t(tpm_matrix_filtered), scale = T)
pca_df <- as.data.frame(pca_result$x)
pca_df$sample <- rownames(pca_df)
pca_df$sample <- rownames(pca_df)
pca_df <- left_join(pca_df, so2$sample_to_covariates, by = "sample")
var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)

pca_df$site <- factor(pca_df$site, levels = c("bronchial", "tracheal"))
pca_df$infection <- factor(pca_df$infection, levels = c("mock", "infected"))
pca_df$time <- factor(pca_df$time, levels = c("6", "24"))
pca_df$group <- interaction(pca_df$site, pca_df$infection, pca_df$time, sep = ":")
group_order <- c("bronchial:infected:24", "bronchial:mock:24",
                  "bronchial:infected:6",  "bronchial:mock:6",
                  "tracheal:infected:24",  "tracheal:mock:24",
                  "tracheal:infected:6",   "tracheal:mock:6")
pca_df$group <- factor(pca_df$group, levels = group_order)
group_levels <- levels(pca_df$group)

print("Doing PCA")

custom_colors <- c(
  "yellow", # b inf 24 
  "cyan", # b m 24
  "orange", # Blb inf 6
  "blue", # b m 6
  "red", # t inf 24
  "navy", # t m 24
  "brown", # t inf 6 Purple
  "grey" ) # t m 6

custom_shapes <- c(
  16, # Solid Circle
  1, # Solid Triangle
  16, # Solid Square
  2, # Solid Diamond
  17,  # Open Circle
  0,  # Open Triangle
  17,  # Open Square
  3  ) # Open Diamond
my_colors <- setNames(custom_colors[1:length(group_levels)], group_levels)
my_shapes <- setNames(custom_shapes[1:length(group_levels)], group_levels) 

pdf("PC1.PC2.pdf", width=11, height=5)
p <-  ggplot(pca_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 5, alpha = 1) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(   x = paste0("PC1 (", var_explained[1], "%)"), 
    y = paste0("PC2 (", var_explained[2], "%)"),
    color = "Group", shape = "Group" ) +
  theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 21))
print(p)
dev.off()
ggsave("PC1.PC2.png", plot = p, width =11, height = 5, dpi = 300)

pdf("PC3.PC4.pdf", width=8, height=5)
p2 <- ggplot(pca_df, aes(x = PC3, y = PC4, color = group, shape = group)) +
  geom_point(size = 5, alpha = 1) +
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = my_shapes) +
  labs(x = paste0("PC3 (", var_explained[3], "%)"), y = paste0("PC4 (", var_explained[4], "%)")) +
  theme_bw() +
  theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20),
        legend.text = element_text(size = 20), legend.title = element_text(size = 21))
print(p2)
dev.off()
ggsave("PC3.PC4.png", plot = p2, width = 8, height = 5, dpi = 300)

pdf("plot_sample_heatmap.pdf", width=11, height=11)
pp <- plot_sample_heatmap(so2,use_filtered = T, color_high = "white",
  color_low = "dodgerblue", x_axis_angle = 90,
  #annotation_cols=setdiff(colnames(so$sample_to_covariates), "sample"),
  cluster_bool = T)
print(pp)
dev.off()
ggsave("plot_sample_heatmap.png", plot = pp, width =10, height =10, dpi = 300)

####### 6 - STRING protein-protein interaction networks ########################

dir.create("STRING", showWarnings = F)

string_db <- STRINGdb$new(version = "12", species = 9823, score_threshold =400)

run_string_ppi <- function(de_df, comparison_name, min_hub_degree = 10) {
  genes_for_string <- data.frame(gene = unique(na.omit(de_df$gene_label)))
  if (nrow(genes_for_string) < 1) {
    cat(comparison_name, ": too few DE genes for a PPI network, skipping\n")
    return(invisible(NULL))
  }

  mapped <- string_db$map(genes_for_string, "gene")
  cat(comparison_name, "- input genes:", nrow(genes_for_string),
      "mapped genes:", nrow(mapped), "\n")

  ids <- unique(mapped$STRING_id)
  ppi <- string_db$get_interactions(ids)
  ppi <- dplyr::filter(ppi, from %in% ids, to %in% ids)

  lookup <- mapped %>%
    dplyr::select(STRING_id, gene) %>%
    dplyr::left_join(de_df %>% dplyr::select(gene_label, logFC) %>% dplyr::distinct(),
                      by = c("gene" = "gene_label")) %>%
    dplyr::mutate(direction = ifelse(logFC > 0, "Upregulated", "Downregulated"))

  edges <- ppi %>%
    dplyr::left_join(lookup, by = c("from" = "STRING_id")) %>%
    dplyr::rename(gene1 = gene) %>%
    dplyr::left_join(lookup, by = c("to" = "STRING_id")) %>%
    dplyr::rename(gene2 = gene)

  g <- graph_from_data_frame(edges[, c("gene1", "gene2")], directed = FALSE)
  node_info <- lookup %>% dplyr::select(gene, direction) %>% dplyr::distinct()
  V(g)$direction <- node_info$direction[match(V(g)$name, node_info$gene)]
  V(g)$degree <- degree(g)

  write.csv(data.frame(gene = V(g)$name, degree = V(g)$degree),
            paste0("STRING/", comparison_name, "_nodes.csv"), row.names = FALSE)

  comp <- components(g)
  g_main <- induced_subgraph(g, V(g)[comp$membership == which.max(comp$csize)])
  V(g_main)$degree <- degree(g_main)

  p <- ggraph(g_main, layout = "fr") +
    geom_edge_link(colour = "grey70", alpha = 0.5, width = 1.2) +
    geom_node_point(aes(size = degree, colour = direction), alpha = 0.8) +
    geom_node_text(aes(label = name), repel = TRUE, size = 5, max.overlaps = Inf) +
    scale_colour_manual(values = c(Upregulated = "#2166AC", Downregulated = "#B2182B")) +
    theme_void() +
    labs(title = paste0(comparison_name, " -- largest connected component"))
  ggsave(paste0("STRING/", comparison_name, "_PPIN.pdf"), p, width = 7, height = 6)
  ggsave(paste0("STRING/", comparison_name, "_PPIN.png"), p, width = 7, height = 6, dpi = 300)

  g_hubs <- induced_subgraph(g_main, V(g_main)[degree(g_main) >= min_hub_degree])
  if (vcount(g_hubs) > 0) {
    p_hubs <- ggraph(g_hubs, layout = "fr") +
      geom_edge_link(colour = "grey70", alpha = 0.5, width = 1.2) +
      geom_node_point(aes(size = degree, colour = direction), alpha = 0.8) +
      geom_node_text(aes(label = name), repel = TRUE, size = 5, max.overlaps = Inf) +
      scale_colour_manual(values = c(Upregulated = "#2166AC", Downregulated = "#B2182B")) +
      theme_void() +
      labs(title = paste0(comparison_name, " -- hub genes (degree >= ", min_hub_degree, ")"))
    ggsave(paste0("STRING/", comparison_name, "_PPIN_hubs.pdf"), p_hubs, width = 6, height = 5)
    ggsave(paste0("STRING/", comparison_name, "_PPIN_hubs.png"), p_hubs, width = 6, height = 5, dpi = 300)
  }
  invisible(g_main)
}

string_gsea_contrasts <- c(primary_contrasts, "Infection_overall")
for (cn in string_gsea_contrasts) {
  run_string_ppi(subset(vol_list[[cn]], significant == "Significant"), cn) }

####### 7 - KEGG GSEA ###########################################################

dir.create("KEGG", showWarnings = FALSE)

run_kegg_gsea <- function(res_df, comparison_name, t2g) {
  res_df$target_id <- rownames(res_df)
  res_df$target_id_noversion <- sub("\\..*$", "", res_df$target_id)

  t2g_lookup <- t2g %>%
    dplyr::mutate(target_id_noversion = sub("\\..*$", "", target_id)) %>%
    dplyr::select(target_id_noversion, ens_gene, ext_gene) %>%
    dplyr::distinct(target_id_noversion, .keep_all = TRUE)

  res_df <- dplyr::left_join(res_df, t2g_lookup, by = "target_id_noversion")

  # one signed score per gene: mean t-statistic across its transcripts
  gene_rank <- res_df %>%
    dplyr::filter(!is.na(ens_gene), ens_gene != "") %>%
    dplyr::group_by(ens_gene) %>%
    dplyr::summarise(rank = mean(t, na.rm = TRUE), .groups = "drop")

  gene_symbols <- t2g %>% dplyr::select(ens_gene, ext_gene) %>% dplyr::distinct(ens_gene, .keep_all = TRUE)
  gene_rank <- dplyr::left_join(gene_rank, gene_symbols, by = "ens_gene")

  entrez_map <- AnnotationDbi::select(
    org.Ss.eg.db, keys = unique(na.omit(gene_rank$ext_gene)),
    keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID")
  ) %>%
    dplyr::rename(ext_gene = SYMBOL, entrezgene_id = ENTREZID) %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::distinct(ext_gene, .keep_all = TRUE)

  gene_rank <- gene_rank %>%
    dplyr::inner_join(entrez_map, by = "ext_gene") %>%
    dplyr::filter(!is.na(entrezgene_id)) %>%
    dplyr::distinct(entrezgene_id, .keep_all = TRUE)

  cat(comparison_name, "- genes with Entrez IDs mapped:", nrow(gene_rank), "\n")

  ranks <- gene_rank$rank
  names(ranks) <- gene_rank$entrezgene_id
  ranks <- sort(ranks, decreasing = TRUE)

  gsea <- gseKEGG(geneList = ranks, organism = "ssc", minGSSize = 10,
                   maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE)
  gsea_df <- as.data.frame(gsea)
  write.csv(gsea_df, paste0("KEGG/", comparison_name, "_KEGG_GSEA.csv"), row.names = FALSE)
  cat(comparison_name, "- enriched pathways:", nrow(gsea_df), "\n")

  if (nrow(gsea_df) > 0) {
    p <- dotplot(gsea, showCategory = 20) +
      labs(title = comparison_name) +
      theme(legend.title = element_text(size = 12), legend.text = element_text(size = 10))
    ggsave(paste0("KEGG/", comparison_name, "_KEGG_GSEA.pdf"), p, width = 9, height = 7)
    ggsave(paste0("KEGG/", comparison_name, "_KEGG_GSEA.png"), p, width = 9, height = 7, dpi = 300)
  } else {
    message(comparison_name, ": no pathways passed pvalueCutoff = 0.05") }
  invisible(gsea_df) }

for (cn in string_gsea_contrasts) {run_kegg_gsea(res_list[[cn]], cn, t2g) }

# # # power testing

# exact power from your real fitted model
sigma_typical <- median(sqrt(fit$s2.post))       # per-gene posterior residual SD
df_total <- fit$df.residual[1] + fit$df.prior    # eBayes-moderated df

power_calc <- function(n_per_group, sigma, df, alpha = 0.05, delta = 1) {
  ncp <- delta / (sigma * sqrt(2 / n_per_group))
  tcrit <- qt(1 - alpha / 2, df)
  pt(tcrit, df, ncp, lower.tail = FALSE) + pt(-tcrit, df, ncp) }

sapply(c(6, 12, 24), power_calc, sigma = sigma_typical, df = df_total)
