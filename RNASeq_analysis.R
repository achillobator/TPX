## making cluster plots for the infected degs 

# Install BiocManager first if you don't have it
install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "ggplot2",
  "ggrepel",
  "patchwork",
  "dplyr",
  "tidyr",
  "tibble",
  "UpSetR",
  "WGCNA"
))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2",
  "edgeR",
  "limma",
  "pheatmap",
  "ComplexHeatmap",
  "circlize",
  "fgsea",
  "clusterProfiler",
  "org.At.tair.db"   # Arabidopsis annotation db for GO enrichment
))

# Core RNA-seq / stats
library(DESeq2)
library(edgeR)          # optional, if you want TMM normalization utilities

# Dimensionality reduction & clustering
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")

# Then reinstall WGCNA
BiocManager::install("WGCNA")

# And try loading again
library(WGCNA)


# Visualization - heatmaps
library(pheatmap)
library(ComplexHeatmap)  # more powerful alternative to pheatmap
library(circlize)        # color scale helper for ComplexHeatmap

# Visualization - general
library(ggplot2)
library(ggrepel)         # non-overlapping labels on PCA/volcano plots
library(patchwork)       # combining multiple ggplots

# Upset plots
library(UpSetR)

# Data wrangling
library(dplyr)
library(tidyr)
library(tibble)          # rownames_to_column / column_to_rownames

# Gene set enrichment (optional but useful)
library(fgsea)
library(clusterProfiler) # GO/KEGG enrichment
library(pheatmap)
library(tidyverse)

getwd()
setwd("~/Documents/UW_nemlab/TPL_TPX_family/Manuscript_files/Figure_2/RNA-Seq_data")

dat <- read.csv("Data/At_mean_norm_counts.csv")

dim(dat)
head(dat[, 1:6])
colnames(dat)

# Clean up the matrix
rownames(dat) <- dat$Gene_ID
dat <- dat[, -c(1,2)]  # drop X and Gene_ID columns

# Build metadata frame
coldata <- data.frame(
  sample    = colnames(dat),
  genotype  = sub("_.*", "", colnames(dat)),
  treatment = ifelse(grepl("_M", colnames(dat)), "Mock", "Botrytis"),
  timepoint = ifelse(grepl("16", colnames(dat)), "T16", "T20"),
  row.names = colnames(dat)
)

coldata$genotype  <- factor(coldata$genotype,  levels = c("WT", "L1", "L14", "T23"))
coldata$treatment <- factor(coldata$treatment, levels = c("Mock", "Botrytis"))
coldata$timepoint <- factor(coldata$timepoint, levels = c("T16", "T20"))

coldata

# Log2+1 transform to stabilize variance
dat_log <- log2(dat + 1)

# Filter low-variance genes (speeds up PCA, reduces noise)
gene_vars <- apply(dat_log, 1, var)
dat_filt  <- dat_log[gene_vars > quantile(gene_vars, 0.25), ]  # keep top 75%

# Run PCA (transpose so samples are rows)
pca_res  <- prcomp(t(dat_filt), scale. = TRUE, center = TRUE)

# Extract variance explained
pct_var  <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

# Build plotting dataframe
pca_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  PC3       = pca_res$x[, 3],
  sample    = rownames(pca_res$x)
) |> left_join(coldata, by = "sample")

# Plot
pca1v3<- ggplot(pca_df, aes(x = PC1, y = PC2,
                   fill     = genotype,
                   shape    = treatment,
                   label    = sample)) +
  geom_point(size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  ))+
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC2 (", pct_var[2], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA — Normalized counts"
  ) +
  theme_classic()+
  theme(
    legend.position   = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )
pca1v3

### trying a PCA with lines:

# Add a line group variable connecting T16 to T20
pca_mod_df <- pca_mod_df %>%
  mutate(line_group = paste(genotype, treatment, sep = "_"))

ggplot(pca_mod_df, aes(x = PC1, y = PC2)) +
  geom_line(aes(group = line_group, color = genotype),
            linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  geom_point(aes(fill = genotype, shape = treatment),
             size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#ff00ff",
    "T23" = "#1bd218"
  )) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#ff00ff",
    "T23" = "#1bd218"
  )) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21)),
    color = "none"
  ) +
  labs(
    x     = paste0("PC1 (", pct_var_mod[1], "%)"),
    y     = paste0("PC2 (", pct_var_mod[2], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA - module genes only (n = 786)"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

# this is on the whole data set not just the modules
pca_df <- pca_df %>%
  mutate(line_group = paste(genotype, treatment, sep = "_"))

ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_line(aes(group = line_group, color = genotype),
            linewidth = 0.6, linetype = "dashed", alpha = 0.7) +
  geom_point(aes(fill = genotype, shape = treatment),
             size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#ff00ff",
    "T23" = "#1bd218"
  )) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#ff00ff",
    "T23" = "#1bd218"
  )) +
  guides(
    fill  = guide_legend(override.aes = list(shape = 21)),
    color = "none"
  ) +
  labs(
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC2 (", pct_var[2], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA - full dataset"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

pdf("PCA1v2_lines_all.pdf", width = 4, height = 4)
pca1v3
dev.off()

## PC1 v PC3
ggplot(pca_df, aes(x = PC1, y = PC3,
                   fill     = genotype,
                   shape    = treatment,
                   label    = sample)) +
  geom_point(size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x     = paste0("PC1 (", pct_var[1], "%)"),
    y     = paste0("PC3 (", pct_var[3], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA — PC1 vs PC3"
  ) +
  theme_classic()

##
library(pheatmap)

# Take top 500 most variable genes
top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])
dat_top <- dat_log[top_var_genes, ]

# Z-score scale by row
dat_scaled <- t(scale(t(dat_top)))

# Annotation dataframe for columns
ann_col <- coldata[, c("genotype", "treatment", "timepoint")]

# Annotation colors
ann_colors <- list(
  genotype  = c("WT" = "#999999", "L1" = "#197ad4", "L14" = "#8f00ff", "T23" = "#1bd218"),
  treatment = c("Mock" = "white", "Botrytis" = "#4a0000"),
  timepoint = c("T16" = "#f5c842", "T20" = "#7b2d8b")
)

pheatmap(dat_scaled,
         annotation_col  = ann_col,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method        = "ward.D2",
         show_rownames   = FALSE,
         show_colnames   = TRUE,
         color           = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main            = "Top 500 variable genes — z-score"
)

library(limma)

# Design matrix with combined group factor
coldata$group <- factor(paste(coldata$genotype, coldata$treatment, coldata$timepoint, sep = "_"))
design <- model.matrix(~ 0 + group, data = coldata)
colnames(design) <- levels(coldata$group)

# Fit model
fit <- lmFit(dat_log, design)

# Define contrasts — each mutant Botrytis vs WT Botrytis (at each timepoint)
contrasts_mat <- makeContrasts(
  L1_vs_WT_B16   = L1_Botrytis_T16  - WT_Botrytis_T16,
  L1_vs_WT_B20   = L1_Botrytis_T20  - WT_Botrytis_T20,
  L14_vs_WT_B16  = L14_Botrytis_T16 - WT_Botrytis_T16,
  L14_vs_WT_B20  = L14_Botrytis_T20 - WT_Botrytis_T20,
  T23_vs_WT_B16  = T23_Botrytis_T16 - WT_Botrytis_T16,
  T23_vs_WT_B20  = T23_Botrytis_T20 - WT_Botrytis_T20,
  levels = design
)

fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2)


get_degs <- function(fit, coef, fdr = 0.05, lfc = 1) {
  tt <- topTable(fit, coef = coef, number = Inf, adjust.method = "BH")
  rownames(tt[tt$adj.P.Val < fdr & abs(tt$logFC) > lfc, ])
}

# Split data and metadata by timepoint
dat_T16 <- dat_log[, coldata$timepoint == "T16"]
dat_T20 <- dat_log[, coldata$timepoint == "T20"]

coldata_T16 <- coldata[coldata$timepoint == "T16", ]
coldata_T20 <- coldata[coldata$timepoint == "T20", ]

# Design matrices
design_T16 <- model.matrix(~ 0 + genotype + treatment, data = coldata_T16)
design_T20 <- model.matrix(~ 0 + genotype + treatment, data = coldata_T20)

colnames(design_T16) <- gsub("genotype|treatment", "", colnames(design_T16))
colnames(design_T20) <- gsub("genotype|treatment", "", colnames(design_T20))

# Fit models
fit_T16 <- eBayes(contrasts.fit(lmFit(dat_T16, design_T16),
                                makeContrasts(
                                  L1_vs_WT  = L1  - WT,
                                  L14_vs_WT = L14 - WT,
                                  T23_vs_WT = T23 - WT,
                                  levels = design_T16
                                )
))

fit_T20 <- eBayes(contrasts.fit(lmFit(dat_T20, design_T20),
                                makeContrasts(
                                  L1_vs_WT  = L1  - WT,
                                  L14_vs_WT = L14 - WT,
                                  T23_vs_WT = T23 - WT,
                                  levels = design_T20
                                )
))

# Extract DEGs
deg_list <- list(
  L1_T16  = get_degs(fit_T16, "L1_vs_WT"),
  L14_T16 = get_degs(fit_T16, "L14_vs_WT"),
  T23_T16 = get_degs(fit_T16, "T23_vs_WT"),
  L1_T20  = get_degs(fit_T20, "L1_vs_WT"),
  L14_T20 = get_degs(fit_T20, "L14_vs_WT"),
  T23_T20 = get_degs(fit_T20, "T23_vs_WT")
)

sapply(deg_list, length)


install.packages("umap")
library(umap)

# Run UMAP on the same filtered, log-transformed matrix used for PCA
umap_res <- umap(t(dat_filt))

# Build plotting dataframe
umap_df <- data.frame(
  UMAP1  = umap_res$layout[, 1],
  UMAP2  = umap_res$layout[, 2],
  sample = rownames(umap_res$layout)
) |> left_join(coldata, by = "sample")

# Plot
ggplot(umap_df, aes(x = UMAP1, y = UMAP2,
                    fill  = genotype,
                    shape = treatment,
                    label = sample)) +
  geom_point(size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x     = "UMAP1",
    y     = "UMAP2",
    fill  = "Genotype",
    shape = "Treatment",
    title = "UMAP — Normalized counts"
  ) +
  theme_classic()

## run timepoits seperately
# Split filtered matrix by timepoint
dat_filt_T16 <- dat_filt[, coldata$timepoint == "T16"]
dat_filt_T20 <- dat_filt[, coldata$timepoint == "T20"]

umap_config <- umap.defaults
umap_config$n_neighbors <- 5  # must be less than n samples (8)
umap_config$min_dist    <- 0.3

umap_T16 <- umap(t(dat_filt_T16), config = umap_config)
umap_T20 <- umap(t(dat_filt_T20), config = umap_config)

# Build plotting dataframes
umap_df_T16 <- data.frame(
  UMAP1  = umap_T16$layout[, 1],
  UMAP2  = umap_T16$layout[, 2],
  sample = rownames(umap_T16$layout)
) |> left_join(coldata, by = "sample")

umap_df_T20 <- data.frame(
  UMAP1  = umap_T20$layout[, 1],
  UMAP2  = umap_T20$layout[, 2],
  sample = rownames(umap_T20$layout)
) |> left_join(coldata, by = "sample")

# Plot function to avoid repeating code
plot_umap <- function(df, title) {
  ggplot(df, aes(x = UMAP1, y = UMAP2,
                 fill  = genotype,
                 shape = treatment,
                 label = sample)) +
    geom_point(size = 4, color = "black", stroke = 0.5) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
    scale_fill_manual(values = c(
      "WT"  = "#999999",
      "L1"  = "#197ad4",
      "L14" = "#8f00ff",
      "T23" = "#1bd218"
    )) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = "UMAP1", y = "UMAP2",
         fill  = "Genotype",
         shape = "Treatment",
         title = title) +
    theme_classic()
}

# Plot side by side
p_T16 <- plot_umap(umap_df_T16, "UMAP — T16 hpi")
p_T20 <- plot_umap(umap_df_T20, "UMAP — T20 hpi")

library(patchwork)
p_T16 + p_T20


## pull in degs
bc <- read.csv("Data/bcvsmock_allmuts_sigdegs.csv")
mut <- read.csv("Data/mutsvsWT_sigdegs_LFC0.5.csv")

head(bc)
head(mut)
colnames(bc)
colnames(mut)

# Check unique groupings
unique(bc[, c("Genotype", "Timepoint")])
unique(mut$contrast)

##build Venn diagrams:
install.packages("ggVennDiagram")
library(ggVennDiagram)

# Standardize contrast names in mut to match bc
mut$Genotype <- dplyr::recode(mut$contrast,
                              "tpltpr1"     = "L1",
                              "tpltpr1tpr4" = "L14",
                              "tpr2tpr3"    = "T23"
)

# --- Per-genotype lists (mock DEGs vs Botrytis DEGs) ---
genotypes <- c("L1", "L14", "T23")

per_geno_lists <- lapply(genotypes, function(g) {
  list(
    Mock_vs_WT  = mut$Gene_ID[mut$Genotype == g],
    Botrytis_T16 = bc$Gene_ID[bc$Genotype == g & bc$Timepoint == "T16"],
    Botrytis_T20 = bc$Gene_ID[bc$Genotype == g & bc$Timepoint == "T20"]
  )
})
names(per_geno_lists) <- genotypes

# --- Across-genotype lists split by condition ---
across_mock <- list(
  L1  = mut$Gene_ID[mut$Genotype == "L1"],
  L14 = mut$Gene_ID[mut$Genotype == "L14"],
  T23 = mut$Gene_ID[mut$Genotype == "T23"]
)

across_bc_T16 <- list(
  L1  = bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T16"],
  L14 = bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T16"],
  T23 = bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T16"]
)

across_bc_T20 <- list(
  L1  = bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T20"],
  L14 = bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T20"],
  T23 = bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T20"]
)

# --- Per-genotype Venns ---
venn_colors <- c("#197ad4", "#b2182b", "#d95f02")  # mock, T16, T20

p_per_geno <- lapply(genotypes, function(g) {
  ggVennDiagram(per_geno_lists[[g]],
                label_alpha = 0,
                label = "count") +
    scale_fill_gradient(low = "white", high = venn_colors[1]) +
    ggtitle(paste(g, "— Mock vs WT / Botrytis T16 / Botrytis T20")) +
    theme(legend.position = "none")
})

p_per_geno[[1]] + p_per_geno[[2]] + p_per_geno[[3]] +
  plot_layout(ncol = 3)

# --- Across-genotype Venns ---
col_L1  <- "#197ad4"
col_L14 <- "#8f00ff"
col_T23 <- "#1bd218"

p_across_mock <- ggVennDiagram(across_mock,
                               label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Mock DEGs — across genotypes") +
  theme(legend.position = "none")

p_across_T16 <- ggVennDiagram(across_bc_T16,
                              label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Botrytis DEGs T16 — across genotypes") +
  theme(legend.position = "none")

p_across_T20 <- ggVennDiagram(across_bc_T20,
                              label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Botrytis DEGs T20 — across genotypes") +
  theme(legend.position = "none")

p_across_mock + p_across_T16 + p_across_T20 +
  plot_layout(ncol = 3)

# Per-genotype Venns
# Rebuild per-genotype plots with plain hyphen
p_per_geno <- lapply(genotypes, function(g) {
  ggVennDiagram(per_geno_lists[[g]],
                label_alpha = 0,
                label = "count") +
    scale_fill_gradient(low = "white", high = venn_colors[1]) +
    ggtitle(paste(g, "- Mock vs WT / Botrytis T16 / Botrytis T20")) +
    theme(legend.position = "none")
})

# Rebuild across-genotype plots
p_across_mock <- ggVennDiagram(across_mock,
                               label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Mock DEGs - across genotypes") +
  theme(legend.position = "none")

p_across_T16 <- ggVennDiagram(across_bc_T16,
                              label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Botrytis DEGs T16 - across genotypes") +
  theme(legend.position = "none")

p_across_T20 <- ggVennDiagram(across_bc_T20,
                              label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "grey60") +
  ggtitle("Botrytis DEGs T20 - across genotypes") +
  theme(legend.position = "none")

# Re-export
pdf("per_genotype_venns.pdf", width = 15, height = 6)
p_per_geno[[1]] + p_per_geno[[2]] + p_per_geno[[3]] +
  plot_layout(ncol = 3)
dev.off()

pdf("across_genotype_venns.pdf", width = 15, height = 6)
p_across_mock + p_across_T16 + p_across_T20 +
  plot_layout(ncol = 3)
dev.off()

### bar graphs
summary <- bc %>%
  group_by(Genotype, Timepoint, DE_Direction) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(Genotype = factor(Genotype, levels = c("WT", "L1", "L14", "T23")))

bargraph <- ggplot(summary, aes(x = Timepoint, y = count, fill = DE_Direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~Genotype, nrow = 1) +
  scale_y_continuous(expand = c(0,0), limit = c(0, 2550)) +
  scale_fill_manual(values = c("Downregulated" = "#e28743", "Upregulated" = "#1e81b0")) +
  labs(y = "Number of DEGs (infected vs mock)",
       fill = "DE Direction") +
  theme_bw() +
  theme(text = element_text(size = 6),
        legend.position = c(0.99, 0.99),
        legend.justification = c("right", "top"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        legend.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        strip.background = element_blank())
bargraph
bargraph
pdf("bargraph_DEGs.pdf", width = 4, height = 4)
bargraph
dev.off()

### pull in modules
modules <- read.csv("Data/mutsvsWT_newclusters.csv")

head(modules)
colnames(modules)
dim(modules)

colnames(modules) <- c("Gene_ID", "cluster")
table(modules$cluster)


library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Step 1 — join modules to normalized counts
dat_counts <- dat %>%
  rownames_to_column("Gene_ID")

degs_counts <- modules %>%
  left_join(dat_counts, by = "Gene_ID")

# Step 2 — pivot to long format
degs_longer <- degs_counts %>%
  pivot_longer(cols = -c(Gene_ID, cluster),
               names_to = c("genotype", "condition"),
               names_sep = "_",
               values_to = "count")

# Split condition into treatment (B/M) and timepoint (16/20)
degs_longer <- degs_longer %>%
  separate(col = condition, into = c("treatment", "timepoint"), sep = 1) %>%
  mutate(
    timepoint = as.numeric(timepoint),
    genotype  = factor(genotype, levels = c("WT", "L1", "L14", "T23")),
    treatment = case_when(
      treatment == "B" ~ "Botrytis",
      treatment == "M" ~ "Mock"
    )
  )

# Step 3 — split by cluster
cl_nums <- as.character(1:11)

cl_num_lists <- lapply(seq_along(cl_nums), function(n)
  filter(modules, cluster == n))
names(cl_num_lists) <- cl_nums

# Step 4 — join counts to each cluster
cl_num_lists_counts <- lapply(seq_along(cl_num_lists), function(x)
  left_join(cl_num_lists[[x]], degs_longer, by = "Gene_ID"))
names(cl_num_lists_counts) <- cl_nums

# Step 5 — z-score scale per gene
cl_num_lists_scaled <- lapply(seq_along(cl_num_lists_counts), function(s)
  cl_num_lists_counts[[s]] %>%
    group_by(Gene_ID) %>%
    mutate(scaled_expression = (count - mean(count)) / sd(count)) %>%
    ungroup())

# Step 6 — summarise per genotype/treatment/timepoint
cl_num_lists_sum <- lapply(seq_along(cl_num_lists_scaled), function(k)
  cl_num_lists_scaled[[k]] %>%
    group_by(genotype, treatment, timepoint) %>%
    summarise(
      mean = mean(scaled_expression, na.rm = TRUE),
      n    = n(),
      sd   = sd(scaled_expression, na.rm = TRUE),
      se   = sd / sqrt(n),
      ci   = 1.96 * se,
      .groups = "drop"
    ))

# Cluster names with n sizes from our table()
cl_sizes <- table(modules$cluster)
k_names <- paste0("Cluster ", 1:11, " (n = ", as.vector(cl_sizes), ")")
names(cl_num_lists_sum) <- k_names

# Step 7 — combine all clusters for facet grid plot
cl_all <- lapply(seq_along(cl_num_lists_sum), function(x)
  cl_num_lists_sum[[x]] %>%
    mutate(cluster_name = k_names[x]))

cl_rowbind <- list_rbind(cl_all) %>%
  mutate(cluster_name = factor(cluster_name, levels = k_names))

# Step 8 — plot all clusters
ggplot(data = cl_rowbind, aes(x = timepoint, y = mean,
                              colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment")

# going back to PCA

# Subset the filtered log matrix to just the 786 module genes
module_genes <- modules$Gene_ID

# Keep only module genes that are present in dat_log
module_genes_present <- module_genes[module_genes %in% rownames(dat_log)]

dat_modules <- dat_log[module_genes_present, ]

# Rerun PCA
pca_mod <- prcomp(t(dat_modules), scale. = TRUE, center = TRUE)

# Variance explained
pct_var_mod <- round(100 * pca_mod$sdev^2 / sum(pca_mod$sdev^2), 1)

# Build plot dataframe
pca_mod_df <- data.frame(
  PC1    = pca_mod$x[, 1],
  PC2    = pca_mod$x[, 2],
  PC3    = pca_mod$x[, 3],
  sample = rownames(pca_mod$x)
) |> left_join(coldata, by = "sample")

# Plot
ggplot(pca_mod_df, aes(x = PC1, y = PC2,
                       fill  = genotype,
                       shape = treatment,
                       label = sample)) +
  geom_point(size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x     = paste0("PC1 (", pct_var_mod[1], "%)"),
    y     = paste0("PC2 (", pct_var_mod[2], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA - module genes only (n = 786)"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

# facet timepoint:
ggplot(pca_mod_df, aes(x = PC1, y = PC2,
                       fill  = genotype,
                       shape = treatment,
                       label = sample)) +
  geom_point(size = 4, color = "black", stroke = 0.5) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~timepoint, nrow = 1) +
  labs(
    x     = paste0("PC1 (", pct_var_mod[1], "%)"),
    y     = paste0("PC2 (", pct_var_mod[2], "%)"),
    fill  = "Genotype",
    shape = "Treatment",
    title = "PCA - module genes only (n = 786)"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("left", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

##
# Extract top gene loadings on PC1 and PC2
loadings <- data.frame(
  Gene_ID = rownames(pca_mod$rotation),
  PC1     = pca_mod$rotation[, 1],
  PC2     = pca_mod$rotation[, 2]
) %>%
  left_join(modules, by = "Gene_ID") %>%
  mutate(cluster = factor(cluster))

# Keep top 20 genes by overall loading magnitude
top_loadings <- loadings %>%
  mutate(magnitude = sqrt(PC1^2 + PC2^2)) %>%
  slice_max(magnitude, n = 20)

# Scaling factor to make arrows visible alongside sample points
scale_factor <- 10

# Plot
ggplot() +
  # Sample points
  geom_point(data = pca_mod_df,
             aes(x = PC1, y = PC2, fill = genotype, shape = treatment),
             size = 4, color = "black", stroke = 0.5) +
  # Loading arrows
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * scale_factor,
                   yend = PC2 * scale_factor,
                   color = cluster),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.6) +
  # Gene labels
  geom_text_repel(data = top_loadings,
                  aes(x = PC1 * scale_factor,
                      y = PC2 * scale_factor,
                      label = Gene_ID,
                      color = cluster),
                  size = 2.5, max.overlaps = 20) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x      = paste0("PC1 (", pct_var_mod[1], "%)"),
    y      = paste0("PC2 (", pct_var_mod[2], "%)"),
    fill   = "Genotype",
    shape  = "Treatment",
    color  = "Cluster",
    title  = "PCA biplot - top 20 gene loadings"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

# Check how many loadings match to modules
nrow(top_loadings)
sum(!is.na(top_loadings$cluster))

# Check if Gene_IDs match between the two
head(rownames(pca_mod$rotation))
head(modules$Gene_ID)

# Any whitespace or character issues?
sum(rownames(pca_mod$rotation) %in% modules$Gene_ID)

# Auto-scale loadings to sample score range
scale_factor <- max(abs(pca_mod_df[, c("PC1", "PC2")])) / 
  max(abs(top_loadings[, c("PC1", "PC2")])) * 0.5

ggplot() +
  geom_point(data = pca_mod_df,
             aes(x = PC1, y = PC2, fill = genotype, shape = treatment),
             size = 4, color = "black", stroke = 0.5) +
  geom_segment(data = top_loadings,
               aes(x = 0, y = 0,
                   xend = PC1 * scale_factor,
                   yend = PC2 * scale_factor,
                   color = cluster),
               arrow = arrow(length = unit(0.2, "cm")),
               linewidth = 0.6) +
  geom_text_repel(data = top_loadings,
                  aes(x = PC1 * scale_factor,
                      y = PC2 * scale_factor,
                      label = Gene_ID,
                      color = cluster),
                  size = 2.5, max.overlaps = 30) +
  scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "L14" = "#8f00ff",
    "T23" = "#1bd218"
  )) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(
    x      = paste0("PC1 (", pct_var_mod[1], "%)"),
    y      = paste0("PC2 (", pct_var_mod[2], "%)"),
    fill   = "Genotype",
    shape  = "Treatment",
    color  = "Cluster",
    title  = "PCA biplot - top 20 gene loadings"
  ) +
  theme_classic() +
  theme(
    legend.position      = c(0.98, 0.98),
    legend.justification = c("right", "top"),
    legend.background    = element_rect(fill = "white", color = "grey80", linewidth = 0.3)
  )

## eigen vector version
# Calculate eigengene (mean scaled expression) per module per sample
eigengenes <- lapply(1:11, function(cl_num) {
  genes <- modules$Gene_ID[modules$cluster == cl_num]
  genes <- genes[genes %in% rownames(dat_log)]
  
  # Scale each gene then average across genes per sample
  scaled <- t(scale(t(dat_log[genes, ])))
  colMeans(scaled, na.rm = TRUE)
})

eigengene_df <- as.data.frame(do.call(rbind, eigengenes))
rownames(eigengene_df) <- paste0("Cluster_", 1:11)

# Transpose so samples are rows
eigengene_df <- as.data.frame(t(eigengene_df)) %>%
  rownames_to_column("sample") %>%
  left_join(pca_mod_df, by = "sample")

# Plot one panel per cluster, color = eigengene value
plots <- lapply(1:11, function(cl_num) {
  col_name <- paste0("Cluster_", cl_num)
  ggplot(eigengene_df, aes(x = PC1, y = PC2,
                           fill = .data[[col_name]],
                           shape = treatment)) +
    geom_point(size = 4, color = "black", stroke = 0.5) +
    scale_shape_manual(values = c("Mock" = 21, "Botrytis" = 24)) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0) +
    labs(
      x     = paste0("PC1 (", pct_var_mod[1], "%)"),
      y     = paste0("PC2 (", pct_var_mod[2], "%)"),
      fill  = "Eigengene",
      shape = "Treatment",
      title = paste0("Cluster ", cl_num,
                     " (n = ", sum(modules$cluster == cl_num), ")")
    ) +
    theme_classic() +
    theme(text = element_text(size = 8),
          legend.position = "right")
})

# Combine all panels
wrap_plots(plots, ncol = 4)

# Compare drivers: genotype vs treatment vs timepoint
library(broom)

eigen_cols <- paste0("Cluster_", 1:11)

factors <- c("genotype", "treatment", "timepoint")

factor_results <- lapply(factors, function(fac) {
  lapply(eigen_cols, function(cl) {
    aov_res <- aov(eigengene_df[[cl]] ~ eigengene_df[[fac]])
    tidy(aov_res) %>%
      filter(term != "Residuals") %>%
      mutate(cluster = cl, factor = fac)
  }) %>% list_rbind()
}) %>%
  list_rbind() %>%
  select(factor, cluster, statistic, p.value) %>%
  mutate(F_stat = round(statistic, 2),
         p.value = signif(p.value, 3)) %>%
  arrange(factor, p.value)

factor_results
print(factor_results, n = 33)

#######
# Genotype-driven clusters
genotype_clusters <- c(2, 11, 4)
genotype_names <- paste0("Cluster ", genotype_clusters, 
                         " (n = ", as.vector(cl_sizes[genotype_clusters]), ")")

cl_geno <- lapply(seq_along(genotype_clusters), function(x) {
  cl_num_lists_sum[[genotype_clusters[x]]] %>%
    mutate(cluster_name = factor(genotype_names[x], levels = genotype_names))
}) %>% list_rbind()

p_genotype <- ggplot(cl_geno, aes(x = timepoint, y = mean,
                                  colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Genotype-driven clusters")

# Timepoint-driven clusters
timepoint_clusters <- c(10, 7, 3)
timepoint_names <- paste0("Cluster ", timepoint_clusters,
                          " (n = ", as.vector(cl_sizes[timepoint_clusters]), ")")

cl_time <- lapply(seq_along(timepoint_clusters), function(x) {
  cl_num_lists_sum[[timepoint_clusters[x]]] %>%
    mutate(cluster_name = factor(timepoint_names[x], levels = timepoint_names))
}) %>% list_rbind()

p_timepoint <- ggplot(cl_time, aes(x = timepoint, y = mean,
                                   colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Timepoint-driven clusters")

p_genotype
p_timepoint

pdf("genotype_driven_clusters.pdf", width = 8, height = 6)
p_genotype
dev.off()

pdf("timepoint_driven_clusters.pdf", width = 8, height = 6)
p_timepoint
dev.off()

exists("dat_scaled")
exists("ann_col")


# Gene sets for each cluster group
geno_genes <- modules$Gene_ID[modules$cluster %in% genotype_clusters]
time_genes  <- modules$Gene_ID[modules$cluster %in% timepoint_clusters]

# Subset scaled matrix
dat_scaled_geno <- dat_scaled[rownames(dat_scaled) %in% geno_genes, ]
dat_scaled_time <- dat_scaled[rownames(dat_scaled) %in% time_genes, ]

# Add cluster annotation for rows
row_ann_geno <- modules %>%
  filter(cluster %in% genotype_clusters) %>%
  mutate(cluster = factor(cluster)) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

row_ann_time <- modules %>%
  filter(cluster %in% timepoint_clusters) %>%
  mutate(cluster = factor(cluster)) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

# Align row annotation to matrix rows
row_ann_geno <- row_ann_geno[rownames(dat_scaled_geno), , drop = FALSE]
row_ann_time <- row_ann_time[rownames(dat_scaled_time), , drop = FALSE]

# Check overlap between scaled matrix rownames and module genes
sum(geno_genes %in% rownames(dat_scaled))
sum(time_genes %in% rownames(dat_scaled))

# Check what rownames look like in dat_scaled
head(rownames(dat_scaled))
head(geno_genes)

# Rebuild scaled matrices from dat_log for each cluster group
dat_scaled_geno <- t(scale(t(dat_log[rownames(dat_log) %in% geno_genes, ])))
dat_scaled_time <- t(scale(t(dat_log[rownames(dat_log) %in% time_genes, ])))

# Check dimensions
dim(dat_scaled_geno)
dim(dat_scaled_time)

# Genotype-driven heatmap
pheatmap(dat_scaled_geno,
         annotation_col    = ann_col,
         annotation_row    = row_ann_geno,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method        = "ward.D2",
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)

# Timepoint-driven heatmap
pheatmap(dat_scaled_time,
         annotation_col    = ann_col,
         annotation_row    = row_ann_time,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method        = "ward.D2",
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)

head(row_ann_geno)
head(row_ann_time)

# Check rownames align
all(rownames(row_ann_geno) == rownames(dat_scaled_geno))
all(rownames(row_ann_time) == rownames(dat_scaled_time))

# Check why row_ann_time is empty
head(modules %>% filter(cluster %in% timepoint_clusters))
sum(modules$cluster %in% timepoint_clusters)

# Check dimension mismatch for geno
nrow(row_ann_geno)
nrow(dat_scaled_geno)

# Check what's in timepoint_clusters
timepoint_clusters

# Rebuild row annotations cleanly
row_ann_geno <- modules %>%
  filter(cluster %in% genotype_clusters) %>%
  filter(Gene_ID %in% rownames(dat_scaled_geno)) %>%
  mutate(cluster = factor(cluster)) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

row_ann_time <- modules %>%
  filter(cluster %in% timepoint_clusters) %>%
  filter(Gene_ID %in% rownames(dat_scaled_time)) %>%
  mutate(cluster = factor(cluster)) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

# Check
nrow(row_ann_geno)
nrow(row_ann_time)

all(rownames(row_ann_geno) %in% rownames(dat_scaled_geno))
all(rownames(row_ann_time) %in% rownames(dat_scaled_time))

# Genotype-driven heatmap
pheatmap(dat_scaled_geno,
         annotation_col    = ann_col,
         annotation_row    = row_ann_geno,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method        = "ward.D2",
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)

# Timepoint-driven heatmap
pheatmap(dat_scaled_time,
         annotation_col    = ann_col,
         annotation_row    = row_ann_time,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method        = "ward.D2",
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)

sum(is.na(dat_scaled_geno))
sum(is.na(dat_scaled_time))

pheatmap(dat_scaled_geno,
         show_rownames = FALSE,
         show_colnames = TRUE
)

# Reset the graphics device
dev.off()

# Try a basic plot to confirm plotting works at all
plot(1:10)

# Add column annotation
pheatmap(dat_scaled_geno,
         annotation_col = ann_col,
         show_rownames  = FALSE,
         show_colnames  = TRUE
)

pheatmap(dat_scaled_geno,
         annotation_col    = ann_col,
         annotation_row    = row_ann_geno,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)

pheatmap(dat_scaled_time,
         annotation_col    = ann_col,
         annotation_row    = row_ann_time,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)


# Genotype-driven clusters (2, 11, 4)
ann_colors_geno <- ann_colors
ann_colors_geno$Cluster <- c(
  "2"  = "#d9d9d9",
  "4"  = "#969696",
  "11" = "#525252"
)

# Timepoint-driven clusters (10, 7, 3)
ann_colors_time <- ann_colors
ann_colors_time$Cluster <- c(
  "3"  = "#d9d9d9",
  "7"  = "#969696",
  "10" = "#525252"
)

# Genotype-driven
pheatmap(dat_scaled_geno,
         annotation_col    = ann_col,
         annotation_row    = row_ann_geno,
         annotation_colors = ann_colors_geno,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)

# Timepoint-driven
pheatmap(dat_scaled_time,
         annotation_col    = ann_col,
         annotation_row    = row_ann_time,
         annotation_colors = ann_colors_time,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)

# Function to cluster within each module and return ordered gene names
order_within_clusters <- function(mat, row_ann) {
  clusters <- unique(row_ann$Cluster)
  ordered_genes <- lapply(clusters, function(cl) {
    genes <- rownames(row_ann)[row_ann$Cluster == cl]
    if (length(genes) < 2) return(genes)
    sub_mat <- mat[genes, ]
    dist_mat <- as.dist(1 - cor(t(sub_mat)))
    hc <- hclust(dist_mat, method = "ward.D2")
    genes[hc$order]
  })
  unlist(ordered_genes)
}

# Get ordered gene names
gene_order_geno <- order_within_clusters(dat_scaled_geno, row_ann_geno)
gene_order_time <- order_within_clusters(dat_scaled_time, row_ann_time)

# Reorder matrices
dat_scaled_geno_ord <- dat_scaled_geno[gene_order_geno, ]
dat_scaled_time_ord <- dat_scaled_time[gene_order_time, ]

# Reorder annotations to match
row_ann_geno_ord <- row_ann_geno[gene_order_geno, , drop = FALSE]
row_ann_time_ord <- row_ann_time[gene_order_time, , drop = FALSE]

# Plot with cluster_rows = FALSE to preserve order
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_ord,
         annotation_colors = ann_colors_geno,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)

pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_ord,
         annotation_colors = ann_colors_time,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)

# Define column order
col_order <- c("WT_M16", "WT_B16", "WT_M20", "WT_B20",
               "L1_M16", "L1_B16", "L1_M20", "L1_B20",
               "T23_M16", "T23_B16", "T23_M20", "T23_B20",
               "L14_M16", "L14_B16", "L14_M20", "L14_B20")

# Reorder matrices
dat_scaled_geno_ord <- dat_scaled_geno_ord[, col_order]
dat_scaled_time_ord <- dat_scaled_time_ord[, col_order]

# Reorder annotation
ann_col_ord <- ann_col[col_order, ]

ann_colors_geno$treatment <- c("Mock" = "#d9d9d9", "Botrytis" = "#4a0000")
ann_colors_time$treatment <- c("Mock" = "#d9d9d9", "Botrytis" = "#4a0000")

pdf("genotype_driven_heatmap.pdf", width = 8, height = 6)
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_ord,
         annotation_colors = ann_colors_geno,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (2, 11, 4)"
)
dev.off()

pdf("timepoint_driven_heatmap.pdf", width = 8, height = 6)
pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_ord,
         annotation_colors = ann_colors_time,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (10, 7, 3)"
)
dev.off()

modules %>% filter(cluster == 2)

# Update genotype clusters
genotype_clusters <- c(6, 11, 4)
genotype_names <- paste0("Cluster ", genotype_clusters,
                         " (n = ", as.vector(cl_sizes[genotype_clusters]), ")")

# Rebuild gene set and scaled matrix
geno_genes <- modules$Gene_ID[modules$cluster %in% genotype_clusters]
dat_scaled_geno <- t(scale(t(dat_log[rownames(dat_log) %in% geno_genes, ])))
dat_scaled_geno <- dat_scaled_geno[complete.cases(dat_scaled_geno), ]

# Rebuild row annotation
row_ann_geno <- modules %>%
  filter(cluster %in% genotype_clusters) %>%
  filter(Gene_ID %in% rownames(dat_scaled_geno)) %>%
  mutate(cluster = factor(cluster)) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

# Reorder within clusters
gene_order_geno <- order_within_clusters(dat_scaled_geno, row_ann_geno)
dat_scaled_geno_ord <- dat_scaled_geno[gene_order_geno, col_order]
row_ann_geno_ord <- row_ann_geno[gene_order_geno, , drop = FALSE]

# Update cluster colors
ann_colors_geno$Cluster <- c(
  "4"  = "#d9d9d9",
  "6"  = "#969696",
  "11" = "#525252"
)

# Plot
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_ord,
         annotation_colors = ann_colors_geno,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (4, 6, 11)"
)

cl_geno <- lapply(seq_along(genotype_clusters), function(x) {
  cl_num_lists_sum[[genotype_clusters[x]]] %>%
    mutate(cluster_name = factor(genotype_names[x], levels = genotype_names))
}) %>% list_rbind()

cl_geno <- cl_geno %>%
  mutate(genotype = factor(genotype, levels = c("WT", "L1", "T23", "L14")))

cl_time <- cl_time %>%
  mutate(genotype = factor(genotype, levels = c("WT", "L1", "T23", "L14")))

p_genotype <- ggplot(cl_geno, aes(x = timepoint, y = mean,
                                  colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Genotype-driven clusters")

p_genotype

cl_time <- cl_time %>%
  mutate(genotype = factor(genotype, levels = c("WT", "L1", "T23", "L14")))

p_timepoint <- ggplot(cl_time, aes(x = timepoint, y = mean,
                                   colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Timepoint-driven clusters")

p_timepoint

# Genotype-driven heatmap — main figure
pdf("genotype_driven_heatmap.pdf", width = 4, height = 6)
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_ord,
         annotation_colors = ann_colors_geno,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (4, 6, 11)"
)
dev.off()

# Genotype-driven line plot — main figure
pdf("genotype_driven_lineplots.pdf", width = 4, height = 6)
p_genotype
dev.off()

# Timepoint-driven heatmap — supplement
pdf("timepoint_driven_heatmap_supp.pdf", width = 4, height = 6)
pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_ord,
         annotation_colors = ann_colors_time,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (3, 7, 10)"
)
dev.off()

# Timepoint-driven line plot — supplement
pdf("timepoint_driven_lineplots_supp.pdf", width = 4, height = 6)
p_timepoint
dev.off()

###
###
###
###
# this is for the genes in common

# Triple intersections from bc file
shared_BT16 <- Reduce(intersect, list(
  bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T16"],
  bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T16"],
  bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T16"]
))

shared_BT20 <- Reduce(intersect, list(
  bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T20"],
  bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T20"],
  bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T20"]
))

# Union of both timepoints
shared_bc_all <- union(shared_BT16, shared_BT20)

# Check sizes
length(shared_BT16)
length(shared_BT20)
length(shared_bc_all)

# Subset and scale the matrix
dat_shared <- dat_log[rownames(dat_log) %in% shared_bc_all, ]
dat_scaled_shared <- t(scale(t(dat_shared)))
dat_scaled_shared <- dat_scaled_shared[complete.cases(dat_scaled_shared), ]

# Reorder columns
dat_scaled_shared <- dat_scaled_shared[, col_order]

# Heatmap
pheatmap(dat_scaled_shared,
         annotation_col    = ann_col_ord,
         annotation_colors = ann_colors_geno,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Shared Botrytis DEGs (n = 1677)"
)

# Run hierarchical clustering to get dendrogram
dist_shared <- as.dist(1 - cor(t(dat_scaled_shared)))
hc_shared <- hclust(dist_shared, method = "ward.D2")

# Plot dendrogram to help decide number of clusters
plot(hc_shared, labels = FALSE, main = "Shared DEGs dendrogram")

# Cut dendrogram into 8 clusters
clusters_shared <- cutree(hc_shared, k = 8)

# Check sizes
table(clusters_shared)

# Build row annotation
row_ann_shared <- data.frame(
  Cluster = factor(clusters_shared),
  row.names = names(clusters_shared)
)

# Order genes by cluster then within-cluster correlation
gene_order_shared <- order_within_clusters(dat_scaled_shared, row_ann_shared)
dat_scaled_shared_ord <- dat_scaled_shared[gene_order_shared, ]
row_ann_shared_ord <- row_ann_shared[gene_order_shared, , drop = FALSE]

# Cluster colors - 8 grays
ann_colors_shared <- ann_colors_geno
ann_colors_shared$Cluster <- c(
  "1" = "#f7f7f7",
  "2" = "#d9d9d9",
  "3" = "#bdbdbd",
  "4" = "#969696",
  "5" = "#737373",
  "6" = "#525252",
  "7" = "#252525",
  "8" = "#000000"
)

# Plot
pheatmap(dat_scaled_shared_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_shared_ord,
         annotation_colors = ann_colors_shared,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Shared Botrytis DEGs (n = 1677)"
)


### 

genotype_clusters <- c(11, 4, 6)
genotype_names <- paste0("Cluster ", genotype_clusters,
                         " (n = ", as.vector(cl_sizes[genotype_clusters]), ")")

cl_geno <- lapply(seq_along(genotype_clusters), function(x) {
  cl_num_lists_sum[[genotype_clusters[x]]] %>%
    mutate(cluster_name = factor(genotype_names[x], levels = genotype_names))
}) %>% list_rbind() %>%
  mutate(genotype = factor(genotype, levels = c("WT", "L1", "T23", "L14")))

p_genotype <- ggplot(cl_geno, aes(x = timepoint, y = mean,
                                  colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Genotype-driven clusters")

p_genotype

###
## GO analysis on each cluster:

library(clusterProfiler)
library(org.At.tair.db)

# Gene lists per cluster
cluster_gene_lists <- lapply(genotype_clusters, function(cl) {
  modules$Gene_ID[modules$cluster == cl]
})
names(cluster_gene_lists) <- paste0("Cluster_", genotype_clusters)

# Background gene set — all genes in your normalized counts matrix
background_genes <- rownames(dat_log)

# Run GO enrichment per cluster
go_results <- lapply(names(cluster_gene_lists), function(cl) {
  enrichGO(
    gene          = cluster_gene_lists[[cl]],
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",          # Biological Process — change to "MF" or "CC" if needed
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_results) <- names(cluster_gene_lists)

# Quick summary of how many terms per cluster
sapply(go_results, function(x) nrow(x@result[x@result$p.adjust < 0.05, ]))


## ALL clusters
# All 11 clusters
all_clusters <- 1:11
all_cluster_gene_lists <- lapply(all_clusters, function(cl) {
  modules$Gene_ID[modules$cluster == cl]
})
names(all_cluster_gene_lists) <- paste0("Cluster_", all_clusters)

# Run GO enrichment for all clusters
go_results_all <- lapply(names(all_cluster_gene_lists), function(cl) {
  enrichGO(
    gene          = all_cluster_gene_lists[[cl]],
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_results_all) <- names(all_cluster_gene_lists)

# Summary of significant terms per cluster
sapply(go_results_all, function(x) nrow(x@result[x@result$p.adjust < 0.05, ]))

install.packages("writexl")
library(writexl)

go_tables <- lapply(names(go_results_all), function(cl) {
  go_results_all[[cl]]@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
})
names(go_tables) <- names(go_results_all)

write_xlsx(go_tables, "GO_enrichment_all_clusters.xlsx")

# Check raw results before filtering
nrow(go_results_all[["Cluster_11"]]@result)
head(go_results_all[["Cluster_11"]]@result)

# Check how many genes matched in the database
length(all_cluster_gene_lists[["Cluster_11"]])
sum(all_cluster_gene_lists[["Cluster_11"]] %in% keys(org.At.tair.db, keytype = "TAIR"))

go_tables <- lapply(names(go_results_all), function(cl) {
  go_results_all[[cl]]@result %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
})
names(go_tables) <- names(go_results_all)

write_xlsx(go_tables, "GO_enrichment_all_clusters_unfiltered.xlsx")


#####

# T16 row annotations
row_anns <- bc %>%
  dplyr::filter(Timepoint == "T16") %>%
  dplyr::select(Gene_ID, DE_Direction, Genotype) %>%
  pivot_wider(names_from = Genotype, values_from = DE_Direction)

row_anns[is.na(row_anns)] <- "Not DE"
row_anns <- as.data.frame(row_anns)
rownames(row_anns) <- row_anns$Gene_ID
row_anns <- dplyr::select(row_anns, Gene_ID, T23, L14, L1, WT)

# T20 row annotations
row_anns_20 <- bc %>%
  dplyr::filter(Timepoint == "T20") %>%
  dplyr::select(Gene_ID, DE_Direction, Genotype, Timepoint) %>%
  transmute(Gene_ID, DE_Direction,
            Genotype = paste0(Genotype, "_", Timepoint)) %>%
  pivot_wider(names_from = Genotype, values_from = DE_Direction)

row_anns_20[is.na(row_anns_20)] <- "Not DE"
row_anns_20 <- as.data.frame(row_anns_20)
rownames(row_anns_20) <- row_anns_20$Gene_ID
row_anns_20 <- dplyr::select(row_anns_20, Gene_ID, T23_T20, L14_T20, L1_T20, WT_T20)

# Combine
all_rows <- full_join(row_anns, row_anns_20, by = "Gene_ID")
all_rows[is.na(all_rows)] <- "Not DE"
rownames(all_rows) <- all_rows$Gene_ID
all_rows <- all_rows %>%
  dplyr::select(T23_T20, T23, L14_T20, L14, L1_T20, L1, WT_T20, WT)

# Annotation colors from collaborator
ann_colours <- list(
  WT       = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  L1       = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  L14      = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  T23      = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  WT_T20   = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  L1_T20   = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  L14_T20  = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5"),
  T23_T20  = c("Upregulated" = "#1e81b0", "Downregulated" = "#e28743", "Not DE" = "#B2BEB5")
)

# Check
dim(all_rows)
head(all_rows)

# Subset all_rows to genes in each heatmap and merge with cluster annotations
row_ann_geno_full <- row_ann_geno_ord %>%
  merge(all_rows[rownames(all_rows) %in% rownames(dat_scaled_geno_ord), ],
        by = "row.names", all.x = TRUE) %>%
  tibble::column_to_rownames("Row.names")
row_ann_geno_full[is.na(row_ann_geno_full)] <- "Not DE"
row_ann_geno_full <- row_ann_geno_full[rownames(dat_scaled_geno_ord), ]

row_ann_time_full <- row_ann_time_ord %>%
  merge(all_rows[rownames(all_rows) %in% rownames(dat_scaled_time_ord), ],
        by = "row.names", all.x = TRUE) %>%
  tibble::column_to_rownames("Row.names")
row_ann_time_full[is.na(row_ann_time_full)] <- "Not DE"
row_ann_time_full <- row_ann_time_full[rownames(dat_scaled_time_ord), ]

# Add ann_colours to existing annotation color lists
ann_colors_geno_full <- ann_colors_geno
ann_colors_geno_full <- c(ann_colors_geno_full, ann_colours)

ann_colors_time_full <- ann_colors_time
ann_colors_time_full <- c(ann_colors_time_full, ann_colours)

# Genotype heatmap
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_full,
         annotation_colors = ann_colors_geno_full,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (11, 4, 6)"
)

# Timepoint heatmap
pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_full,
         annotation_colors = ann_colors_time_full,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (3, 7, 10)"
)

# Create single DEG status column
deg_any <- all_rows %>%
  mutate(DEG_status = case_when(
    if_any(everything(), ~ . == "Upregulated")   ~ "Upregulated",
    if_any(everything(), ~ . == "Downregulated") ~ "Downregulated",
    TRUE ~ "Not DE"
  )) %>%
  dplyr::select(DEG_status)

# Merge with cluster annotations
row_ann_geno_simple <- row_ann_geno_ord %>%
  merge(deg_any[rownames(deg_any) %in% rownames(dat_scaled_geno_ord), , drop = FALSE],
        by = "row.names", all.x = TRUE) %>%
  tibble::column_to_rownames("Row.names")
row_ann_geno_simple[is.na(row_ann_geno_simple)] <- "Not DE"
row_ann_geno_simple <- row_ann_geno_simple[rownames(dat_scaled_geno_ord), ]

row_ann_time_simple <- row_ann_time_ord %>%
  merge(deg_any[rownames(deg_any) %in% rownames(dat_scaled_time_ord), , drop = FALSE],
        by = "row.names", all.x = TRUE) %>%
  tibble::column_to_rownames("Row.names")
row_ann_time_simple[is.na(row_ann_time_simple)] <- "Not DE"
row_ann_time_simple <- row_ann_time_simple[rownames(dat_scaled_time_ord), ]

# Simplified annotation colors
ann_colors_geno_simple <- ann_colors_geno
ann_colors_geno_simple$DEG_status <- c(
  "Upregulated"   = "#1e81b0",
  "Downregulated" = "#e28743",
  "Not DE"        = "#B2BEB5"
)

ann_colors_time_simple <- ann_colors_time
ann_colors_time_simple$DEG_status <- c(
  "Upregulated"   = "#1e81b0",
  "Downregulated" = "#e28743",
  "Not DE"        = "#B2BEB5"
)

ann_colors_geno_simple$Cluster <- c(
  "11" = "#d9d9d9",
  "4"  = "#969696",
  "6"  = "#525252"
)

ann_colors_time_simple$Cluster <- c(
  "3"  = "#d9d9d9",
  "7"  = "#969696",
  "10" = "#525252"
)

# Genotype heatmap
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_simple,
         annotation_colors = ann_colors_geno_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (11, 4, 6)"
)

# Timepoint heatmap
pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_simple,
         annotation_colors = ann_colors_time_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (3, 7, 10)"
)


pdf("genotype_driven_heatmap.pdf", width = 4, height = 6)
pheatmap(dat_scaled_geno_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_geno_simple,
         annotation_colors = ann_colors_geno_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (11, 4, 6)"
)
dev.off()

pdf("timepoint_driven_heatmap_supp.pdf", width = 4, height = 6)
pheatmap(dat_scaled_time_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_time_simple,
         annotation_colors = ann_colors_time_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (3, 7, 10)"
)
dev.off()

#####
# Ok this is experimental but I am moving to cytoscape to look at network changes in response to botrytis. 
# Export DEG lists for Cytoscape import
# Shared Botrytis DEGs (center of Venn)
write.csv(data.frame(Gene_ID = shared_bc_all),
          "cytoscape_shared_botrytis_degs.csv", row.names = FALSE)

# Per genotype DEG lists from bc file
for (geno in c("WT", "L1", "L14", "T23")) {
  genes <- bc$Gene_ID[bc$Genotype == geno]
  write.csv(data.frame(Gene_ID = unique(genes)),
            paste0("cytoscape_", geno, "_botrytis_degs.csv"),
            row.names = FALSE)
}

# Module cluster assignments for node annotation
write.csv(modules, "cytoscape_module_assignments.csv", row.names = FALSE)

# Full DEG status table for node coloring
deg_node_attrs <- bc %>%
  dplyr::select(Gene_ID, Genotype, Timepoint, DE_Direction, log2FoldChange) %>%
  dplyr::mutate(contrast = paste(Genotype, Timepoint, sep = "_"))

write.csv(deg_node_attrs, "cytoscape_node_attributes.csv", row.names = FALSE)



#####
# Load the exported node table
cyto_nodes <- read.csv("DEGs_0.6--clusteredMCLif4 default  node.csv", check.names = FALSE)

# Check column names
colnames(cyto_nodes)

# Rename for convenience
cyto_nodes <- cyto_nodes %>%
  dplyr::rename(
    Gene_ID  = `query term`,
    cluster  = `__mclCluster`,
    symbol   = `display name`
  )

# Check cluster sizes — filter out singletons
cluster_sizes_mcl <- table(cyto_nodes$cluster)
cluster_sizes_mcl[cluster_sizes_mcl >= 5]

# How many usable clusters
sum(cluster_sizes_mcl >= 5)

# Get cluster IDs with n >= 5
valid_clusters <- names(cluster_sizes_mcl[cluster_sizes_mcl >= 5])

# Build gene lists per MCL cluster
mcl_gene_lists <- lapply(valid_clusters, function(cl) {
  cyto_nodes$Gene_ID[cyto_nodes$cluster == cl & !is.na(cyto_nodes$cluster)]
})
names(mcl_gene_lists) <- paste0("MCL_", valid_clusters)

# Background — all genes in your normalized counts
background_genes <- rownames(dat_log)

# Run GO enrichment per cluster
go_mcl <- lapply(names(mcl_gene_lists), function(cl) {
  enrichGO(
    gene          = mcl_gene_lists[[cl]],
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_mcl) <- names(mcl_gene_lists)

# Summary of significant terms per cluster
sapply(go_mcl, function(x) nrow(x@result[x@result$p.adjust < 0.05, ]))

# Build cluster stats separately
cluster_summary <- cyto_nodes %>%
  dplyr::filter(cluster %in% valid_clusters) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_genes = n(),
    n_up    = sum(DE_Direction == "Upregulated", na.rm = TRUE),
    n_down  = sum(DE_Direction == "Downregulated", na.rm = TRUE),
    .groups = "drop"
  )

go_term_counts <- data.frame(
  cluster    = as.integer(valid_clusters),
  n_go_terms = sapply(go_mcl, function(x) 
    nrow(x@result[x@result$p.adjust < 0.05, ]))
)

cluster_summary <- cluster_summary %>%
  left_join(go_term_counts, by = "cluster") %>%
  dplyr::arrange(as.numeric(cluster))

cluster_summary

# Export GO tables with cluster stats as first sheet
go_mcl_tables <- lapply(names(go_mcl), function(cl) {
  go_mcl[[cl]]@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
})
names(go_mcl_tables) <- names(go_mcl)

# Add cluster summary as first sheet
go_mcl_tables <- c(list(Summary = as.data.frame(cluster_summary)), go_mcl_tables)

write_xlsx(go_mcl_tables, "GO_enrichment_MCL_clusters.xlsx")

print(cluster_summary, n = 43)

#1. Semantic similarity clustering (best option) 
#Use simplify() from clusterProfiler which removes redundant GO terms based on semantic similarity, keeping only the most representative term per group:

go_mcl_simplified <- lapply(go_mcl, simplify, 
                            cutoff = 0.7, 
                            by = "p.adjust",
                            select_fun = min)

BiocManager::install("rrvgo")
library(rrvgo)

# First simplify each cluster's GO results using rrvgo
# We need the GO term IDs and scores for each cluster

go_mcl_reduced <- lapply(names(go_mcl), function(cl) {
  res <- go_mcl[[cl]]@result %>%
    dplyr::filter(p.adjust < 0.05)
  
  if (nrow(res) == 0) return(NULL)
  
  # Calculate semantic similarity matrix
  simMatrix <- tryCatch(
    calculateSimMatrix(
      res$ID,
      orgdb  = "org.At.tair.db",
      ont    = "BP",
      method = "Rel"
    ),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(
      calculateSimMatrix(
        res$ID,
        orgdb  = "org.At.tair.db",
        ont    = "BP",
        method = "Rel"
      )
    )
  )
  
  # Guard against NULL, empty, or too-small matrices
  if (is.null(simMatrix)) return(NULL)
  if (!is.matrix(simMatrix)) return(NULL)
  if (nrow(simMatrix) < 2) return(NULL)
  
  # Reduce to representative terms
  scores <- setNames(-log10(res$p.adjust), res$ID)
  scores <- scores[names(scores) %in% rownames(simMatrix)]
  
  if (length(scores) < 2) return(NULL)
  
  tryCatch({
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = 0.7,
                                    orgdb     = "org.At.tair.db")
    reducedTerms %>% mutate(cluster = cl)
  }, error = function(e) NULL)
})
names(go_mcl_reduced) <- names(go_mcl)

# Remove NULLs
go_mcl_reduced <- Filter(Negate(is.null), go_mcl_reduced)

sapply(go_mcl_reduced, function(x) length(unique(x$parentTerm)))

#### switch away from MCL clustersL
# Conserved — DEG in WT AND at least one mutant
wt_degs <- bc$Gene_ID[bc$Genotype == "WT"]
mutant_degs <- bc$Gene_ID[bc$Genotype != "WT"]

conserved <- intersect(wt_degs, mutant_degs)

# Diverged — DEG in at least one mutant but NOT in WT
diverged <- setdiff(mutant_degs, wt_degs)

length(conserved)
length(diverged)

# GO enrichment on conserved vs diverged
go_conserved <- enrichGO(
  gene          = unique(conserved),
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_diverged <- enrichGO(
  gene          = unique(diverged),
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(go_conserved@result[go_conserved@result$p.adjust < 0.05, ])
nrow(go_diverged@result[go_diverged@result$p.adjust < 0.05, ])

# rrvgo on conserved
simMatrix_conserved <- calculateSimMatrix(
  go_conserved@result$ID[go_conserved@result$p.adjust < 0.05],
  orgdb  = "org.At.tair.db",
  ont    = "BP",
  method = "Rel"
)

scores_conserved <- setNames(
  -log10(go_conserved@result$p.adjust[go_conserved@result$p.adjust < 0.05]),
  go_conserved@result$ID[go_conserved@result$p.adjust < 0.05]
)
scores_conserved <- scores_conserved[names(scores_conserved) %in% rownames(simMatrix_conserved)]

reduced_conserved <- reduceSimMatrix(simMatrix_conserved,
                                     scores_conserved,
                                     threshold = 0.7,
                                     orgdb     = "org.At.tair.db") %>%
  mutate(response = "Conserved")

# rrvgo on diverged
simMatrix_diverged <- calculateSimMatrix(
  go_diverged@result$ID[go_diverged@result$p.adjust < 0.05],
  orgdb  = "org.At.tair.db",
  ont    = "BP",
  method = "Rel"
)

scores_diverged <- setNames(
  -log10(go_diverged@result$p.adjust[go_diverged@result$p.adjust < 0.05]),
  go_diverged@result$ID[go_diverged@result$p.adjust < 0.05]
)
scores_diverged <- scores_diverged[names(scores_diverged) %in% rownames(simMatrix_diverged)]

reduced_diverged <- reduceSimMatrix(simMatrix_diverged,
                                    scores_diverged,
                                    threshold = 0.7,
                                    orgdb     = "org.At.tair.db") %>%
  mutate(response = "Diverged")

# Check representative terms
length(unique(reduced_conserved$parentTerm))
length(unique(reduced_diverged$parentTerm))

# Combine both sets
reduced_all <- bind_rows(reduced_conserved, reduced_diverged)

# Get p.adjust values for dot sizing from original GO results
conserved_padj <- go_conserved@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(ID, p.adjust, Count)

diverged_padj <- go_diverged@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(ID, p.adjust, Count)

all_padj <- bind_rows(
  conserved_padj %>% mutate(response = "Conserved"),
  diverged_padj  %>% mutate(response = "Diverged")
)

# Join to reduced terms - use parentTerm as the plotting label
plot_df <- reduced_all %>%
  left_join(all_padj, by = c("go" = "ID", "response")) %>%
  dplyr::group_by(parentTerm, response) %>%
  dplyr::summarise(
    mean_padj = mean(p.adjust, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(neg_log10_padj = -log10(mean_padj))

# Plot
ggplot(plot_df, aes(x = response, y = parentTerm,
                    size = total_count,
                    color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Conserved vs Diverged Botrytis response"
  ) +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 8),
    axis.text.x  = element_text(size = 10),
    panel.grid   = element_blank(),
    legend.position = "right"
  )

pdf("GO_conserved_vs_diverged_supplement.pdf", width = 8, height = 12)
ggplot(plot_df, aes(x = response, y = parentTerm,
                    size = total_count,
                    color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Conserved vs Diverged Botrytis response"
  ) +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    panel.grid   = element_blank(),
    legend.position = "right"
  )
dev.off()


## refine size
# Get top 10 per response category by significance
plot_df_top <- plot_df %>%
  dplyr::group_by(response) %>%
  dplyr::slice_min(mean_padj, n = 10) %>%
  dplyr::ungroup()

# Plot
ggplot(plot_df_top, aes(x = response, y = reorder(parentTerm, neg_log10_padj),
                        size = total_count,
                        color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Conserved vs Diverged Botrytis response"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 10),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
### changing axis

# WT Botrytis DEGs
wt_bc <- unique(bc$Gene_ID[bc$Genotype == "WT"])

# Per-mutant Botrytis DEGs
l1_bc  <- unique(bc$Gene_ID[bc$Genotype == "L1"])
l14_bc <- unique(bc$Gene_ID[bc$Genotype == "L14"])
t23_bc <- unique(bc$Gene_ID[bc$Genotype == "T23"])

# Conserved in WT — DEG in WT AND at least one mutant
conserved_wt <- intersect(wt_bc, union(l1_bc, union(l14_bc, t23_bc)))

# Mutant-specific — DEG in that mutant but NOT in WT
l1_specific  <- setdiff(l1_bc,  wt_bc)
l14_specific <- setdiff(l14_bc, wt_bc)
t23_specific <- setdiff(t23_bc, wt_bc)

# Check sizes
length(conserved_wt)
length(l1_specific)
length(l14_specific)
length(t23_specific)

gene_sets <- list(
  "Conserved (WT)" = conserved_wt,
  "L1 specific"    = l1_specific,
  "L14 specific"   = l14_specific,
  "m23 specific"   = t23_specific
)

go_by_genotype <- lapply(names(gene_sets), function(gs) {
  enrichGO(
    gene          = gene_sets[[gs]],
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
})
names(go_by_genotype) <- names(gene_sets)

# Check significant terms per set
sapply(go_by_genotype, function(x) 
  nrow(x@result[x@result$p.adjust < 0.05, ]))

go_reduced_genotype <- lapply(names(go_by_genotype), function(gs) {
  res <- go_by_genotype[[gs]]@result %>%
    dplyr::filter(p.adjust < 0.05)
  
  if (nrow(res) == 0) return(NULL)
  
  simMatrix <- tryCatch(
    calculateSimMatrix(
      res$ID,
      orgdb  = "org.At.tair.db",
      ont    = "BP",
      method = "Rel"
    ),
    error = function(e) NULL,
    warning = function(w) suppressWarnings(
      calculateSimMatrix(
        res$ID,
        orgdb  = "org.At.tair.db",
        ont    = "BP",
        method = "Rel"
      )
    )
  )
  
  if (is.null(simMatrix)) return(NULL)
  if (!is.matrix(simMatrix)) return(NULL)
  if (nrow(simMatrix) < 2) return(NULL)
  
  scores <- setNames(-log10(res$p.adjust), res$ID)
  scores <- scores[names(scores) %in% rownames(simMatrix)]
  
  if (length(scores) < 2) return(NULL)
  
  tryCatch({
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = 0.7,
                                    orgdb     = "org.At.tair.db") %>%
      mutate(response = gs)
    reducedTerms
  }, error = function(e) NULL)
})
names(go_reduced_genotype) <- names(go_by_genotype)
go_reduced_genotype <- Filter(Negate(is.null), go_reduced_genotype)

# Check representative terms per set
sapply(go_reduced_genotype, function(x) length(unique(x$parentTerm)))

# Get p.adjust and Count from original GO results
all_padj_genotype <- lapply(names(go_by_genotype), function(gs) {
  go_by_genotype[[gs]]@result %>%
    dplyr::filter(p.adjust < 0.05) %>%
    dplyr::select(ID, p.adjust, Count) %>%
    mutate(response = gs)
}) %>% list_rbind()

# Combine reduced terms
reduced_all_genotype <- list_rbind(go_reduced_genotype)

# Join and summarise per parent term
plot_df_genotype <- reduced_all_genotype %>%
  left_join(all_padj_genotype, by = c("go" = "ID", "response")) %>%
  dplyr::group_by(parentTerm, response) %>%
  dplyr::summarise(
    mean_padj   = mean(p.adjust, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(neg_log10_padj = -log10(mean_padj))

# Top 10 per category by significance
plot_df_top_genotype <- plot_df_genotype %>%
  dplyr::group_by(response) %>%
  dplyr::slice_min(mean_padj, n = 10) %>%
  dplyr::ungroup() %>%
  mutate(response = factor(response,
                           levels = c("Conserved (WT)", "L1 specific",
                                      "m23 specific", "L14 specific")))

# Plot
ggplot(plot_df_top_genotype,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Botrytis response by genotype"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

pdf("GO_botrytis_by_genotype.pdf", width = 5, height = 8)
ggplot(plot_df_top_genotype,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Botrytis response by genotype"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

# Print all parent terms for each category
print(plot_df_top_genotype %>% 
        dplyr::select(response, parentTerm, mean_padj, total_count) %>%
        dplyr::arrange(response, mean_padj), n = 40)


keep_terms <- c(
  "response to wounding",
  "immune system process",
  "glutathione metabolic process",
  "indolalkylamine biosynthetic process",
  "plant-type hypersensitive response",
  "response to oxidative stress",
  "response to fatty acid",
  "sulfur compound biosynthetic process",
  "S-glycoside biosynthetic process",
  "plastid transcription",
  "response to hypoxia",
  "protein phosphorylation",
  "leaf senescence",
  "immune response",
  "response to nitrogen compound",
  "starch metabolic process",
  "photosynthesis",
  "plastid organization",
  "response to light intensity",
  "pigment metabolic process",
  "chlorophyll metabolic process",
  "monosaccharide metabolic process",
  "DNA replication",
  "defense response to bacterium",
  "inorganic anion transport",
  "regulation of DNA endoreduplication",
  "response to starvation",
  "glutamine metabolic process"
)

plot_df_manual <- plot_df_genotype %>%
  dplyr::filter(parentTerm %in% keep_terms) %>%
  mutate(response = factor(response,
                           levels = c("Conserved (WT)", "L1 specific",
                                      "m23 specific", "L14 specific")))

ggplot(plot_df_manual,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

pdf("GO_botrytis_by_genotype_sub1.pdf", width = 4.5, height = 6.5)
ggplot(plot_df_manual,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

### going back and trying to see whats unique about L14 at T16. 
# L14 specific at T16 — DEG in L14 but not in any other genotype at T16
l14_T16 <- bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T16"]
l1_T16  <- bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T16"]
t23_T16 <- bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T16"]
wt_T16  <- bc$Gene_ID[bc$Genotype == "WT"  & bc$Timepoint == "T16"]

# Unique to L14 at T16
l14_T16_unique <- setdiff(l14_T16, 
                          union(l1_T16, union(t23_T16, wt_T16)))

length(l14_T16_unique)

go_l14_T16_unique <- enrichGO(
  gene          = l14_T16_unique,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(go_l14_T16_unique@result[go_l14_T16_unique@result$p.adjust < 0.05, ])

# View top terms
go_l14_T16_unique@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, GeneRatio, p.adjust, geneID) %>%
  head(20)

# Check what was in the across_bc_T16 list
length(across_bc_T16[["L14"]])

# L14 unique from the Venn — not in L1 or T23 (but may include WT)
l14_T16_venn <- setdiff(l14_T16, union(l1_T16, t23_T16))
length(l14_T16_venn)

# L14 unique excluding only mutants, keeping WT overlap
l14_T16_nomuts <- setdiff(l14_T16, union(l1_T16, t23_T16))
length(l14_T16_nomuts)


go_l14_T16_unique <- enrichGO(
  gene          = l14_T16_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(go_l14_T16_unique@result[go_l14_T16_unique@result$p.adjust < 0.05, ])

# View top terms
go_l14_T16_unique@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, GeneRatio, p.adjust, geneID) %>%
  head(20)

# Relax to raw pvalue
go_l14_T16_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(Description, GeneRatio, pvalue, p.adjust, geneID) %>%
  head(20)

# What are the L14 unique T16 genes?
l14_T16_venn

# L14 unique across both timepoints
l14_T20 <- bc$Gene_ID[bc$Genotype == "L14" & bc$Timepoint == "T20"]
l1_T20  <- bc$Gene_ID[bc$Genotype == "L1"  & bc$Timepoint == "T20"]
t23_T20 <- bc$Gene_ID[bc$Genotype == "T23" & bc$Timepoint == "T20"]

l14_T20_venn <- setdiff(l14_T20, union(l1_T20, t23_T20))

# Combined across both timepoints
l14_both_unique <- union(l14_T16_venn, l14_T20_venn)
length(l14_both_unique)

go_l14_both_unique <- enrichGO(
  gene          = l14_both_unique,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(go_l14_both_unique@result[go_l14_both_unique@result$p.adjust < 0.05, ])

# View top terms
go_l14_both_unique@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::select(Description, GeneRatio, p.adjust, geneID) %>%
  head(20)

go_l14_both_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(Description, GeneRatio, pvalue, p.adjust, geneID) %>%
  head(20)

go_l14_T16_table <- go_l14_T16_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(Description, GeneRatio, pvalue, p.adjust, geneID)

go_l14_both_table <- go_l14_both_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(Description, GeneRatio, pvalue, p.adjust, geneID)

write_xlsx(
  list(
    "L14_unique_T16"   = go_l14_T16_table,
    "L14_unique_T16_T20" = go_l14_both_table
  ),
  "GO_L14_unique_DEGs.xlsx"
)

# Helper function to safely map symbols to TAIR IDs
safe_map_tair <- function(geneID_string) {
  symbols <- strsplit(geneID_string, "/")[[1]]
  symbols <- symbols[symbols != "NA" & !is.na(symbols) & symbols != ""]
  if (length(symbols) == 0) return(NA_character_)
  tair_ids <- tryCatch({
    mapIds(org.At.tair.db,
           keys    = symbols,
           column  = "TAIR",
           keytype = "SYMBOL",
           multiVals = "first")
  }, error = function(e) {
    setNames(rep(NA_character_, length(symbols)), symbols)
  })
  paste(tair_ids, collapse = "/")
}

# Apply to each GO result
get_tair_ids <- function(go_result) {
  go_result@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(ID, Description, GeneRatio, pvalue, p.adjust, geneID) %>%
    mutate(TAIR_IDs = sapply(geneID, safe_map_tair))
}

go_l14_T16_table  <- get_tair_ids(go_l14_T16_unique)
go_l14_both_table <- get_tair_ids(go_l14_both_unique)

write_xlsx(
  list(
    "L14_unique_T16"     = go_l14_T16_table,
    "L14_unique_T16_T20" = go_l14_both_table
  ),
  "GO_L14_unique_DEGs_withTAIR.xlsx"
)

# L14 unique at T20
l14_T20_venn <- setdiff(l14_T20, union(l1_T20, t23_T20))
length(l14_T20_venn)

go_l14_T20_unique <- enrichGO(
  gene          = l14_T20_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

nrow(go_l14_T20_unique@result[go_l14_T20_unique@result$p.adjust < 0.05, ])

# The one significant term
go_l14_T20_unique@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(Description, GeneRatio, p.adjust, geneID)

# And relaxed for exploration
go_l14_T20_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(Description, GeneRatio, pvalue, p.adjust, geneID) %>%
  head(20)

go_l14_T16_table <- get_tair_ids(go_l14_T16_unique)
go_l14_T20_table <- get_tair_ids(go_l14_T20_unique)

write_xlsx(
  list(
    "L14_unique_T16" = go_l14_T16_table,
    "L14_unique_T20" = go_l14_T20_table
  ),
  "GO_L14_unique_DEGs_withTAIR.xlsx"
)
### extend to L1 and m23 next
# L1 unique at T16 and T20
l1_T16_venn <- setdiff(l1_T16, union(l14_T16, t23_T16))
l1_T20_venn <- setdiff(l1_T20, union(l14_T20, t23_T20))

# m23 unique at T16 and T20
t23_T16_venn <- setdiff(t23_T16, union(l1_T16, l14_T16))
t23_T20_venn <- setdiff(t23_T20, union(l1_T20, l14_T20))

# Check sizes
length(l1_T16_venn)
length(l1_T20_venn)
length(t23_T16_venn)
length(t23_T20_venn)

go_l1_T16 <- enrichGO(
  gene          = l1_T16_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_l1_T20 <- enrichGO(
  gene          = l1_T20_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_t23_T16 <- enrichGO(
  gene          = t23_T16_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_t23_T20 <- enrichGO(
  gene          = t23_T20_venn,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Check significant terms
sapply(list(
  L1_T16  = go_l1_T16,
  L1_T20  = go_l1_T20,
  m23_T16 = go_t23_T16,
  m23_T20 = go_t23_T20
), function(x) nrow(x@result[x@result$p.adjust < 0.05, ]))

go_l1_T16_table  <- get_tair_ids(go_l1_T16)
go_l1_T20_table  <- get_tair_ids(go_l1_T20)
go_t23_T16_table <- go_t23_T16@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(ID, Description, GeneRatio, pvalue, p.adjust, geneID) %>%
  mutate(TAIR_IDs = sapply(geneID, safe_map_tair))
go_t23_T20_table <- get_tair_ids(go_t23_T20)

write_xlsx(
  list(
    "L1_T16"  = go_l1_T16_table,
    "L1_T20"  = go_l1_T20_table,
    "m23_T16" = go_t23_T16_table,
    "m23_T20" = go_t23_T20_table
  ),
  "GO_L1_m23_unique_DEGs_withTAIR.xlsx"
)

# Run rrvgo on all six sets using relaxed pvalue
go_list_unique <- list(
  "L1 T16"  = go_l1_T16,
  "L1 T20"  = go_l1_T20,
  "m23 T16" = go_t23_T16,
  "m23 T20" = go_t23_T20,
  "L14 T16" = go_l14_T16_unique,
  "L14 T20" = go_l14_T20_unique
)

# Reduce with rrvgo using relaxed pvalue
go_reduced_unique <- lapply(names(go_list_unique), function(gs) {
  res <- go_list_unique[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05)
  
  if (nrow(res) == 0) return(NULL)
  
  simMatrix <- tryCatch(
    calculateSimMatrix(
      res$ID,
      orgdb  = "org.At.tair.db",
      ont    = "BP",
      method = "Rel"
    ),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      calculateSimMatrix(
        res$ID,
        orgdb  = "org.At.tair.db",
        ont    = "BP",
        method = "Rel"
      )
    )
  )
  
  if (is.null(simMatrix)) return(NULL)
  if (!is.matrix(simMatrix)) return(NULL)
  if (nrow(simMatrix) < 2) return(NULL)
  
  scores <- setNames(-log10(res$pvalue), res$ID)
  scores <- scores[names(scores) %in% rownames(simMatrix)]
  if (length(scores) < 2) return(NULL)
  
  tryCatch({
    reduceSimMatrix(simMatrix,
                    scores,
                    threshold = 0.7,
                    orgdb     = "org.At.tair.db") %>%
      mutate(response = gs)
  }, error = function(e) NULL)
})
names(go_reduced_unique) <- names(go_list_unique)
go_reduced_unique <- Filter(Negate(is.null), go_reduced_unique)

# Check representative terms per set
sapply(go_reduced_unique, function(x) length(unique(x$parentTerm)))

# Get pvalue and Count from original GO results
all_padj_unique <- lapply(names(go_list_unique), function(gs) {
  go_list_unique[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::select(ID, pvalue, Count) %>%
    mutate(response = gs)
}) %>% list_rbind()

# Combine reduced terms
reduced_all_unique <- list_rbind(go_reduced_unique)

# Join and summarise per parent term
plot_df_unique <- reduced_all_unique %>%
  left_join(all_padj_unique, by = c("go" = "ID", "response")) %>%
  dplyr::group_by(parentTerm, response) %>%
  dplyr::summarise(
    mean_pvalue = mean(pvalue, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    neg_log10_pvalue = -log10(mean_pvalue),
    response = factor(response, levels = c(
      "L1 T16", "L1 T20",
      "m23 T16", "m23 T20",
      "L14 T16", "L14 T20"
    ))
  )

# Top 10 per category
plot_df_unique_top <- plot_df_unique %>%
  dplyr::group_by(response) %>%
  dplyr::slice_min(mean_pvalue, n = 10) %>%
  dplyr::ungroup()

# Plot
ggplot(plot_df_unique_top,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_pvalue),
           size = total_count,
           color = neg_log10_pvalue)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.value)",
    size  = "Gene count"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

# Split by genotype and plot separately
genotype_levels <- list(
  "L1"  = c("L1 T16", "L1 T20"),
  "m23" = c("m23 T16", "m23 T20"),
  "L14" = c("L14 T16", "L14 T20")
)

plots_by_genotype <- lapply(names(genotype_levels), function(geno) {
  df <- plot_df_unique_top %>%
    dplyr::filter(response %in% genotype_levels[[geno]]) %>%
    mutate(response = factor(response, levels = genotype_levels[[geno]]))
  
  ggplot(df,
         aes(x = response,
             y = reorder(parentTerm, neg_log10_pvalue),
             size = total_count,
             color = neg_log10_pvalue)) +
    geom_point() +
    scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x     = NULL,
      y     = NULL,
      color = "-log10(p.value)",
      size  = "Gene count",
      title = geno
    ) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 8),
      axis.text.x     = element_text(size = 10),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
})
names(plots_by_genotype) <- names(genotype_levels)

plots_by_genotype[["L1"]]
plots_by_genotype[["m23"]]
plots_by_genotype[["L14"]]

pdf("GO_L1_unique_dotplot.pdf", width = 6, height = 8)
plots_by_genotype[["L1"]]
dev.off()

pdf("GO_m23_unique_dotplot.pdf", width = 6, height = 8)
plots_by_genotype[["m23"]]
dev.off()

pdf("GO_L14_unique_dotplot.pdf", width = 6, height = 8)
plots_by_genotype[["L14"]]
dev.off()

## sub analysis on the 51 L14 module
# Subset normalized counts to L14 T16 unique genes
dat_l14_T16 <- dat_log[rownames(dat_log) %in% l14_T16_venn, ]
dat_scaled_l14_T16 <- t(scale(t(dat_l14_T16)))
dat_scaled_l14_T16 <- dat_scaled_l14_T16[complete.cases(dat_scaled_l14_T16), ]

# Reorder columns
dat_scaled_l14_T16 <- dat_scaled_l14_T16[, col_order]

# Check dimensions
dim(dat_scaled_l14_T16)
pheatmap(dat_scaled_l14_T16,
         annotation_col    = ann_col_ord,
         annotation_colors = ann_colors_geno_simple,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 6,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "L14 unique T16 Botrytis DEGs (n = 51)"
)

# Plot and save the clustering object
p_l14_T16 <- pheatmap(dat_scaled_l14_T16,
                      annotation_col    = ann_col_ord,
                      annotation_colors = ann_colors_geno_simple,
                      clustering_distance_rows = "correlation",
                      clustering_method        = "ward.D2",
                      cluster_rows      = TRUE,
                      cluster_cols      = FALSE,
                      show_rownames     = TRUE,
                      show_colnames     = TRUE,
                      fontsize_row      = 6,
                      color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
                      main              = "L14 unique T16 Botrytis DEGs (n = 51)"
)

plot(p_l14_T16$tree_row)

# Plot dendrogram to decide number of clusters
plot(p_l14_T16$tree_row, main = "L14 T16 unique genes dendrogram")

dist_l14_T16 <- as.dist(1 - cor(t(dat_scaled_l14_T16)))
hc_l14_T16 <- hclust(dist_l14_T16, method = "ward.D2")
plot(hc_l14_T16, labels = rownames(dat_scaled_l14_T16), 
     main = "L14 T16 unique genes", cex = 0.6)

# Cut into 4 clusters
clusters_l14_T16 <- cutree(hc_l14_T16, k = 4)
table(clusters_l14_T16)

# Build row annotation
row_ann_l14_T16 <- data.frame(
  Cluster = factor(clusters_l14_T16),
  row.names = names(clusters_l14_T16)
)

# Order within clusters
gene_order_l14_T16 <- order_within_clusters(dat_scaled_l14_T16, row_ann_l14_T16)
dat_scaled_l14_T16_ord <- dat_scaled_l14_T16[gene_order_l14_T16, ]
row_ann_l14_T16_ord <- row_ann_l14_T16[gene_order_l14_T16, , drop = FALSE]

# Cluster colors
ann_colors_l14_T16 <- ann_colors_geno_simple
ann_colors_l14_T16$Cluster <- c(
  "1" = "#d9d9d9",
  "2" = "#969696",
  "3" = "#525252",
  "4" = "#252525"
)

# Plot
pheatmap(dat_scaled_l14_T16_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_l14_T16_ord,
         annotation_colors = ann_colors_l14_T16,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 6,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "L14 unique T16 Botrytis DEGs (n = 51)"
)

pdf("L14_T16_unique_heatmap.pdf", width = 8, height = 8)
pheatmap(dat_scaled_l14_T16_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_l14_T16_ord,
         annotation_colors = ann_colors_l14_T16,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 6,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "L14 unique T16 Botrytis DEGs (n = 51)"
)
dev.off()

# Join cluster assignments to normalized counts
l14_T16_longer <- dat_log[rownames(dat_log) %in% l14_T16_venn, ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  left_join(data.frame(Gene_ID = names(clusters_l14_T16),
                       cluster = clusters_l14_T16), by = "Gene_ID") %>%
  pivot_longer(cols = -c(Gene_ID, cluster),
               names_to = c("genotype", "condition"),
               names_sep = "_",
               values_to = "count") %>%
  separate(col = condition, into = c("treatment", "timepoint"), sep = 1) %>%
  mutate(
    timepoint = as.numeric(timepoint),
    genotype  = factor(genotype, levels = c("WT", "L1", "T23", "L14")),
    treatment = case_when(
      treatment == "B" ~ "Botrytis",
      treatment == "M" ~ "Mock"
    )
  )

# Scale per gene
l14_T16_scaled <- l14_T16_longer %>%
  group_by(Gene_ID) %>%
  mutate(scaled_expression = (count - mean(count)) / sd(count)) %>%
  ungroup()

# Summarise per cluster/genotype/treatment/timepoint
l14_T16_sum <- l14_T16_scaled %>%
  group_by(cluster, genotype, treatment, timepoint) %>%
  summarise(
    mean = mean(scaled_expression, na.rm = TRUE),
    n    = n(),
    sd   = sd(scaled_expression, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci   = 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(cluster_name = paste0("Cluster ", cluster,
                               " (n = ", table(clusters_l14_T16)[cluster], ")"),
         cluster_name = factor(cluster_name))

# Plot
ggplot(l14_T16_sum, aes(x = timepoint, y = mean,
                        colour = genotype, fill = genotype,
                        group = interaction(genotype, treatment))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci),
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(cluster_name ~ treatment, switch = "y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype",
       title = "L14 unique T16 Botrytis DEGs")

pdf("L14_T16_unique_lineplots.pdf", width = 4, height = 8)
ggplot(l14_T16_sum, aes(x = timepoint, y = mean,
                        colour = genotype, fill = genotype,
                        group = interaction(genotype, treatment))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci),
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(cluster_name ~ treatment, switch = "y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype",
       title = "L14 unique T16 Botrytis DEGs")
dev.off()


# Get L14 T16 reduced terms
l14_T16_plot_df <- go_reduced_unique[["L14 T16"]] %>%
  left_join(
    go_l14_T16_unique@result %>%
      dplyr::select(ID, pvalue, Count),
    by = c("go" = "ID")
  ) %>%
  dplyr::group_by(parentTerm) %>%
  dplyr::summarise(
    mean_pvalue = mean(pvalue, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(neg_log10_pvalue = -log10(mean_pvalue)) %>%
  dplyr::arrange(neg_log10_pvalue)

ggplot(l14_T16_plot_df,
       aes(x = reorder(parentTerm, neg_log10_pvalue),
           y = neg_log10_pvalue,
           fill = neg_log10_pvalue,
           size = total_count)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  coord_flip() +
  labs(
    x    = NULL,
    y    = "-log10(p.value)",
    fill = "-log10(p.value)",
    size = "Gene count",
    title = "L14 unique T16 Botrytis DEGs — GO enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggplot(l14_T16_plot_df,
       aes(x = neg_log10_pvalue,
           y = reorder(parentTerm, neg_log10_pvalue),
           fill = neg_log10_pvalue,
           size = total_count)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = "-log10(p.value)",
    y    = NULL,
    fill = "-log10(p.value)",
    size = "Gene count",
    title = "L14 unique T16 Botrytis DEGs — GO enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

# Check the genes in the plot vs the 51 gene list
cat("Genes in l14_T16_venn:\n")
l14_T16_venn

cat("\nGO results source — genes in any enriched term:\n")
go_l14_T16_unique@result %>%
  dplyr::filter(pvalue < 0.05) %>%
  dplyr::pull(geneID) %>%
  paste(collapse = "/") %>%
  strsplit("/") %>%
  unlist() %>%
  unique()

pdf("L14_T16_unique_GO_dotplot.pdf", width = 8, height = 3)
ggplot(l14_T16_plot_df,
       aes(x = neg_log10_pvalue,
           y = reorder(parentTerm, neg_log10_pvalue),
           fill = neg_log10_pvalue,
           size = total_count)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = "-log10(p.value)",
    y    = NULL,
    fill = "-log10(p.value)",
    size = "Gene count",
    title = "L14 unique T16 Botrytis DEGs - GO enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

write.csv(data.frame(Gene_ID = l14_T16_venn),
          "L14_specific_botrytis_T16_degs.csv", row.names = FALSE)

genes_of_interest <- c("AT1G72890", "AT3G22275", "AT3G55150", "AT4G34530")

# Check valid keytypes
keytypes(org.At.tair.db)

# Then try with explicit argument name
AnnotationDbi::mapIds(org.At.tair.db,
                      keys      = genes_of_interest,
                      column    = "SYMBOL",
                      keytype   = "TAIR",
                      multiVals = "first")

## heat map subset
# Subset to just these 4 genes
genes_4 <- c("AT1G72890", "AT3G22275", "AT3G55150", "AT4G34530")

dat_4 <- dat_log[rownames(dat_log) %in% genes_4, ]
dat_scaled_4 <- t(scale(t(dat_4)))

# Reorder columns
dat_scaled_4 <- dat_scaled_4[, col_order]

# Row annotation with gene symbols
row_ann_4 <- data.frame(
  Symbol = c("TN6", "JAZ13", "ATEXO70H1", "CIB1"),
  row.names = genes_4
)

pheatmap(dat_scaled_4,
         annotation_col    = ann_col_ord,
         annotation_colors = ann_colors_geno_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = TRUE,
         show_colnames     = TRUE,
         fontsize_row      = 10,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "L14 T16 unique Botrytis DEGs — selected genes"
)

# Pivot to long format
dat_4_longer <- dat_log[rownames(dat_log) %in% genes_4, ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  mutate(Symbol = case_when(
    Gene_ID == "AT1G72890" ~ "TN6",
    Gene_ID == "AT3G22275" ~ "JAZ13",
    Gene_ID == "AT3G55150" ~ "ATEXO70H1",
    Gene_ID == "AT4G34530" ~ "CIB1"
  )) %>%
  pivot_longer(cols = -c(Gene_ID, Symbol),
               names_to = c("genotype", "condition"),
               names_sep = "_",
               values_to = "count") %>%
  separate(col = condition, into = c("treatment", "timepoint"), sep = 1) %>%
  mutate(
    timepoint = as.numeric(timepoint),
    genotype  = factor(genotype, levels = c("WT", "L1", "T23", "L14")),
    treatment = case_when(
      treatment == "B" ~ "Botrytis",
      treatment == "M" ~ "Mock"
    )
  )

# Plot — raw expression, one line per genotype
ggplot(dat_4_longer, aes(x = timepoint, y = count,
                         colour = genotype,
                         group = interaction(genotype, treatment),
                         linetype = treatment)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  scale_linetype_manual(values = c("Mock" = "dashed", "Botrytis" = "solid")) +
  facet_wrap(~ Symbol, nrow = 2, scales = "free_y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Normalized counts") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", linetype = "Treatment")

ggplot(dat_4_longer, aes(x = timepoint, y = count,
                         colour = genotype,
                         group = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(Symbol ~ treatment, scales = "free_y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Normalized counts") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype")

pdf("JAZ13_TN6_ATEXO70H1_CIB1_lineplots.pdf", width = 4, height = 3)
ggplot(dat_4_longer, aes(x = timepoint, y = count,
                         colour = genotype,
                         group = genotype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(treatment ~ Symbol, scales = "free_y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Normalized counts") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype")
dev.off()

#### now lets do this again with m23
length(t23_T20_venn)

# Subset normalized counts to m23 T20 unique genes
dat_t23_T20 <- dat_log[rownames(dat_log) %in% t23_T20_venn, ]
dat_scaled_t23_T20 <- t(scale(t(dat_t23_T20)))
dat_scaled_t23_T20 <- dat_scaled_t23_T20[complete.cases(dat_scaled_t23_T20), ]

# Reorder columns
dat_scaled_t23_T20 <- dat_scaled_t23_T20[, col_order]

# Check dimensions
dim(dat_scaled_t23_T20)

# Cluster
dist_t23_T20 <- as.dist(1 - cor(t(dat_scaled_t23_T20)))
hc_t23_T20 <- hclust(dist_t23_T20, method = "ward.D2")

# Plot heatmap
pheatmap(dat_scaled_t23_T20,
         annotation_col    = ann_col_ord,
         annotation_colors = ann_colors_geno_simple,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "m23 unique T20 Botrytis DEGs (n = 1119)"
)

# Cut into 4 clusters
clusters_t23_T20 <- cutree(hc_t23_T20, k = 4)
table(clusters_t23_T20)

# Build row annotation
row_ann_t23_T20 <- data.frame(
  Cluster = factor(clusters_t23_T20),
  row.names = names(clusters_t23_T20)
)

# Order within clusters
gene_order_t23_T20 <- order_within_clusters(dat_scaled_t23_T20, row_ann_t23_T20)
dat_scaled_t23_T20_ord <- dat_scaled_t23_T20[gene_order_t23_T20, ]
row_ann_t23_T20_ord <- row_ann_t23_T20[gene_order_t23_T20, , drop = FALSE]

# Cluster colors
ann_colors_t23_T20 <- ann_colors_geno_simple
ann_colors_t23_T20$Cluster <- c(
  "1" = "#d9d9d9",
  "2" = "#969696",
  "3" = "#525252",
  "4" = "#252525"
)

# Plot
pheatmap(dat_scaled_t23_T20_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_t23_T20_ord,
         annotation_colors = ann_colors_t23_T20,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "m23 unique T20 Botrytis DEGs (n = 1119)"
)

pdf("m23_T20_unique_heatmap.pdf", width = 4, height = 8)
pheatmap(dat_scaled_t23_T20_ord,
         annotation_col    = ann_col_ord,
         annotation_row    = row_ann_t23_T20_ord,
         annotation_colors = ann_colors_t23_T20,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "m23 unique T20 Botrytis DEGs (n = 1119)"
)
dev.off()

# Join cluster assignments to normalized counts
t23_T20_longer <- dat_log[rownames(dat_log) %in% t23_T20_venn, ] %>%
  as.data.frame() %>%
  rownames_to_column("Gene_ID") %>%
  left_join(data.frame(Gene_ID = names(clusters_t23_T20),
                       cluster = clusters_t23_T20), by = "Gene_ID") %>%
  pivot_longer(cols = -c(Gene_ID, cluster),
               names_to = c("genotype", "condition"),
               names_sep = "_",
               values_to = "count") %>%
  separate(col = condition, into = c("treatment", "timepoint"), sep = 1) %>%
  mutate(
    timepoint = as.numeric(timepoint),
    genotype  = factor(genotype, levels = c("WT", "L1", "T23", "L14")),
    treatment = case_when(
      treatment == "B" ~ "Botrytis",
      treatment == "M" ~ "Mock"
    )
  )

# Scale per gene
t23_T20_scaled <- t23_T20_longer %>%
  group_by(Gene_ID) %>%
  mutate(scaled_expression = (count - mean(count)) / sd(count)) %>%
  ungroup()

# Summarise per cluster/genotype/treatment/timepoint
t23_T20_sum <- t23_T20_scaled %>%
  group_by(cluster, genotype, treatment, timepoint) %>%
  summarise(
    mean = mean(scaled_expression, na.rm = TRUE),
    n    = n(),
    sd   = sd(scaled_expression, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci   = 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(cluster_name = paste0("Cluster ", cluster,
                               " (n = ", table(clusters_t23_T20)[cluster], ")"),
         cluster_name = factor(cluster_name))

# Plot
ggplot(t23_T20_sum, aes(x = timepoint, y = mean,
                        colour = genotype, fill = genotype,
                        group = interaction(genotype, treatment))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci),
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(cluster_name ~ treatment, switch = "y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype",
       title = "m23 unique T20 Botrytis DEGs")

pdf("m23_T20_unique_lineplots.pdf", width = 4, height = 8)
ggplot(t23_T20_sum, aes(x = timepoint, y = mean,
                        colour = genotype, fill = genotype,
                        group = interaction(genotype, treatment))) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci),
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  scale_fill_manual(values = c(
    "WT"  = "#999999",
    "L1"  = "#197ad4",
    "T23" = "#1bd218",
    "L14" = "#8f00ff"
  )) +
  facet_grid(cluster_name ~ treatment, switch = "y") +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement  = "outside",
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype",
       title = "m23 unique T20 Botrytis DEGs")
dev.off()

# rrvgo reduction
simMatrix_t23_T20 <- calculateSimMatrix(
  go_t23_T20@result$ID[go_t23_T20@result$p.adjust < 0.05],
  orgdb  = "org.At.tair.db",
  ont    = "BP",
  method = "Rel"
)

scores_t23_T20 <- setNames(
  -log10(go_t23_T20@result$p.adjust[go_t23_T20@result$p.adjust < 0.05]),
  go_t23_T20@result$ID[go_t23_T20@result$p.adjust < 0.05]
)
scores_t23_T20 <- scores_t23_T20[names(scores_t23_T20) %in% rownames(simMatrix_t23_T20)]

reduced_t23_T20 <- reduceSimMatrix(simMatrix_t23_T20,
                                   scores_t23_T20,
                                   threshold = 0.7,
                                   orgdb     = "org.At.tair.db")

# Check number of representative terms
length(unique(reduced_t23_T20$parentTerm))

# Get p.adjust and Count from original GO results
t23_T20_padj <- go_t23_T20@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(ID, p.adjust, Count)

# Join and summarise per parent term
plot_df_t23_T20 <- reduced_t23_T20 %>%
  left_join(t23_T20_padj, by = c("go" = "ID")) %>%
  dplyr::group_by(parentTerm) %>%
  dplyr::summarise(
    mean_padj   = mean(p.adjust, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(neg_log10_padj = -log10(mean_padj))

# Plot
ggplot(plot_df_t23_T20,
       aes(x = neg_log10_padj,
           y = reorder(parentTerm, neg_log10_padj),
           fill = neg_log10_padj,
           size = total_count)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = "-log10(p.adjust)",
    y    = NULL,
    fill = "-log10(p.adjust)",
    size = "Gene count",
    title = "m23 unique T20 Botrytis DEGs - GO enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

pdf("m23_T20_unique_GO_dotplot.pdf", width = 5.2, height = 6)
ggplot(plot_df_t23_T20,
       aes(x = neg_log10_padj,
           y = reorder(parentTerm, neg_log10_padj),
           fill = neg_log10_padj,
           size = total_count)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = "-log10(p.adjust)",
    y    = NULL,
    fill = "-log10(p.adjust)",
    size = "Gene count",
    title = "m23 unique T20 Botrytis DEGs - GO enrichment"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 8),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

# Get full DEG info for m23 T20 unique genes from bc file
m23_T20_degs <- bc %>%
  dplyr::filter(Gene_ID %in% t23_T20_venn,
                Genotype == "T23",
                Timepoint == "T20") %>%
  dplyr::arrange(padj)

# Export
write_xlsx(
  list("m23_unique_T20_Botrytis_DEGs" = as.data.frame(m23_T20_degs)),
  "m23_unique_T20_Botrytis_DEGs.xlsx"
)

# Quick summary
nrow(m23_T20_degs)
table(m23_T20_degs$DE_Direction)

# Split by direction
m23_T20_up   <- m23_T20_degs$Gene_ID[m23_T20_degs$DE_Direction == "Upregulated"]
m23_T20_down <- m23_T20_degs$Gene_ID[m23_T20_degs$DE_Direction == "Downregulated"]

# GO enrichment
go_m23_T20_up <- enrichGO(
  gene          = m23_T20_up,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

go_m23_T20_down <- enrichGO(
  gene          = m23_T20_down,
  universe      = background_genes,
  OrgDb         = org.At.tair.db,
  keyType       = "TAIR",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Check significant terms
sapply(list(Up = go_m23_T20_up, Down = go_m23_T20_down),
       function(x) nrow(x@result[x@result$p.adjust < 0.05, ]))

# rrvgo on upregulated
simMatrix_m23_up <- calculateSimMatrix(
  go_m23_T20_up@result$ID[go_m23_T20_up@result$p.adjust < 0.05],
  orgdb  = "org.At.tair.db",
  ont    = "BP",
  method = "Rel"
)

scores_m23_up <- setNames(
  -log10(go_m23_T20_up@result$p.adjust[go_m23_T20_up@result$p.adjust < 0.05]),
  go_m23_T20_up@result$ID[go_m23_T20_up@result$p.adjust < 0.05]
)
scores_m23_up <- scores_m23_up[names(scores_m23_up) %in% rownames(simMatrix_m23_up)]

reduced_m23_up <- reduceSimMatrix(simMatrix_m23_up,
                                  scores_m23_up,
                                  threshold = 0.7,
                                  orgdb     = "org.At.tair.db") %>%
  mutate(direction = "Upregulated")

# rrvgo on downregulated
simMatrix_m23_down <- calculateSimMatrix(
  go_m23_T20_down@result$ID[go_m23_T20_down@result$p.adjust < 0.05],
  orgdb  = "org.At.tair.db",
  ont    = "BP",
  method = "Rel"
)

scores_m23_down <- setNames(
  -log10(go_m23_T20_down@result$p.adjust[go_m23_T20_down@result$p.adjust < 0.05]),
  go_m23_T20_down@result$ID[go_m23_T20_down@result$p.adjust < 0.05]
)
scores_m23_down <- scores_m23_down[names(scores_m23_down) %in% rownames(simMatrix_m23_down)]

reduced_m23_down <- reduceSimMatrix(simMatrix_m23_down,
                                    scores_m23_down,
                                    threshold = 0.7,
                                    orgdb     = "org.At.tair.db") %>%
  mutate(direction = "Downregulated")

# Check representative terms
length(unique(reduced_m23_up$parentTerm))
length(unique(reduced_m23_down$parentTerm))

# Get p.adjust and Count from original GO results
m23_up_padj <- go_m23_T20_up@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(ID, p.adjust, Count) %>%
  mutate(direction = "Upregulated")

m23_down_padj <- go_m23_T20_down@result %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::select(ID, p.adjust, Count) %>%
  mutate(direction = "Downregulated")

all_padj_m23 <- bind_rows(m23_up_padj, m23_down_padj)

# Combine reduced terms and summarise
reduced_m23_all <- bind_rows(reduced_m23_up, reduced_m23_down)

plot_df_m23_dir <- reduced_m23_all %>%
  left_join(all_padj_m23, by = c("go" = "ID", "direction")) %>%
  dplyr::group_by(parentTerm, direction) %>%
  dplyr::summarise(
    mean_padj   = mean(p.adjust, na.rm = TRUE),
    total_count = sum(Count, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(
    neg_log10_padj = -log10(mean_padj),
    direction = factor(direction, levels = c("Upregulated", "Downregulated"))
  )

# Plot
ggplot(plot_df_m23_dir,
       aes(x = direction,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           fill = neg_log10_padj)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = NULL,
    y    = NULL,
    fill = "-log10(p.adjust)",
    size = "Gene count",
    title = "m23 unique T20 Botrytis DEGs - GO enrichment by direction"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

pdf("m23_T20_unique_GO_dotplot_by_direction.pdf", width = 5, height = 6)
ggplot(plot_df_m23_dir,
       aes(x = direction,
           y = reorder(parentTerm, neg_log10_padj),
           size = total_count,
           fill = neg_log10_padj)) +
  geom_point(shape = 21, color = "black", stroke = 0.3) +
  scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x    = NULL,
    y    = NULL,
    fill = "-log10(p.adjust)",
    size = "Gene count",
    title = "m23 unique T20 Botrytis DEGs - GO enrichment by direction"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

write_xlsx(
  list(
    "m23_T20_up_GO"   = go_m23_T20_up@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID),
    "m23_T20_down_GO" = go_m23_T20_down@result %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, geneID)
  ),
  "GO_m23_T20_unique_by_direction.xlsx"
)

### ============================================================
###  Reorder genotype axis: WT -> L1 -> L14 -> T23
###  Paste this block at the END of PCA_type_analysis.R
###  Produces new PDFs with the updated genotype order.
### ============================================================


## ── 1. NEW COLUMN ORDER FOR HEATMAPS ─────────────────────────
## Swap L14 and T23 relative to the col_order used earlier.

col_order_v2 <- c(
  "WT_M16",  "WT_B16",  "WT_M20",  "WT_B20",
  "L1_M16",  "L1_B16",  "L1_M20",  "L1_B20",
  "L14_M16", "L14_B16", "L14_M20", "L14_B20",
  "T23_M16", "T23_B16", "T23_M20", "T23_B20"
)

# Reorder column annotation
ann_col_ord_v2 <- ann_col[col_order_v2, ]

# Reorder expression matrices
dat_scaled_geno_ord_v2 <- dat_scaled_geno_ord[, col_order_v2]
dat_scaled_time_ord_v2 <- dat_scaled_time_ord[, col_order_v2]

# Reorder shared DEG matrix (if dat_scaled_shared_ord exists)
if (exists("dat_scaled_shared_ord")) {
  dat_scaled_shared_ord_v2 <- dat_scaled_shared_ord[, col_order_v2]
}


## ── 2. HEATMAPS — GENOTYPE-DRIVEN ────────────────────────────

pdf("genotype_driven_heatmap_v2.pdf", width = 4, height = 6)
pheatmap(dat_scaled_geno_ord_v2,
         annotation_col    = ann_col_ord_v2,
         annotation_row    = row_ann_geno_simple,   # row annotations unchanged
         annotation_colors = ann_colors_geno_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Genotype-driven clusters (11, 4, 6)"
)
dev.off()


## ── 3. HEATMAPS — TIMEPOINT-DRIVEN ───────────────────────────

pdf("timepoint_driven_heatmap_supp_v2.pdf", width = 4, height = 6)
pheatmap(dat_scaled_time_ord_v2,
         annotation_col    = ann_col_ord_v2,
         annotation_row    = row_ann_time_simple,
         annotation_colors = ann_colors_time_simple,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Timepoint-driven clusters (3, 7, 10)"
)
dev.off()


## ── 4. HEATMAP — SHARED BOTRYTIS DEGs (optional) ─────────────
## Only runs if dat_scaled_shared_ord was built earlier in the script.

if (exists("dat_scaled_shared_ord_v2")) {
  pdf("shared_botrytis_degs_heatmap_v2.pdf", width = 4, height = 6)
  pheatmap(dat_scaled_shared_ord_v2,
           annotation_col    = ann_col_ord_v2,
           annotation_row    = row_ann_shared_ord,
           annotation_colors = ann_colors_shared,
           cluster_rows      = FALSE,
           cluster_cols      = FALSE,
           show_rownames     = FALSE,
           show_colnames     = TRUE,
           color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           main              = "Shared Botrytis DEGs (n = 1677)"
  )
  dev.off()
}


## ── 5. LINE PLOTS — GENOTYPE-DRIVEN CLUSTERS ─────────────────
## Rebuild cl_geno with the corrected genotype factor order.

genotype_levels_v2 <- c("WT", "L1", "L14", "T23")

cl_geno_v2 <- lapply(seq_along(genotype_clusters), function(x) {
  cl_num_lists_sum[[genotype_clusters[x]]] %>%
    mutate(cluster_name = factor(genotype_names[x], levels = genotype_names))
}) %>%
  list_rbind() %>%
  mutate(genotype = factor(genotype, levels = genotype_levels_v2))

p_genotype_v2 <- ggplot(cl_geno_v2, aes(x = timepoint, y = mean,
                                        colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Genotype-driven clusters")

pdf("genotype_driven_lineplots_v2.pdf", width = 4, height = 6)
p_genotype_v2
dev.off()


## ── 6. LINE PLOTS — TIMEPOINT-DRIVEN CLUSTERS ────────────────

cl_time_v2 <- cl_time %>%
  mutate(genotype = factor(genotype, levels = genotype_levels_v2))

p_timepoint_v2 <- ggplot(cl_time_v2, aes(x = timepoint, y = mean,
                                         colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci), alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Timepoint-driven clusters")

pdf("timepoint_driven_lineplots_supp_v2.pdf", width = 4, height = 6)
p_timepoint_v2
dev.off()


## ── 7. GO DOT PLOTS — REORDER GENOTYPE AXIS ──────────────────
## The x-axis here is the gene set (Conserved / mutant-specific),
## so "L14 specific" moves before "m23 specific" to match the
## WT -> L1 -> L14 -> T23 order your collaborators requested.

go_levels_v2 <- c("Conserved (WT)", "L1 specific", "L14 specific", "m23 specific")

# Top-10-per-category plot
plot_df_top_genotype_v2 <- plot_df_genotype %>%
  dplyr::group_by(response) %>%
  dplyr::slice_min(mean_padj, n = 10) %>%
  dplyr::ungroup() %>%
  mutate(response = factor(response, levels = go_levels_v2))

pdf("GO_botrytis_by_genotype_v2.pdf", width = 5, height = 8)
ggplot(plot_df_top_genotype_v2,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size  = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count",
    title = "GO enrichment — Botrytis response by genotype"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

# Manual / curated term subset plot
plot_df_manual_v2 <- plot_df_genotype %>%
  dplyr::filter(parentTerm %in% keep_terms) %>%
  mutate(response = factor(response, levels = go_levels_v2))

pdf("GO_botrytis_by_genotype_sub1_v2.pdf", width = 4.5, height = 6.5)
ggplot(plot_df_manual_v2,
       aes(x = response,
           y = reorder(parentTerm, neg_log10_padj),
           size  = total_count,
           color = neg_log10_padj)) +
  geom_point() +
  scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
  scale_size_continuous(range = c(2, 8)) +
  labs(
    x     = NULL,
    y     = NULL,
    color = "-log10(p.adjust)",
    size  = "Gene count"
  ) +
  theme_bw() +
  theme(
    axis.text.y     = element_text(size = 8),
    axis.text.x     = element_text(size = 10),
    panel.grid      = element_blank(),
    legend.position = "right"
  )
dev.off()

message("Done — new PDFs written with genotype order WT -> L1 -> L14 -> T23.")


#########################
## ── Clusters 1, 5, 8 — heatmap and line plots ─────────────────

clusters_185 <- c(1, 5, 8)
names_185    <- paste0("Cluster ", clusters_185,
                       " (n = ", as.vector(cl_sizes[clusters_185]), ")")

## ── Line plot ─────────────────────────────────────────────────

cl_185 <- lapply(seq_along(clusters_185), function(x) {
  cl_num_lists_sum[[clusters_185[x]]] %>%
    mutate(cluster_name = factor(names_185[x], levels = names_185))
}) %>%
  list_rbind() %>%
  mutate(genotype = factor(genotype, levels = genotype_levels_v2))

p_185 <- ggplot(cl_185, aes(x = timepoint, y = mean,
                            colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - ci, ymax = mean + ci),
              alpha = 0.3, color = NA) +
  scale_color_manual(values = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  scale_fill_manual(values  = c("Botrytis" = "#4a0000", "Mock" = "#4a90d9")) +
  theme_bw() +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, switch = "y") +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 12),
    strip.placement  = "outside",
    strip.background = element_blank()
  ) +
  labs(color = "Treatment", fill = "Treatment",
       title = "Clusters 1, 5, 8")

p_185

pdf("clusters_185_lineplots.pdf", width = 4, height = 6)
p_185
dev.off()


## ── Heatmap ───────────────────────────────────────────────────

# Build scaled matrix
genes_185      <- modules$Gene_ID[modules$cluster %in% clusters_185]
dat_scaled_185 <- t(scale(t(dat_log[rownames(dat_log) %in% genes_185, ])))
dat_scaled_185 <- dat_scaled_185[complete.cases(dat_scaled_185), ]

# Row annotation
row_ann_185 <- modules %>%
  filter(cluster %in% clusters_185) %>%
  filter(Gene_ID %in% rownames(dat_scaled_185)) %>%
  mutate(cluster = factor(cluster, levels = c("1", "5", "8"))) %>%
  column_to_rownames("Gene_ID") %>%
  rename(Cluster = cluster)

# Order genes — cluster 1, then 5, then 8 top to bottom,
# with within-cluster correlation ordering
gene_order_185 <- unlist(lapply(c("1", "5", "8"), function(cl) {
  genes <- rownames(row_ann_185)[as.character(row_ann_185$Cluster) == cl]
  if (length(genes) < 2) return(genes)
  sub_mat  <- dat_scaled_185[genes, ]
  dist_mat <- as.dist(1 - cor(t(sub_mat)))
  hc       <- hclust(dist_mat, method = "ward.D2")
  genes[hc$order]
}))

dat_scaled_185_ord <- dat_scaled_185[gene_order_185, col_order_v2]
row_ann_185_ord    <- row_ann_185[gene_order_185, , drop = FALSE]

# Annotation colors
ann_colors_185 <- ann_colors
ann_colors_185$Cluster <- c(
  "1" = "#d9d9d9",
  "5" = "#969696",
  "8" = "#525252"
)
ann_colors_185$treatment <- c("Mock" = "#d9d9d9", "Botrytis" = "#4a0000")

# Plot
pheatmap(dat_scaled_185_ord,
         annotation_col    = ann_col_ord_v2,
         annotation_row    = row_ann_185_ord,
         annotation_colors = ann_colors_185,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Clusters 1, 5, 8"
)

pdf("clusters_185_heatmap.pdf", width = 4, height = 6)
pheatmap(dat_scaled_185_ord,
         annotation_col    = ann_col_ord_v2,
         annotation_row    = row_ann_185_ord,
         annotation_colors = ann_colors_185,
         cluster_rows      = FALSE,
         cluster_cols      = FALSE,
         show_rownames     = FALSE,
         show_colnames     = TRUE,
         color             = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         main              = "Clusters 1, 5, 8"
)
dev.off()


## ── GO by genotype, timepoint, and direction ─────────────────────────────
## Produces 16 enrichments:
## 4 categories × 2 timepoints × 2 directions

library(clusterProfiler)
library(org.At.tair.db)

# ── Step 1: Build per-timepoint, per-direction gene sets ──────────────────

timepoints <- c("T16", "T20")
directions <- c("Upregulated", "Downregulated")

# Helper — pull gene IDs from bc with optional direction filter
get_bc_genes <- function(genotype, timepoint, direction = NULL) {
  sub <- bc[bc$Genotype == genotype & bc$Timepoint == timepoint, ]
  if (!is.null(direction)) sub <- sub[sub$DE_Direction == direction, ]
  unique(sub$Gene_ID)
}

# Build all 16 gene sets as a flat named list
gene_sets_dir <- list()

for (tp in timepoints) {
  
  # Full gene lists per genotype at this timepoint
  wt_all  <- get_bc_genes("WT",  tp)
  l1_all  <- get_bc_genes("L1",  tp)
  l14_all <- get_bc_genes("L14", tp)
  t23_all <- get_bc_genes("T23", tp)
  
  # Category membership (same logic as before, now per timepoint)
  conserved_all <- intersect(wt_all, union(l1_all, union(l14_all, t23_all)))
  l1_spec_all   <- setdiff(l1_all,  wt_all)
  l14_spec_all  <- setdiff(l14_all, wt_all)
  t23_spec_all  <- setdiff(t23_all, wt_all)
  
  for (dir in directions) {
    
    # Direction-filtered gene lists
    wt_dir  <- get_bc_genes("WT",  tp, dir)
    l1_dir  <- get_bc_genes("L1",  tp, dir)
    l14_dir <- get_bc_genes("L14", tp, dir)
    t23_dir <- get_bc_genes("T23", tp, dir)
    
    # Intersect category membership with direction
    # Conserved: use WT direction as reference
    # Specific:  use the relevant mutant's direction
    dir_short <- ifelse(dir == "Upregulated", "Up", "Down")
    
    gene_sets_dir[[paste("Conserved (WT)", tp, dir_short, sep = "_")]] <-
      intersect(conserved_all, wt_dir)
    
    gene_sets_dir[[paste("L1 specific",  tp, dir_short, sep = "_")]] <-
      intersect(l1_spec_all,  l1_dir)
    
    gene_sets_dir[[paste("L14 specific", tp, dir_short, sep = "_")]] <-
      intersect(l14_spec_all, l14_dir)
    
    gene_sets_dir[[paste("T23 specific", tp, dir_short, sep = "_")]] <-
      intersect(t23_spec_all, t23_dir)
  }
}

# Check sizes before running GO
sapply(gene_sets_dir, length)


# ── Step 2: Run GO enrichment on all 16 sets ─────────────────────────────

go_by_genotype_dir <- lapply(names(gene_sets_dir), function(gs) {
  genes <- gene_sets_dir[[gs]]
  
  if (length(genes) == 0) {
    message(paste("Skipping", gs, "— no genes"))
    return(NULL)
  }
  
  enrichGO(
    gene          = genes,
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,     # keep all terms — filter at plot stage
    qvalueCutoff  = 1,
    readable      = TRUE
  )
})
names(go_by_genotype_dir) <- names(gene_sets_dir)

# Remove any NULLs (empty gene sets)
go_by_genotype_dir <- Filter(Negate(is.null), go_by_genotype_dir)

# Check how many significant terms per set at p < 0.05
sapply(go_by_genotype_dir, function(x)
  nrow(x@result[x@result$pvalue < 0.05, ]))

## ── Dot plots — GO by genotype, timepoint, direction ─────────────────────

# Category order on x-axis (matches genotype order in rest of figures)
category_levels <- c("Conserved (WT)", "L1 specific", "L14 specific", "T23 specific")

# Build a single long dataframe from all 16 enrichments
go_dir_long <- lapply(names(go_by_genotype_dir), function(gs) {
  res <- go_by_genotype_dir[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice_head(n = 10)     # top 10 terms per set
  
  if (nrow(res) == 0) return(NULL)
  
  # Parse set name back into components
  parts <- strsplit(gs, "_")[[1]]
  # names are like "Conserved (WT)_T16_Up" — last two parts are timepoint and direction
  # category is everything before the last two underscores
  direction <- parts[length(parts)]
  timepoint <- parts[length(parts) - 1]
  category  <- paste(parts[1:(length(parts) - 2)], collapse = "_")
  
  res %>%
    dplyr::select(Description, pvalue, Count) %>%
    dplyr::mutate(
      neg_log10_p = -log10(pvalue),
      category    = category,
      timepoint   = timepoint,
      direction   = direction
    )
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    category = factor(category, levels = category_levels),
    direction = factor(direction, levels = c("Up", "Down"))
  )


# ── Plot function — one panel per timepoint × direction ───────────────────

make_go_dotplot <- function(tp, dir) {
  
  df <- go_dir_long %>%
    dplyr::filter(timepoint == tp, direction == dir)
  
  if (nrow(df) == 0) {
    message(paste("No terms for", tp, dir))
    return(NULL)
  }
  
  # Order terms by mean significance across categories
  term_order <- df %>%
    dplyr::group_by(Description) %>%
    dplyr::summarise(mean_nlp = mean(neg_log10_p), .groups = "drop") %>%
    dplyr::arrange(mean_nlp) %>%
    dplyr::pull(Description)
  
  df <- df %>%
    dplyr::mutate(Description = factor(Description, levels = term_order))
  
  ggplot(df, aes(x = category,
                 y = Description,
                 size = Count,
                 color = neg_log10_p)) +
    geom_point() +
    scale_color_gradient(low = "#4a90d9", high = "#b2182b") +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x     = NULL,
      y     = NULL,
      color = "-log10(p)",
      size  = "Gene count",
      title = paste0("GO enrichment — ", tp, " ", dir, "regulated")
    ) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(size = 8, angle = 30, hjust = 1),
      axis.text.y     = element_text(size = 8),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

# ── Generate and export all 4 plots ──────────────────────────────────────

p_T16_Up   <- make_go_dotplot("T16", "Up")
p_T16_Down <- make_go_dotplot("T16", "Down")
p_T20_Up   <- make_go_dotplot("T20", "Up")
p_T20_Down <- make_go_dotplot("T20", "Down")

p_T16_Up
p_T16_Down
p_T20_Up
p_T20_Down

pdf("GO_by_genotype_direction_T16_Up.pdf",   width = 7, height = 8)
p_T16_Up
dev.off()

pdf("GO_by_genotype_direction_T16_Down.pdf", width = 7, height = 8)
p_T16_Down
dev.off()

pdf("GO_by_genotype_direction_T20_Up.pdf",   width = 7, height = 8)
p_T20_Up
dev.off()

pdf("GO_by_genotype_direction_T20_Down.pdf", width = 7, height = 8)
p_T20_Down
dev.off()

## ── rrvgo semantic reduction on all 16 GO sets ───────────────────────────

library(rrvgo)

# Run rrvgo on each set — same tryCatch pattern as earlier in pipeline
go_dir_reduced <- lapply(names(go_by_genotype_dir), function(gs) {
  
  res <- go_by_genotype_dir[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05)
  
  if (nrow(res) == 0) {
    message(paste("No significant terms for", gs))
    return(NULL)
  }
  
  simMatrix <- tryCatch(
    calculateSimMatrix(
      res$ID,
      orgdb  = "org.At.tair.db",
      ont    = "BP",
      method = "Rel"
    ),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      calculateSimMatrix(
        res$ID,
        orgdb  = "org.At.tair.db",
        ont    = "BP",
        method = "Rel"
      )
    )
  )
  
  if (is.null(simMatrix) || !is.matrix(simMatrix) || nrow(simMatrix) < 2) {
    message(paste("rrvgo failed for", gs))
    return(NULL)
  }
  
  scores <- setNames(-log10(res$pvalue), res$ID)
  scores <- scores[names(scores) %in% rownames(simMatrix)]
  
  if (length(scores) < 2) return(NULL)
  
  tryCatch({
    reducedTerms <- reduceSimMatrix(
      simMatrix,
      scores,
      threshold = 0.7,
      orgdb     = "org.At.tair.db"
    )
    reducedTerms %>% dplyr::mutate(set_name = gs)
  }, error = function(e) {
    message(paste("reduceSimMatrix failed for", gs))
    NULL
  })
})
names(go_dir_reduced) <- names(go_by_genotype_dir)
go_dir_reduced <- Filter(Negate(is.null), go_dir_reduced)

# Check how many representative parent terms per set
sapply(go_dir_reduced, function(x) length(unique(x$parentTerm)))


## ── Build long dataframe from reduced terms ───────────────────────────────
## For each parent term, summarise mean pvalue and total gene count
## across all child terms in that group

# First pull pvalue and Count from original results for joining
all_padj_dir <- lapply(names(go_by_genotype_dir), function(gs) {
  go_by_genotype_dir[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::select(ID, pvalue, Count) %>%
    dplyr::mutate(set_name = gs)
}) %>%
  dplyr::bind_rows()

# Join reduced terms to pvalue/count, then summarise per parent term
go_dir_reduced_long <- dplyr::bind_rows(go_dir_reduced) %>%
  dplyr::left_join(all_padj_dir, by = c("go" = "ID", "set_name")) %>%
  dplyr::group_by(set_name, parentTerm) %>%
  dplyr::summarise(
    mean_pvalue = mean(pvalue, na.rm = TRUE),
    total_count = sum(Count,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(neg_log10_p = -log10(mean_pvalue)) %>%
  # Parse set name into components
  dplyr::mutate(
    direction = dplyr::case_when(
      grepl("_Up$",   set_name) ~ "Up",
      grepl("_Down$", set_name) ~ "Down"
    ),
    timepoint = dplyr::case_when(
      grepl("_T16_", set_name) ~ "T16",
      grepl("_T20_", set_name) ~ "T20"
    ),
    category = gsub("_T16_Up|_T16_Down|_T20_Up|_T20_Down", "", set_name)
  ) %>%
  dplyr::mutate(
    category  = factor(category,  levels = category_levels),
    direction = factor(direction, levels = c("Up", "Down"))
  )


## ── Dot plot function — same structure, now using parent terms ────────────

make_go_dotplot_reduced <- function(tp, dir) {
  
  df <- go_dir_reduced_long %>%
    dplyr::filter(timepoint == tp, direction == dir) %>%
    dplyr::group_by(category) %>%
    dplyr::slice_min(mean_pvalue, n = 15) %>%
    dplyr::ungroup()
  
  if (nrow(df) == 0) {
    message(paste("No reduced terms for", tp, dir))
    return(NULL)
  }
  
  term_order <- df %>%
    dplyr::group_by(parentTerm) %>%
    dplyr::summarise(mean_nlp = mean(neg_log10_p), .groups = "drop") %>%
    dplyr::arrange(mean_nlp) %>%
    dplyr::pull(parentTerm)
  
  df <- df %>%
    dplyr::mutate(parentTerm = factor(parentTerm, levels = term_order))
  
  ggplot(df, aes(x = category, y = parentTerm, size = total_count, fill = neg_log10_p)) +
    geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.7) +
    scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x     = NULL,
      y     = NULL,
      fill = "-log10(p)",
      size  = "Gene count",
      title = paste0("GO enrichment — ", tp, " ", dir, "regulated")
    ) +
    theme_bw() +
    theme(
      axis.text.x     = element_text(size = 8, angle = 30, hjust = 1),
      axis.text.y     = element_text(size = 8),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

## ── Generate and export ───────────────────────────────────────────────────

p_T16_Up_red   <- make_go_dotplot_reduced("T16", "Up")
p_T16_Down_red <- make_go_dotplot_reduced("T16", "Down")
p_T20_Up_red   <- make_go_dotplot_reduced("T20", "Up")
p_T20_Down_red <- make_go_dotplot_reduced("T20", "Down")

p_T16_Up_red
p_T16_Down_red
p_T20_Up_red
p_T20_Down_red

pdf("GO_reduced_T16_Up.pdf",   width = 7, height = 8)
p_T16_Up_red
dev.off()

pdf("GO_reduced_T16_Down.pdf", width = 7, height = 8)
p_T16_Down_red
dev.off()

pdf("GO_reduced_T20_Up.pdf",   width = 7, height = 8)
p_T20_Up_red
dev.off()

pdf("GO_reduced_T20_Down.pdf", width = 7, height = 8)
p_T20_Down_red
dev.off()

## ── Export all 16 GO sets to a single Excel file ─────────────────────────

library(writexl)

go_dir_tables <- lapply(names(go_by_genotype_dir), function(gs) {
  go_by_genotype_dir[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)
})
names(go_dir_tables) <- c(
  "Cons_WT_T16_Up",   "L1_T16_Up",   "L14_T16_Up",   "T23_T16_Up",
  "Cons_WT_T16_Down", "L1_T16_Down", "L14_T16_Down", "T23_T16_Down",
  "Cons_WT_T20_Up",   "L1_T20_Up",   "L14_T20_Up",   "T23_T20_Up",
  "Cons_WT_T20_Down", "L1_T20_Down", "L14_T20_Down", "T23_T20_Down"
)
names(go_dir_tables) <- names(go_by_genotype_dir)
write_xlsx(go_dir_tables, "GO_by_genotype_direction.xlsx")

## ── Export genotype-specific DEG lists for Cytoscape ─────────────────────

library(writexl)

# Convert gene_sets_dir to a long dataframe with classification columns
gene_sets_long <- lapply(names(gene_sets_dir), function(gs) {
  genes <- gene_sets_dir[[gs]]
  if (length(genes) == 0) return(NULL)
  
  parts     <- strsplit(gs, "_")[[1]]
  direction <- parts[length(parts)]
  timepoint <- parts[length(parts) - 1]
  category  <- paste(parts[1:(length(parts) - 2)], collapse = "_")
  
  data.frame(
    Gene_ID   = genes,
    category  = category,
    timepoint = timepoint,
    direction = direction,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

# Join back to bc to recover log2FoldChange and stats
# Match on Gene_ID + the relevant genotype + timepoint
# For conserved: use WT stats as reference
# For specific sets: use the relevant mutant stats

gene_sets_annotated <- gene_sets_long %>%
  dplyr::mutate(
    genotype_source = dplyr::case_when(
      category == "Conserved (WT)" ~ "WT",
      category == "L1 specific"    ~ "L1",
      category == "L14 specific"   ~ "L14",
      category == "T23 specific"   ~ "T23"
    )
  ) %>%
  dplyr::left_join(
    bc %>% dplyr::select(Gene_ID, Genotype, Timepoint,
                         log2FoldChange, pvalue, padj, DE_Direction),
    by = c("Gene_ID"        = "Gene_ID",
           "genotype_source" = "Genotype",
           "timepoint"       = "Timepoint")
  ) %>%
  dplyr::select(Gene_ID, category, timepoint, direction,
                log2FoldChange, pvalue, padj, DE_Direction, genotype_source) %>%
  dplyr::arrange(timepoint, category, direction, dplyr::desc(abs(log2FoldChange)))

# Check structure
head(gene_sets_annotated)
dim(gene_sets_annotated)


## ── Export 1: Single CSV for Cytoscape node attribute import ─────────────
write.csv(gene_sets_annotated,
          "cytoscape_genotype_specific_DEGs.csv",
          row.names = FALSE)


## ── Export 2: Excel with one sheet per timepoint × direction ─────────────

excel_sheets <- list(
  "T16_Upregulated"   = gene_sets_annotated %>%
    dplyr::filter(timepoint == "T16", direction == "Up"),
  "T16_Downregulated" = gene_sets_annotated %>%
    dplyr::filter(timepoint == "T16", direction == "Down"),
  "T20_Upregulated"   = gene_sets_annotated %>%
    dplyr::filter(timepoint == "T20", direction == "Up"),
  "T20_Downregulated" = gene_sets_annotated %>%
    dplyr::filter(timepoint == "T20", direction == "Down")
)

# Quick summary of how many genes per sheet
sapply(excel_sheets, nrow)

write_xlsx(excel_sheets, "cytoscape_genotype_specific_DEGs.xlsx")


## ── GO on clusters 1, 5, 8 split by direction ────────────────────────────

# Re-add Genotype column to mut if not already present from earlier in session
if (!"Genotype" %in% colnames(mut)) {
  mut$Genotype <- dplyr::recode(mut$contrast,
                                "tpltpr1"     = "L1",
                                "tpltpr1tpr4" = "L14",
                                "tpr2tpr3"    = "T23")
}

# Assign direction to each gene using mean log2FC across all mutant contrasts

clusters_185 <- c(1, 5, 8)

## ── GO enrichment on clusters 1, 5, 8 ────────────────────────────────────

clusters_185 <- c(1, 5, 8)

# Build one gene list per cluster
go_185_lists <- lapply(clusters_185, function(cl) {
  modules$Gene_ID[modules$cluster == cl]
})
names(go_185_lists) <- paste0("Cluster_", clusters_185)

sapply(go_185_lists, length)

# Run GO enrichment
go_185_results <- lapply(names(go_185_lists), function(gs) {
  enrichGO(
    gene          = go_185_lists[[gs]],
    universe      = background_genes,
    OrgDb         = org.At.tair.db,
    keyType       = "TAIR",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
})
names(go_185_results) <- names(go_185_lists)

# Check significant terms per cluster
sapply(go_185_results, function(x)
  nrow(x@result[x@result$pvalue < 0.05, ]))


## ── rrvgo reduction ───────────────────────────────────────────────────────

go_185_reduced <- lapply(names(go_185_results), function(gs) {
  
  res <- go_185_results[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05)
  
  if (nrow(res) == 0) {
    message(paste("No significant terms for", gs))
    return(NULL)
  }
  
  simMatrix <- tryCatch(
    calculateSimMatrix(
      res$ID,
      orgdb  = "org.At.tair.db",
      ont    = "BP",
      method = "Rel"
    ),
    error   = function(e) NULL,
    warning = function(w) suppressWarnings(
      calculateSimMatrix(
        res$ID,
        orgdb  = "org.At.tair.db",
        ont    = "BP",
        method = "Rel"
      )
    )
  )
  
  if (is.null(simMatrix) || !is.matrix(simMatrix) || nrow(simMatrix) < 2) {
    message(paste("rrvgo failed for", gs))
    return(NULL)
  }
  
  scores <- setNames(-log10(res$pvalue), res$ID)
  scores <- scores[names(scores) %in% rownames(simMatrix)]
  
  if (length(scores) < 2) return(NULL)
  
  tryCatch({
    reduceSimMatrix(
      simMatrix,
      scores,
      threshold = 0.7,
      orgdb     = "org.At.tair.db"
    ) %>%
      dplyr::mutate(set_name = gs)
  }, error = function(e) {
    message(paste("reduceSimMatrix failed for", gs))
    NULL
  })
})
names(go_185_reduced) <- names(go_185_results)
go_185_reduced <- Filter(Negate(is.null), go_185_reduced)

sapply(go_185_reduced, function(x) length(unique(x$parentTerm)))


## ── Build long dataframe ──────────────────────────────────────────────────

all_padj_185 <- lapply(names(go_185_results), function(gs) {
  go_185_results[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::select(ID, pvalue, Count) %>%
    dplyr::mutate(set_name = gs)
}) %>%
  dplyr::bind_rows()

go_185_reduced_long <- dplyr::bind_rows(go_185_reduced) %>%
  dplyr::left_join(all_padj_185, by = c("go" = "ID", "set_name")) %>%
  dplyr::group_by(set_name, parentTerm) %>%
  dplyr::summarise(
    mean_pvalue = mean(pvalue,  na.rm = TRUE),
    total_count = sum(Count,    na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(
    neg_log10_p = -log10(mean_pvalue),
    cluster     = factor(set_name, levels = paste0("Cluster_", clusters_185))
  )


## ── One dot plot per cluster ──────────────────────────────────────────────

make_go_185_dotplot <- function(cl) {
  
  df <- go_185_reduced_long %>%
    dplyr::filter(cluster == cl) %>%
    dplyr::slice_min(mean_pvalue, n = 15) %>%
    dplyr::mutate(parentTerm = reorder(parentTerm, neg_log10_p))
  
  if (nrow(df) == 0) {
    message(paste("No terms for", cl))
    return(NULL)
  }
  
  ggplot(df, aes(x = neg_log10_p,
                 y = parentTerm,
                 size  = total_count,
                 fill  = neg_log10_p)) +
    geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.7) +
    scale_fill_gradient(low = "#4a90d9", high = "#b2182b") +
    scale_size_continuous(range = c(2, 8)) +
    labs(
      x     = "-log10(p)",
      y     = NULL,
      fill  = "-log10(p)",
      size  = "Gene count",
      title = cl
    ) +
    theme_bw() +
    theme(
      axis.text.y     = element_text(size = 8),
      axis.text.x     = element_text(size = 8),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

p_cl1 <- make_go_185_dotplot("Cluster_1")
p_cl5 <- make_go_185_dotplot("Cluster_5")
p_cl8 <- make_go_185_dotplot("Cluster_8")

p_cl1
p_cl5
p_cl8

pdf("GO_Cluster1.pdf", width = 6, height = 6)
p_cl1
dev.off()

pdf("GO_Cluster5.pdf", width = 6, height = 6)
p_cl5
dev.off()

pdf("GO_Cluster8.pdf", width = 6, height = 6)
p_cl8
dev.off()


## ── Export to Excel ───────────────────────────────────────────────────────

go_185_tables <- lapply(names(go_185_results), function(gs) {
  go_185_results[[gs]]@result %>%
    dplyr::filter(pvalue < 0.05) %>%
    dplyr::arrange(pvalue) %>%
    dplyr::select(Description, GeneRatio, BgRatio, pvalue, p.adjust, Count, geneID)
})
names(go_185_tables) <- names(go_185_results)

write_xlsx(go_185_tables, "GO_clusters_185.xlsx")
pdf("GO_clusters_185_Down.pdf", width = 6, height = 8)
p_185_Down
dev.off()
