## making cluster plots for the infected degs 

library(pheatmap)
library(tidyverse)

getwd()
setwd("03_results/14_REDO_bcvsmock_tps")

upreg_uncon_degs <- read.csv("bcvsmock_upreg_unconserved_sigdegs.csv")
downreg_uncon_degs <- read.csv("bcvsmock_downreg_unconserved_sigdegs.csv")

upreg_uncon_degs <- select(upreg_uncon_degs, Gene_ID)
downreg_uncon_degs <- select(downreg_uncon_degs, Gene_ID)
all_uncond_degs  <-  bind_rows(upreg_uncon_degs, downreg_uncon_degs) %>% 
  unique()

norm_counts <- read.csv("../../01_data/USETHIS_At_tplmuts_mean_norm_counts.csv") %>% 
  select(-1)

deg_counts <- all_uncond_degs %>% left_join(norm_counts, by = join_by(Gene_ID))

rownames(deg_counts) <- deg_counts$Gene_ID
deg_counts <- select(deg_counts, -1)
deg_counts <- deg_counts %>% select( WT_M16, WT_B16, WT_M20, WT_B20, L1_M16, L1_B16, L1_M20, L1_B20, L14_M16, L14_B16, L14_M20, L14_B20, T23_M16, T23_B16, T23_M20, T23_B20)
deg_counts_mat <- as.matrix(deg_counts)

deg_counts_mat <- t(scale(t(deg_counts_mat)))
deg_counts_mat <- na.omit(deg_counts_mat)

p <- pheatmap(deg_counts_mat, 
                cluster_cols = FALSE, 
                cluster_rows = TRUE,
                color = colorRampPalette(c("blue", "white", "red"))(100),
                show_rownames = FALSE, 
                cutree_rows = 6,
                annotation_row = all_rows,
                annotation_colors = ann_colours,
                annotation_legend = FALSE,
                treeheight_row = 0,
                gaps_col = c(2, 4, 6, 8, 10, 12, 14))
p

##saving clusters
cl <- cutree(p$tree_row,6)
ann <- data.frame(cl)
rownames(ann) <- rownames(deg_counts_mat)
write.csv(ann, "bcvsmock_T16T20_HUGE_pheatmap_clusters_unconserved.csv")
#row anns
degs <- read.csv("bcvsmock_allmuts_sigdegs.csv")
degs_T16 <- degs %>% 
  filter(Timepoint == "T16")
degs_T16 <- degs_T16 %>% 
  select(Gene_ID, DE_Direction, Genotype)
row_anns <- degs_T16 %>% 
  pivot_wider(names_from = Genotype, 
              values_from = DE_Direction)
row_anns[is.na(row_anns)] <- "Not DE"
row_anns <- as.data.frame(row_anns)
rownames(row_anns) <- row_anns$Gene_ID


row_anns <- select(row_anns,Gene_ID, T23, L14, L1, WT )

##
degs_T20 <- degs %>% 
  filter(Timepoint == "T20")
degs_T20 <- degs_T20 %>% 
  select(Gene_ID, DE_Direction, Genotype, Timepoint) %>% 
  transmute(Gene_ID, DE_Direction, Genotype = paste0(Genotype, "_", Timepoint))
row_anns_20 <- degs_T20 %>% 
  pivot_wider(names_from = Genotype, 
              values_from = DE_Direction)
row_anns_20[is.na(row_anns_20)] <- "Not DE"
row_anns_20 <- as.data.frame(row_anns_20)
rownames(row_anns_20) <- row_anns_20$Gene_ID


row_anns_20 <- select(row_anns_20, Gene_ID, T23_T20, L14_T20, L1_T20, WT_T20)


## full row_anns 
all_rows <- full_join(row_anns, row_anns_20, by = join_by(Gene_ID))
all_rows[is.na(all_rows)] <-  "Not DE"
rownames(all_rows) <- all_rows$Gene_ID

all_rows <- all_rows %>% select(-1)

all_rows <- select(all_rows, c("T23_T20", "T23", "L14_T20", "L14","L1_T20", "L1", "WT_T20", "WT"))

ann_colours <- list(
  WT = c("Upregulated" = "#1e81b0",
         "Downregulated" = "#e28743",
         "Not DE" = "#B2BEB5"),
  L1 = c("Upregulated" = "#1e81b0",
         "Downregulated" = "#e28743",
         "Not DE" = "#B2BEB5"),
  L14 = c("Upregulated" = "#1e81b0",
          "Downregulated" = "#e28743",
          "Not DE" = "#B2BEB5"),
  T23 = c("Upregulated" = "#1e81b0",
          "Downregulated" = "#e28743",
          "Not DE" = "#B2BEB5"),
  WT_T20 = c("Upregulated" = "#1e81b0",
         "Downregulated" = "#e28743",
         "Not DE" = "#B2BEB5"),
  L1_T20 = c("Upregulated" = "#1e81b0",
         "Downregulated" = "#e28743",
         "Not DE" = "#B2BEB5"),
  L14_T20 = c("Upregulated" = "#1e81b0",
          "Downregulated" = "#e28743",
          "Not DE" = "#B2BEB5"),
  T23_T20 = c("Upregulated" = "#1e81b0",
          "Downregulated" = "#e28743",
          "Not DE" = "#B2BEB5")
)


#### finished huge venn ]
png("../../04_figures/REDO_bcvsmock_T16T20_unconserved_degs_pheatmap_HUGE.png",
    width = 3000,      # width in pixels
    height = 2500,     # height in pixels
    res = 300)
pheatmap(deg_counts_mat, 
         cluster_cols = FALSE, 
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE, 
         cutree_rows = 6,
         annotation_row = all_rows,
         annotation_colors = ann_colours,
         annotation_legend = FALSE,
         treeheight_row = 0,
         gaps_col = c(2, 4, 6, 8, 10, 12, 14))
dev.off()

##
ann_summ <- read.csv("03_results/10_bcvsmock_timepoint_inclusive/bcvsmock_T16T20_HUGE_pheatmap_clusters_unconserved.csv") %>% 
  group_by(clall) %>% 
  summarise(n = n())
