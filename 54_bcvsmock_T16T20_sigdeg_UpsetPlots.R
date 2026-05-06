## redo UPSET plots
## making upset plots but for the bcvsmock sepTPs DEGs 
library(tidyverse)
library(UpSetR)

getwd()
setwd("/03_results/14_REDO_bcvsmock_tps")
# data loading ------------------------------------------------------------
files <- list.files(pattern = "\\_at_sig_degs.csv", full.names = TRUE) #this searches for all .csv files in the folder
files

df_names <- files %>% 
  str_replace_all(c("./" ="",
                    "_at_sig_degs.csv" = "")) #this removes the / and .csv from all the file names, leaving just the proper names

print(df_names)
df_names <- as.vector(df_names)

data <- lapply(files, FUN = function(x)read.csv(x)) %>% #read in the csvs
  set_names(df_names)
## separating upregulated genes
data_upreg <- lapply(seq_along(data), function(y) {
  data[[y]] <- data[[y]] %>% 
    filter(DE_Direction == "Upregulated") %>% 
    mutate(contrast = paste0(Genotype, ", ", Timepoint))
})
names(data_upreg) <- df_names

data_upreg <- lapply(seq_along(data_upreg), function(i){
  data_upreg[[i]] <- data_upreg[[i]] %>% 
    mutate(ones = 1) %>% 
    pivot_wider(names_from = contrast, values_from = ones) %>% 
    as.data.frame() %>% select(2, 12)
})
names(data_upreg) <- paste0(df_names, "_upreg")

## separating downregulated genes
data_downreg <- lapply(seq_along(data), function(y) {
  data[[y]] <- data[[y]] %>% 
    filter(DE_Direction == "Downregulated") %>% 
    mutate(contrast = paste0(Genotype, ", ", Timepoint))
})
names(data_downreg) <- df_names

data_downreg <- lapply(seq_along(data_downreg), function(i){
  data_downreg[[i]] <- data_downreg[[i]] %>% 
    mutate(ones = 1) %>% 
    pivot_wider(names_from = contrast, values_from = ones) %>% 
    as.data.frame() %>% select(2, 12)
})
names(data_downreg) <- paste0(df_names, "_downreg")

## making a binary table
all_data_downreg <- data_downreg %>% reduce(full_join) %>% replace(is.na(.), 0)
all_data_upreg <- data_upreg %>% reduce(full_join) %>% replace(is.na(.), 0)

rownames(all_data_downreg) <- all_data_downreg$Gene_ID
all_data_downreg <- select(all_data_downreg, -1)

rownames(all_data_upreg) <- all_data_upreg$Gene_ID
all_data_upreg <- select(all_data_upreg, -1)

## plotting the upset plots
# first the upregulated genes

## combination matrix
#T16, upreg
T16_upreg <- select(all_data_upreg, c(1, 3, 5, 7))

png("../../04_figures/REDO_bcvsmock_sepTPs_upsetplot_T16upregDEGs.png",
    width = 3000,      # width in pixels
    height = 2500,     # height in pixels
    res = 300) 
upset(T16_upreg, order.by = "freq",
      text.scale = c(2, 2, 1, 1, 2, 3),
      point.size = 4,
      line.size = 1, 
      main.bar.color = "#1e81b0", 
      shade.color = "#1e81b0", 
      shade.alpha = 0.1, 
      sets.bar.color ="#1e81b0")


dev.off()

#T20
T20_upreg <- select(all_data_upreg, c(2, 4, 6, 8))

png("../../04_figures/REDO_bcvsmock_sepTPs_upsetplot_T20upregDEGs.png",
    width = 3000,      # width in pixels
    height = 2500,     # height in pixels
    res = 300) 
upset(T20_upreg, order.by = "freq",
      text.scale = c(2, 2, 1, 1, 2, 3),
      point.size = 4,
      line.size = 1, 
      main.bar.color = "#1e81b0", 
      shade.color = "#1e81b0", 
      shade.alpha = 0.1, 
      sets.bar.color ="#1e81b0")
dev.off()

#T16, downreg
T16_downreg <- select(all_data_downreg, c(1, 3, 5, 7))

png("../../04_figures/REDO_bcvsmock_sepTPs_upsetplot_T16downregDEGs.png",
    width = 3000,      # width in pixels
    height = 2500,     # height in pixels
    res = 300) 
upset(T16_downreg, order.by = "freq",
      text.scale = c(2, 2, 2, 2, 2, 3),
      point.size = 4,
      line.size = 1, 
      main.bar.color = "#e28743", 
      shade.color = "#e28743", 
      shade.alpha = 0.1, 
      sets.bar.color ="#e28743")
dev.off()

#T20, downreg
T20_downreg <- select(all_data_downreg, c(2, 4, 6, 8))

png("../../04_figures/REDO_bcvsmock_sepTPs_upsetplot_T20downregDEGs.png",
    width = 3000,      # width in pixels
    height = 2500,     # height in pixels
    res = 300) 
upset(T20_downreg, order.by = "freq",
      text.scale = c(2, 2, 2, 2, 2, 3),
      point.size = 4,
      line.size = 1, 
      main.bar.color = "#e28743", 
      shade.color = "#e28743", 
      shade.alpha = 0.1, 
      sets.bar.color ="#e28743")
dev.off()
