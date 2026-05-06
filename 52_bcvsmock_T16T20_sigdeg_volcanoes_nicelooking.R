#REDo bcvsmock volcano plots
library(tidyverse)

getwd()
setwd("")


### dataaaa 
### reading all the data into this workspace
files <- list.files(pattern = "\\_allLFCs.csv$", full.names = TRUE) #this searches for all .csv files in the folder
files

df_names <- files %>% 
  str_replace_all(c("./" ="",
                    "_allLFCs.csv" = "")) #this removes the / and .csv from all the file names, leaving just the proper names

print(df_names)
df_names <- as.vector(df_names)

data <- lapply(files, FUN = function(x)read.csv(x)) %>% #read in the csvs
  set_names(df_names)


data <- lapply(seq_along(data), function(i){
  data[[i]] <- data[[i]][grep("^AT", data[[i]]$Gene_ID),]
})
names(data) <- df_names

all_contrast <- list_rbind(data)

data2 <- all_contrast %>% select(-1) %>% 
  mutate(log10_padj = -log10(padj))

###trying to rearrange data to show all the points 
## you want to have "Not DE" to be at the top of the table, and then the others
## the points at the top of the df will be plotted first, underneath all the other points

data2 <- data2 %>% mutate(order = case_when(
  DE_Direction == "Downregulated" ~ "first",
  DE_Direction == "Upregulated" ~ "first",
  DE_Direction == "Not DE" ~ "second"
)) %>% 
  arrange(desc(order))

all_volcano <- ggplot(data2, aes(x = log2FoldChange, y = log10_padj)) +
  geom_point(aes(color = ifelse(log2FoldChange <= -0.5 & log10_padj > -log10(0.01), "Downregulated", 
                                ifelse(log2FoldChange >= 0.5 & log10_padj > -log10(0.01), "Upregulated", "Not significant")))) +
  scale_color_manual(values = c("#e28743", "#1e81b0", "#B2BEB5"), 
                     breaks = c("Downregulated", "Upregulated", "Not DE")) +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") +
  geom_vline(xintercept = -0.5, 
             linetype = "dashed") +
  geom_vline(xintercept = 0.5, 
             linetype = "dashed") +
  labs(x = "log2FC", y = "-log10(p-value)") +
  facet_wrap(~ Genotype + Timepoint, nrow = 4) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 16))
all_volcano

ggsave("../../04_figures/REDO_bcvsmock_sigdeg_volcanoes_slim.png")
