
library(tidyverse)

getwd()
setwd("03_results/14_REDO_bcvsmock_tps")
files <- list.files(pattern = "\\_at_sig_degs.csv", full.names = TRUE) #this searches for all .csv files in the folder
files

df_names <- files %>% 
  str_replace_all(c("./" ="",
                    "_at_sig_degs.csv" = "")) #this removes the / and .csv from all the file names, leaving just the proper names

print(df_names)
df_names <- as.vector(df_names)

data <- lapply(files, FUN = function(x)read.csv(x)) %>% #read in the csvs
  set_names(df_names)

all_data <- list_rbind(data)
write.csv(all_data, "bcvsmock_allmuts_sigdegs.csv")

summary <- all_data %>% 
  group_by(Genotype, Timepoint, DE_Direction) %>% 
  summarise(count = n())

## making a bargraph 
bargraph <- ggplot(summary, aes(x = Timepoint, y = count, fill = DE_Direction)) +
  geom_col(position = "dodge") +
  facet_wrap(~Genotype, nrow = 1) +
  scale_y_continuous(expand = c(0,0), limit = c(0, 2550)) +
  scale_fill_manual(values = c("Downregulated" = "#e28743", "Upregulated" = "#1e81b0")) +
  labs(y = "Number of DEGs (infected vs mock)", 
       fill = "DE Direction") +
  theme_bw() +
  theme(text = element_text(size = 18),
        legend.position = c(0.99, 0.99), 
        legend.justification = c("right", "top"), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12), 
        legend.background = element_rect(fill = "transparent")
bargraph

ggsave("../../04_figures/REDO_bcvsmock_sigDEGs_bargraph1.png")

