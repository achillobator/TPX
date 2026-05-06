## redo bcvsmock clusters scaled plots
## plotting clusters in line graph 
library(tidyverse)

setwd("03_results/14_REDO_bcvsmock_tps/")

upreg_degs <- read.csv("bcvsmock_upreg_unconserved_sigdegs.csv")
downreg_degs <- read.csv("bcvsmock_downreg_unconserved_sigdegs.csv")

upreg_degs_genes <- select(upreg_degs, Gene_ID)
downreg_degs_genes <- select(downreg_degs, Gene_ID)
all_uncond_degs  <-  bind_rows(upreg_degs_genes, downreg_degs_genes) %>% 
  unique()

norm_counts <- read.csv("../../01_data/USETHIS_At_tplmuts_mean_norm_counts.csv")

degs_counts <- all_uncond_degs %>% left_join(norm_counts, by = join_by(Gene_ID)) %>% select(-2)

degs_longer <- degs_counts %>% 
  pivot_longer(cols = -1,
               names_to = c("genotype", "condition"), 
               names_sep = "_", 
               values_to = "count")
degs_longer <- degs_longer %>% 
  separate(col = condition, into = c("treatment", "timepoint"), sep = 1)

## getting clusters 
cl <- read.csv("bcvsmock_T16T20_HUGE_pheatmap_clusters_unconserved.csv")
cl <- cl %>% rename(Gene_ID = X, cl = cl)

#####trying to make a freaking loop 
cl_nums <- c("1", "2", "3", "4", "5", "6")

cl_num_lists <- lapply(seq_along(cl_nums), function(n) filter(cl, cl == n))
names(cl_num_lists) <- cl_nums

degs_longer$timepoint <- as.numeric(degs_longer$timepoint)

cl_num_lists_counts <- lapply(seq_along(cl_num_lists), function(x) left_join(cl_num_lists[[x]], degs_longer, by = join_by(Gene_ID)))
names(cl_num_lists_counts) <- cl_nums

cl_num_lists_scaled <- lapply(seq_along(cl_num_lists_counts), 
                              function(s) cl_num_lists_counts[[s]] %>% 
                                group_by(Gene_ID) %>% 
                                mutate(scaled_expression = (count - mean(count))/sd(count)) %>% 
                                ungroup())

cl_num_lists_sum <- lapply(seq_along(cl_num_lists_scaled),
                           function(k) cl_num_lists_scaled[[k]] %>% 
                             group_by(genotype, treatment, timepoint) %>% 
                             summarise(
                               mean = mean(scaled_expression),
                               n = length(scaled_expression),
                               sd = sd(scaled_expression),
                               se = sd/sqrt(n), 
                               ci = 1.96*se #95% confidence interval
                             ) %>% 
                             ungroup())

k_names <- c("Cluster 1 (n = 411)", "Cluster 2 (n = 1574)", "Cluster 3 (n = 317)", "Cluster 4 (n = 110)", "Cluster 5 (n = 785)", "Cluster 6 (n = 1154)")

names(cl_num_lists_sum) <- k_names

### saving all the lists 
# defining the file path to save all my lists
file_path <- "T16T20_clusters_scaled_expression/"

# save output function 
for(k in seq_along(cl_num_lists_sum)) {
  df_name <- names(cl_num_lists_sum)[k]
  file_name <- paste0(file_path, df_name, ".csv")
  write.csv(cl_num_lists_sum[[k]], file_name)
}

## now saving plots
image_path <- "cluster_scaled_exp_plots/"

for(m in seq_along(cl_num_lists_sum)) {
  df_name <- names(cl_num_lists_sum)[m]
  
  ggplot( data = cl_num_lists_sum[[m]], aes(x=timepoint, y=mean, colour = treatment, fill = treatment)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci), alpha = 0.2, colour = NA) +
    theme_bw() + 
    scale_x_continuous(breaks = seq(16,20,4))+
    xlab("Hours post inoculation (hpi)") +
    ylab("Scaled Expression") +
    facet_wrap(~genotype) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
          text = element_text(size = 18)) +
    labs(title = names(cl_num_lists_sum)[m], 
         color = "Treatment",
         fill = "Treatment")
  
  file_name <- paste0(image_path, df_name, "_1.png")
  ggsave(file_name)
}

## can i make a plot with facet wrap and ALL the clusters?
cl_all <- lapply(seq_along(cl_num_lists_sum), function(x) 
  cl_num_lists_sum[[x]] %>% 
    mutate(cluster_name = k_names[x]))
names(cl_all) <- k_names

cl_rowbind <- list_rbind(cl_all)

cl_rowbind <- cl_rowbind %>% 
  mutate(treatment = case_when(
    treatment == "B" ~ "Botrytis",
    treatment == "M" ~ "Mock"
  ))

ggplot( data = cl_rowbind, aes(x=timepoint, y=mean, colour = treatment, fill = treatment)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin=mean-ci, ymax=mean+ci), alpha = 0.3, color = NA) +
  theme_bw() + 
  scale_x_continuous(breaks = seq(16,20,4))+
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  facet_grid(cluster_name ~ genotype, 
             switch = "y") +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5),
        text = element_text(size = 16),
        strip.placement = "outside") +
  labs(color = "Treatment",
       fill = "Treatment")

ggsave("cluster_scaled_exp_plots/bigplot3.png", width = 7.14, height = 12, unit = "in", dpi = 1000)
