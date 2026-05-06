# Join modules to normalized counts
mut_counts <- modules %>%
  left_join(dat_counts, by = "Gene_ID")

# Pivot to long format
mut_longer <- mut_counts %>%
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
mut_scaled <- mut_longer %>%
  group_by(Gene_ID) %>%
  mutate(scaled_expression = (count - mean(count)) / sd(count)) %>%
  ungroup()

# Summarise per cluster/genotype/treatment/timepoint
mut_sum <- mut_scaled %>%
  group_by(cluster, genotype, treatment, timepoint) %>%
  summarise(
    mean = mean(scaled_expression, na.rm = TRUE),
    n    = n(),
    sd   = sd(scaled_expression, na.rm = TRUE),
    se   = sd / sqrt(n),
    ci   = 1.96 * se,
    .groups = "drop"
  )

# Add cluster labels with n sizes
cl_sizes_mut <- table(modules$cluster)
mut_sum <- mut_sum %>%
  mutate(cluster_name = paste0("Cluster ", cluster,
                               " (n = ", cl_sizes_mut[cluster], ")"),
         cluster_name = factor(cluster_name,
                               levels = paste0("Cluster ", 1:11,
                                               " (n = ", as.vector(cl_sizes_mut), ")")))

# Plot with facet_wrap in 2 rows
ggplot(mut_sum, aes(x = timepoint, y = mean,
                    colour = genotype, fill = genotype)) +
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
  facet_wrap(~ cluster_name, nrow = 2) +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi) - Mock samples") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype")

# Check for NAs
sum(is.na(mut_sum$mean))

# Check that all genotypes have both timepoints for all clusters
mut_sum %>%
  dplyr::count(cluster, genotype, treatment, timepoint) %>%
  dplyr::filter(n != 1)


ggplot(mut_sum, aes(x = timepoint, y = mean,
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
  facet_wrap(~ cluster_name, nrow = 2) +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi) - Mock samples") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype")


ggplot(mut_sum %>% dplyr::filter(treatment == "Mock"),
       aes(x = timepoint, y = mean,
           colour = genotype, fill = genotype,
           group = genotype)) +
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
  facet_wrap(~ cluster_name, nrow = 2) +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype")



pdf("line-graph_mockvsmuts.pdf", width = 7, height = 3)
ggplot(mut_sum %>% dplyr::filter(treatment == "Mock"),
       aes(x = timepoint, y = mean,
           colour = genotype, fill = genotype,
           group = genotype)) +
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
  facet_wrap(~ cluster_name, nrow = 2) +
  scale_x_continuous(breaks = seq(16, 20, 4)) +
  xlab("Hours post inoculation (hpi)") +
  ylab("Scaled Expression") +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    text             = element_text(size = 10),
    strip.background = element_blank(),
    panel.grid       = element_blank()
  ) +
  labs(color = "Genotype", fill = "Genotype")
dev.off()
