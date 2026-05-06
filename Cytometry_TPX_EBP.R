
library(dplyr)
library(ggthemes)
library(ggplot2)
library(flowCore)
library(flowTime)
library(ggpubr)
library(RColorBrewer)
library(patchwork)

#Read in the data using flowtime and write the annotation CSV
read.plateSet <- function(path = getwd(), pattern = "\\."){ 
  files <- grep(pattern = pattern, x = list.files(path), value = TRUE)
  for(file in files) {
    plate <- read.flowSet(path = paste(path, file, sep = "/"), alter.names = TRUE)
    plate_num <- which(files == file) 
    sampleNames(plate) <- paste0(plate_num, sampleNames(plate))
    pData(plate)$name <- sampleNames(plate)
    pData(plate)$folder <- file
    if('flow_set' %in% ls()) flow_set <- rbind2(flow_set, plate) else flow_set <- plate
  }
  return(flow_set)
}

setwd("/Users/benjamindowning/Desktop/Work/RawData/Cytometry/20250912_KLD")

flow_set <- read.plateSet(path = "20250912_162202/", pattern = "20250912")

# Check how sample names look
sampleNames(flow_set)

# Clean up annotation names if needed
annotation$name <- gsub("_", "", annotation$name)
annotation$name <- gsub("-", "", annotation$name)
annotation$name <- trimws(annotation$name)  # remove leading/trailing spaces

annotation <- read.csv("20250912_annotation.csv")

annotation$name <- gsub(x = annotation$name, pattern = "_", replacement = "")

flow_set <- annotateFlowSet(flow_set, annotation)

write.flowSet(flow_set, outdir = "20240912_KLD")

#write the annotation file
flow_set <- read.flowSet(path = "20240912_KLD", phenoData = "annotation.txt")

#dat_sum combines the annotation with the actual data
dat_sum <- summarizeFlow(flow_set, channel = NA, ploidy = "haploid", only = "yeast")
dat_sum2 <- summarizeFlow(flow_set, channel = NA, ploidy = "haploid", only = "singlets")

#nice data table :)
dat_sum

median_cl_boot <- function(x, conf = 0.95) {
  lconf <- (1 - conf)/2
  uconf <- 1 - lconf
  require(boot)
  bmedian <- function(x, ind) median(x[ind])
  bt <- boot(x, bmedian, 1000)
  bb <- boot.ci(bt, type = "perc")
  data.frame(y = median(x), ymin = quantile(bt$t, lconf), ymax = quantile(bt$t, uconf))
}
                
                dat_sum <- within(dat_sum, { 
                  normv <- (FL2.Amean/FSC.Amean)
                })
                
                dat_sum1 <- subset(dat_sum, strain != "BLANK")
              
                
                # LARGE CHECK Plot from dat_sum1
                plot48 <- ggplot(dat_sum1, aes(x = strain, y = FL2.Amean, fill = strain)) +
                  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # suppress outliers to avoid double plotting
                  geom_jitter(aes(color = strain), width = 0.2, size = 1, alpha = 0.8) +  # show replicates
                  facet_wrap(~Facet, ncol = 4) +
                  labs(title = "9.12.25", subtitle = "5 + 10 µL", y = "FL2.Amean", x = "strain") +
                  theme_classic(base_size = 10) +
                  theme(plot.title = element_text(hjust = 0.5),
                        plot.subtitle = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle = 45, hjust = 1)) +
                  scale_fill_discrete(guide = "none") +
                  scale_color_discrete(guide = "none")
                plot48
                
				####### 03.12.2026
                
                
                
                dat_sum_f2.3 <- subset(dat_sum1, strain %in% c("4496-2", "4494-2", "4497-3", "4495-3"))
                dat_sum_f2.3_control <- subset(dat_sum1, strain %in% c("4470","4496-2", "4494-2", "4497-3", "4495-3"))
                
                dat_sum_f4 <- subset(dat_sum1, strain %in% c("4496-2","4492-4","4497-3","4493-1"))
                dat_sum_f4_control<- subset(dat_sum1, Sample %in% c("4470","4496-2","4492-4","4497-3","4493-1"))
                
                # median + 95% CI (quantile-based)
                median_ci <- function(x) {
                  x <- x[is.finite(x)]
                  data.frame(
                    y    = median(x, na.rm = TRUE),
                    ymin = quantile(x, 0.025, na.rm = TRUE),
                    ymax = quantile(x, 0.975, na.rm = TRUE)
                  )
                }
                
                desired_order_f23 <- c(
                  "2xHA : TPL : IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR2/3] : IAA3",
                  "2xHA : TPL | IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR2/3] | IAA3"
                )
                
                desired_order_f4 <- c(
                  "2xHA : TPL : IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR4] : IAA3",
                  "2xHA : TPL | IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR4] | IAA3"
                )
                
                ###2/3 samples
                plot_dat <- dat_sum_f2.3_control
                
                plot_dat$Sample <- factor(plot_dat$Sample, levels = desired_order)
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "2xHA : TPL : IAA3 (+) (FL)" = 1.02,
                  "2xHA : TPL [EBP_TPR2/3] : IAA3" = 1.04,
                  "2xHA : TPL | IAA3 (+) (FL)" = 1.06,
                  "2xHA : TPL [EBP_TPR2/3] | IAA3" = 1.08
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Sample)]
                
                # manually set y-axis here
                y_min <- 2000
                y_max <- 4500
                
                plotsc23.1 <- ggplot(
                  plot_dat,
                  aes(x = x_manual, y = FL2.Amean)
                ) +
                  stat_summary(
                    fun.data = median_ci,
                    geom = "errorbar",
                    width = 0.015,
                    linewidth = 1.2,
                    color = "#999999ff"
                  ) +
                  stat_summary(
                    fun = median,
                    geom = "point",
                    shape = 18,
                    size = 9,
                    stroke = 0.5,
                    color = "#000000ff"
                  ) +
                  scale_x_continuous(
                    breaks = unname(x_pos),
                    labels = names(x_pos)
                  ) +
                  coord_cartesian(ylim = c(y_min, y_max)) +
                  labs(
                    title = ".",
                    subtitle = ".",
                    y = ".",
                    x = "."
                  ) +
                  theme_classic(base_size = 12) +
                  theme(axis.text.x = element_text(angle = 0, hjust = 1))
                
                plotsc23.1
                
                # ===== optional export block: turn ON/OFF here =====
                save_plot <- TRUE
                
                # set dimensions manually here
                plot_width  <- 3
                plot_height <- 4
                
                # choose units: "in", "cm", or "mm"
                plot_units <- "in"
                
                # png resolution
                plot_dpi <- 600
                
                if (save_plot) {
                  ggsave(
                    filename = "plotsc23_1.pdf",
                    plot = plotsc23.1,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc23_1.png",
                    plot = plotsc23.1,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }
                
                
                ##2/3 w/ Control
                plot_dat <- dat_sum_f2.3_control
                
                plot_dat$Sample <- factor(plot_dat$Sample, levels = desired_order)
                
                # manually set y-axis here
                y_min <- 7000
                y_max <- 9500
                
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "2xHA : TPL : IAA3 (+) (FL)" = 1.02,
                  "2xHA : TPL [EBP_TPR2/3] : IAA3" = 1.04,
                  "2xHA : TPL | IAA3 (+) (FL)" = 1.06,
                  "2xHA : TPL [EBP_TPR2/3] | IAA3" = 1.08
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Sample)]
                
                plotsc23.1 <- ggplot(
                  plot_dat,
                  aes(x = x_manual, y = FL2.Amean)
                ) +
                  stat_summary(
                    fun.data = median_ci,
                    geom = "errorbar",
                    width = 0.015,
                    linewidth = 1.2,
                    color = "#999999ff"
                  ) +
                  stat_summary(
                    fun = median,
                    geom = "point",
                    shape = 18,
                    size = 9,
                    stroke = 0.5,
                    color = "#000000ff"
                  ) +
                  scale_x_continuous(
                    breaks = unname(x_pos),
                    labels = names(x_pos)
                  ) +
                  coord_cartesian(ylim = c(y_min, y_max)) +
                  labs(
                    title = ".",
                    subtitle = ".",
                    y = ".",
                    x = "."
                  ) +
                  theme_classic(base_size = 12) +
                  theme(axis.text.x = element_text(angle = 0, hjust = 1))
                
                plotsc23.1
                
                # ===== optional export block: turn ON/OFF here =====
                save_plot <- TRUE
                
                # set dimensions manually here
                plot_width  <- 3
                plot_height <- 4
                
                # choose units: "in", "cm", or "mm"
                plot_units <- "in"
                
                # png resolution
                plot_dpi <- 600
                
                if (save_plot) {
                  ggsave(
                    filename = "plotsc23_control_1.pdf",
                    plot = plotsc23.1,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc23_control_1.png",
                    plot = plotsc23.1,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }
            
          ###########
                ###4 samples
               
                dat_sum_f2.3 <- subset(dat_sum1, strain %in% c("4496-2", "4494-2", "4497-3", "4495-3"))
                dat_sum_f2.3_control <- subset(dat_sum1, strain %in% c("4470","4496-2", "4494-2", "4497-3", "4495-3"))
                
                dat_sum_f4 <- subset(dat_sum1, strain %in% c("4496-2","4492-4","4497-3","4493-1"))
                dat_sum_f4_control<- subset(dat_sum1, strain %in% c("4470","4496-2","4492-4","4497-3","4493-1"))
                
                # median + 95% CI (quantile-based)
                median_ci <- function(x) {
                  x <- x[is.finite(x)]
                  data.frame(
                    y    = median(x, na.rm = TRUE),
                    ymin = quantile(x, 0.025, na.rm = TRUE),
                    ymax = quantile(x, 0.975, na.rm = TRUE)
                  )
                }
                
                desired_order_f23 <- c(
                  "2xHA : TPL : IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR2/3] : IAA3",
                  "2xHA : TPL | IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR2/3] | IAA3"
                )
                
                desired_order_f4 <- c(
                  "2xHA : TPL : IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR4] : IAA3",
                  "2xHA : TPL | IAA3 (+) (FL)",
                  "2xHA : TPL [EBP_TPR4] | IAA3"
                )
              
                plot_dat <- dat_sum_f4_control
                
                plot_dat$Sample <- factor(plot_dat$Sample, levels = desired_order)
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "2xHA : TPL : IAA3 (+) (FL)" = 1.02,
                  "2xHA : TPL [EBP_TPR4] : IAA3" = 1.04,
                  "2xHA : TPL | IAA3 (+) (FL)" = 1.06,
                  "2xHA : TPL [EBP_TPR4] | IAA3" = 1.08
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Sample)]
                
                # manually set y-axis here
                y_min <- 2000
                y_max <- 4500
                
                plotsc4.5 <- ggplot(
                  plot_dat,
                  aes(x = x_manual, y = FL2.Amean)
                ) +
                  stat_summary(
                    fun.data = median_ci,
                    geom = "errorbar",
                    width = 0.015,
                    linewidth = 1.2,
                    color = "#999999ff"
                  ) +
                  stat_summary(
                    fun = median,
                    geom = "point",
                    shape = 18,
                    size = 9,
                    stroke = 0.5,
                    color = "#000000ff"
                  ) +
                  scale_x_continuous(
                    breaks = unname(x_pos),
                    labels = names(x_pos)
                  ) +
                  coord_cartesian(ylim = c(y_min, y_max)) +
                  labs(
                    title = ".",
                    subtitle = ".",
                    y = ".",
                    x = "."
                  ) +
                  theme_classic(base_size = 12) +
                  theme(axis.text.x = element_text(angle = 0, hjust = 1))
                
                plotsc4.5
                
                # ===== optional export block: turn ON/OFF here =====
                save_plot <- TRUE
                
                # set dimensions manually here
                plot_width  <- 3
                plot_height <- 4
                
                # choose units: "in", "cm", or "mm"
                plot_units <- "in"
                
                # png resolution
                plot_dpi <- 600
                
                if (save_plot) {
                  ggsave(
                    filename = "plotsc4_1.pdf",
                    plot = plotsc4.5,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc4_1.png",
                    plot = plotsc4.5,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }
                
                
                ##4  w/ Control
                plot_dat <- dat_sum_f4_control
                
                plot_dat$Sample <- factor(plot_dat$Sample, levels = desired_order)
                
                # manually set y-axis here
                y_min <- 7000
                y_max <- 9500
                
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "2xHA : TPL : IAA3 (+) (FL)" = 1.02,
                  "2xHA : TPL [EBP_TPR4] : IAA3" = 1.04,
                  "2xHA : TPL | IAA3 (+) (FL)" = 1.06,
                  "2xHA : TPL [EBP_TPR4] | IAA3" = 1.08
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Sample)]
                
                plotsc4.6 <- ggplot(
                  plot_dat,
                  aes(x = x_manual, y = FL2.Amean)
                ) +
                  stat_summary(
                    fun.data = median_ci,
                    geom = "errorbar",
                    width = 0.015,
                    linewidth = 1.2,
                    color = "#999999ff"
                  ) +
                  stat_summary(
                    fun = median,
                    geom = "point",
                    shape = 18,
                    size = 9,
                    stroke = 0.5,
                    color = "#000000ff"
                  ) +
                  scale_x_continuous(
                    breaks = unname(x_pos),
                    labels = names(x_pos)
                  ) +
                  coord_cartesian(ylim = c(y_min, y_max)) +
                  labs(
                    title = ".",
                    subtitle = ".",
                    y = ".",
                    x = "."
                  ) +
                  theme_classic(base_size = 12) +
                  theme(axis.text.x = element_text(angle = 0, hjust = 1))
                
                plotsc4.6
                
                # ===== optional export block: turn ON/OFF here =====
                save_plot <- TRUE
                
                # set dimensions manually here
                plot_width  <- 3
                plot_height <- 4
                
                # choose units: "in", "cm", or "mm"
                plot_units <- "in"
                
                # png resolution
                plot_dpi <- 600
                
                if (save_plot) {
                  ggsave(
                    filename = "plotsc4_control_1.pdf",
                    plot = plotsc4.6,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc4_control_1.png",
                    plot = plotsc4.6,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }