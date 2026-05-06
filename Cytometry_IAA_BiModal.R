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

setwd("/Users/benjamindowning/Desktop/Work/RawData/Cytometry/20260306_IAA")

flow_set <- read.plateSet(path = "20260306_IAA/", pattern = "20260306")

# Check how sample names look
sampleNames(flow_set)

# Clean up annotation names if needed
annotation$name <- gsub("_", "", annotation$name)
annotation$name <- gsub("-", "", annotation$name)
annotation$name <- trimws(annotation$name)  # remove leading/trailing spaces

annotation <- read.csv("20260306_annotation.csv")

annotation$name <- gsub(x = annotation$name, pattern = "_", replacement = "")

flow_set <- annotateFlowSet(flow_set, annotation)

write.flowSet(flow_set, outdir = "20260306_SS")

#write the annotation file
flow_set <- read.flowSet(path = "20260306_SS", phenoData = "annotation.txt")

#dat_sum combines the annotation with the actual data
dat_sum <- summarizeFlow(flow_set, channel = NA, ploidy = "haploid", only = "yeast")
dat_sum2 <- summarizeFlow(flow_set, channel = NA, ploidy = "haploid", only = "singlets")


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
                  normv <- (FL2.Amean/FSC.Amean)})
                
                
                dat_sum1 <- subset(dat_sum, strain != "BLANK")   
                
                ########03.12.26
                dat_sum_iaa2<- subset(dat_sum1, strain %in% c("4470", "5389", "5391", "5393", "5478-2", "5479-4", "5480-2", "5483-2", "5484-1", "5485-2"))
                
                # 1. Median + 95% CI function
                median_ci <- function(x) {
                  x <- x[is.finite(x)]
                  data.frame(
                    y    = median(x, na.rm = TRUE),
                    ymin = quantile(x, 0.025, na.rm = TRUE),
                    ymax = quantile(x, 0.975, na.rm = TRUE)
                  )
                }
                
                ##BImodal fit 
                plot_dat <- dat_sum_iaa2
                
                plot_dat$Alias <- factor(plot_dat$Alias, levels = desired_order)
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "TPL | IAA3" = 1.02,
                  "TPL | IAA3 (CG)" = 1.04,
                  "TPL | IAA3 (TA)" = 1.06,
                  "TPR4 | IAA3" = 1.08,
                  "TPR4 | IAA3 (CG)" = 1.10,
                  "TPR4 | IAA3 (TA)" = 1.12,
                  "TPR2/3 | IAA3" = 1.14,
                  "TPR2/3 | IAA3 (CG)" = 1.16,
                  "TPR2/3 | IAA3 (TA)" = 1.18
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Alias)]
                
                y_min <- 2000
                y_max <- 4000
                
                plotsc4 <- ggplot(
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
                    size = 7,
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
                
                plotsc4      
                
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
                    filename = "plotsc_bimode_1.pdf",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc4_bimode_1.png",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }
                
                
                ##BImodal fit part 2 
                plot_dat <- dat_sum_iaa2
                
                plot_dat$Alias <- factor(plot_dat$Alias, levels = desired_order)
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "TPL | IAA3" = 1.02,
                  "TPL | IAA3 (CG)" = 1.04,
                  "TPL | IAA3 (TA)" = 1.06,
                  "TPR4 | IAA3" = 1.08,
                  "TPR4 | IAA3 (CG)" = 1.10,
                  "TPR4 | IAA3 (TA)" = 1.12,
                  "TPR2/3 | IAA3" = 1.14,
                  "TPR2/3 | IAA3 (CG)" = 1.16,
                  "TPR2/3 | IAA3 (TA)" = 1.18
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Alias)]
                
                y_min <- 4000
                y_max <- 6000
                
                plotsc4 <- ggplot(
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
                    size = 7,
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
                
                plotsc4      
                
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
                    filename = "plotsc_bimode_2.pdf",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc4_bimode_2.png",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }
                
                
                ##BImodal fit w/ contorl
                plot_dat <- dat_sum_iaa2
                
                plot_dat$Alias <- factor(plot_dat$Alias, levels = desired_order)
                
                # manually choose x positions
                x_pos <- c(
                  "IAA3 (-)" = 1.00,
                  "TPL | IAA3" = 1.02,
                  "TPL | IAA3 (CG)" = 1.04,
                  "TPL | IAA3 (TA)" = 1.06,
                  "TPR4 | IAA3" = 1.08,
                  "TPR4 | IAA3 (CG)" = 1.10,
                  "TPR4 | IAA3 (TA)" = 1.12,
                  "TPR2/3 | IAA3" = 1.14,
                  "TPR2/3 | IAA3 (CG)" = 1.16,
                  "TPR2/3 | IAA3 (TA)" = 1.18
                )
                
                plot_dat$x_manual <- x_pos[as.character(plot_dat$Alias)]
                
                y_min <- 7000
                y_max <- 9000
                
                plotsc4 <- ggplot(
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
                    size = 7,
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
                
                plotsc4          
                
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
                    filename = "plotsc_bimode_control_1.pdf",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units
                  )
                  
                  ggsave(
                    filename = "plotsc4_bimode_control_1.png",
                    plot = plotsc4,
                    width = plot_width,
                    height = plot_height,
                    units = plot_units,
                    dpi = plot_dpi
                  )
                }