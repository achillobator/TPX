library(ggthemes)
library(ggplot2)
library("viridis") 
library(dplyr)
library(forcats)

dat <- read.csv("GAL4UAS_TC_LRD_2026203_trim.csv")

dat <- subset(dat, PRL >  40)

dat <- subset(dat, Genotype != '5913-24')
dat <- subset(dat, Genotype != '5939-1')
dat <- subset(dat, Genotype != '5913-20')


#This piece takes a data file and orders it based on the category after the money sign
dat$geno <- factor(dat$geno, levels=c("Wt","mCherry:RxL21","HA:RxL21"),ordered=TRUE)
B_colors <- c("black","red","blue")
#########
#This plots all the data initially to ensure the data has been loaded in and is using primary root lenght (RPL)
plot1<- ggplot(dat, aes(x=geno, y = PRL, color = Genotype, group = geno)) + 
  geom_boxplot(lwd=1, outlier.size = 0, alpha = 0.7) +
  geom_jitter(size =2, width = 0.1, alpha = 0.7) +
  ylab('Root length (mm)') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_size = 10)+
  #scale_colour_manual(values = B_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position="none")
plot1


plot2<- ggplot(dat, aes(x=geno, y = LRD, color = geno, group = geno)) + 
  geom_boxplot(lwd=1, outlier.size = 0, alpha = 0.7) +
  geom_jitter(size =2, width = 0.1, alpha = 0.7) +
  ylab('Lateral Root Density') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_size = 14)+
  scale_colour_manual(values = B_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position="none")
plot2
plot2+facet_wrap(~day, nrow=1)

plot2<- ggplot(dat, aes(x=Genotype, y = LRD, color = Genotype)) + 
  geom_boxplot(lwd=1, outlier.size = 0, alpha = 0.7) +
  geom_jitter(size =2, width = 0.1, alpha = 0.7) +
  ylab('Lateral Root Density') +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_size = 14)+
  #scale_colour_manual(values = A_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position="none")
plot2+facet_wrap(~day, nrow=1)
plot2+facet_grid(geno~day)

dat_13 <- subset(dat, day == 13)

plot13<- ggplot(dat_13, aes(x=geno, y = LRD, color = geno)) + 
  geom_boxplot(lwd=1, outlier.size = 0, alpha = 0.7) +
  geom_jitter(size =2, width = 0.1, alpha = 0.7) +
  ylab('Lateral Root Density') +
  xlab("Genotype") + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) +
  theme_classic(base_size = 14)+
  scale_colour_manual(values = A_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.position="none")
plot13

plot13 <- ggplot(dat_13, aes(x=geno, y = LRD, color = geno)) + 
  geom_boxplot(lwd=1, outlier.size = 0, alpha = 0.7) +
  geom_jitter(size =2, width = 0.1, alpha = 0.7) +
  ylab('Lateral Root Density') +
  xlab("Genotype") + 
  theme_classic(base_size = 14)+
  scale_colour_manual(values = A_colors)+
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position="none") +
  annotate("text", x = 1, y = max(dat_13$LRD) * 1.05, label = "a", size = 5) +
  annotate("text", x = 2, y = max(dat_13$LRD) * 1.05, label = "b", size = 5) +
  annotate("text", x = 3, y = max(dat_13$LRD) * 1.05, label = "b", size = 5)
plot13


B_colors <- c("black","red","blue")
jitter_colors <- c(
  "Wt" = "#000000",        # black
  
  # 6 shades of red (light to dark)
  "#FF9999",
  "#FF6666",
  "#FF3333",
  "#CC0000",
  "#990000",
  "#660000",
  
  # 6 shades of blue (light to dark)
  "#99BBFF",
  "#6699FF",
  "#3366FF",
  "#0033CC",
  "#002299",
  "#001166"
)

install.packages("ggnewscale")
library(ggnewscale)

plot13b <- ggplot(dat_13, aes(x=geno, y = LRD)) + 
  geom_boxplot(aes(color = geno), lwd=1, outlier.size = 0, alpha = 0.7) +
  scale_colour_manual(values = B_colors) +  # boxplot colors
  ggnewscale::new_scale_color() +           # reset color scale
  geom_jitter(aes(color = Genotype), size=2, width = 0.1, alpha = 0.7) +
  scale_colour_manual(values = jitter_colors)+
  ylab('Lateral Root Density') +
  xlab("Genotype") + 
  theme_classic(base_size = 14)+
  annotate("text", x = 1, y = max(dat_13$LRD) * 1.05, label = "a", size = 5) +
  annotate("text", x = 2, y = max(dat_13$LRD) * 1.05, label = "b", size = 5) +
  annotate("text", x = 3, y = max(dat_13$LRD) * 1.05, label = "b", size = 5)

plot13b

B_colors <- c("Wt" = "black", "mCherry:RxL21" = "red", "HA:RxL21" = "blue")  # your geno group colors

jitter_colors <- c(
  "5057"    = "#000000",
  "5913-17" = "#FF6666",
  "5913-18" = "#FF3333",
  "5913-21" = "#990000",
  "5913-23" = "#660000",
  "5913-4"  = "#FF9999",
  "5921-2"  = "#3366FF",
  "5921-3"  = "#6699FF",
  "5921-4"  = "#0033CC",
  "5921-5"  = "#001166",
  "5938-3"  = "#CC0000",
  "5939-2"  = "#002299",
  "5939-3"  = "#001166"
)
####

library(dplyr)
dat_E<- group_by(dat, Genotype, day, geno) %>% 
  summarise(
    count = n(), 
    mean = mean(LRD, na.rm = TRUE),
    sd = sd(LRD, na.rm = TRUE),
    ymin = (mean-sd),
    ymax = (mean+sd),
    count = n()
  )

plot_time <- ggplot(data = dat_E, mapping = aes(x = day, y = mean, color=geno)) +
  geom_point(size=.5) + 
  geom_line() + 
  ylab("Mean LRD") + 
  xlab("Time (Days)") + 
  labs(color = "Genotype", title = "LRD TC") + 
  theme_classic(base_family = 'Arial Bold', base_size = 12) +
  geom_errorbar(aes(ymin=mean-(sd/sqrt(count)), ymax=mean+(sd/sqrt(count))), size = 0.2, width = 0.05) + 
  theme(legend.position="none")+
  scale_colour_manual(values = A_colors)+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background  = element_blank()) 
plot_time
plot_time+facet_wrap(~geno, nrow = 5)

plot_time <- ggplot(data = dat_E, mapping = aes(x = day, y = mean, color = geno, group = Genotype)) +
  geom_point(size=.5) + 
  geom_line(linewidth = 0.8) + 
  ylab("Mean LRD") + 
  xlab("Time (Days)") + 
  labs(color = "Genotype", title = "LRD TC") + 
  theme_classic(base_family = 'Arial Bold', base_size = 12) +
  geom_errorbar(aes(ymin=mean-(sd/sqrt(count)), ymax=mean+(sd/sqrt(count))), size = 0.3, width = 0.05) + 
  theme(legend.position="none")+
  scale_colour_manual(values = A_colors)+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background = element_blank())
plot_time

#########

# Pull out just the Wt rows
wt_rows <- dat_E %>% filter(geno == "Wt")  # or however Wt is coded in geno

# Create a copy for each non-Wt geno group
wt_duplicated <- bind_rows(
  wt_rows %>% mutate(geno = "geno1"),  # replace with your actual geno names
  wt_rows %>% mutate(geno = "geno2")
)

# Bind back to the original data (excluding original Wt if you want, or keep it)
dat_E_plot <- bind_rows(dat_E, wt_duplicated)

A_colors <- c("Wt" = "black", "mCherry:RxL21" = "red", "HA:RxL21" = "blue")

plot_time <- ggplot(data = dat_E_plot, mapping = aes(x = day, y = mean, color = color_group, group = Genotype)) +
  geom_point(size=.5) + 
  geom_line(linewidth = 0.8) + 
  ylab("Mean LRD") + 
  xlab("Time (Days)") + 
  labs(color = "Genotype", title = "UAS LRD TC") + 
  theme_classic(base_size = 12) +
  geom_errorbar(aes(ymin=mean-(sd/sqrt(count)), ymax=mean+(sd/sqrt(count))), size = 0.3, width = 0.05) + 
  theme(legend.position="none")+
  scale_colour_manual(values = A_colors)+
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), strip.background = element_blank())
plot_time+facet_wrap(~geno, nrow = 1)

plot_time + facet_wrap(~geno, nrow = 1)
plot_time + facet_wrap(~Genotype, nrow = 1)

wt_duplicated <- do.call(bind_rows, lapply(non_wt_genos, function(x) wt_rows %>% mutate(geno = x, color_group = "Wt")))

dat_E_plot <- bind_rows(
  dat_E %>% filter(geno != "Wt") %>% mutate(color_group = geno),
  wt_duplicated
)

# Check what values are actually in geno
unique(dat_E_plot$geno)

# Check what names are in your color vector
names(A_colors)
## old below here


##########
dat7 <- dat
dat7 <- subset(dat, day == 13)
#dat_7 <- subset(dat7, comment !=  "fungus")
#write.csv(dat7, "dat7.csv")

#Run ANNOVA
attach(dat7)
data(dat7)
str(dat7$dat7)
tapply(dat7$geno, dat7$LRD, dat7$mean) 
tapply(dat7$geno, dat7$LRD, dat7$var)
tapply(dat7$geno, dat7$LRD, dat7$length_mm)
boxplot(LRD ~ geno, data=dat7)

lm.out = with(dat7, lm(LRD ~ geno))
aov.out = aov(LRD ~ geno, data=dat7)
oneway.test(LRD ~ geno, data=dat7)
is.factor(dat7$LRD)
is.factor(dat7$geno)

aov.out
summary(aov.out)
TukeyHSD(aov.out)

summary.lm(aov.out)

model=lm( dat7$LRD ~ dat7$geno )
ANOVA=aov(model)
TUKEY <- TukeyHSD(x=ANOVA, 'dat7$construct', conf.level=0.95)
plot(TUKEY , las=1 , col="brown")

library(multcompView)
install.packages('rcompanion')
library(rcompanion)
library(lsmeans)
library(multcomp)

marginal = lsmeans(model, ~ geno)
pairs(marginal, adjust="tukey")
CLD = cld(marginal, alpha = 0.05, Letters = letters, adjust  = "tukey")
CLD

#
#
#



