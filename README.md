# TPX
Conservation and Divergence in TPX family - Code Repository

#Representative Code for flow cytometry of auxin response circuits in yeast (Saccharomyces cereviseae): 

>Cytometry_IAA_BiModal.R
This script analyzes flow cytometry data from IAA3-related yeast mutants. It loads FCS files from plates into a FlowSet, reads in annotation metadata with sample information, and summarizes the data for singlet analysis. The script normalizes FL2 fluorescence (a proxy for protein expression) and filters to specific strain genotypes. It creates three publication-quality plots comparing median FL2.Amean levels across strains expressing IAA3 with various TPL and TPR domain constructs, each spanning different y-axis ranges (suggesting bimodal populations at different expression levels). Error bars show 95% confidence intervals, and results export to both PDF and PNG.
>
>Cytometry_TPX_EBP.R
This script follows a similar workflow but analyzes yeast strains expressing TPL with engineered EBP (EBP_TPR2/3 or EBP_TPR4) substitutions. After loading and annotating flow data, it generates an initial faceted boxplot to assess overall data quality across experimental conditions. It then creates subset-specific analyses: one comparing TPR2/3 domain variants and another for TPR4 domain variants. For each, it produces two plots per domain (low and high expression ranges) showing how EBP domain substitutions affect FL2 signal in the presence of IAA3. Like the IAA3 script, it uses median+CI visualization and exports to PDF/PNG.
Common workflow: Both scripts follow a reproducible template—load flow data → annotate → normalize → subset by strain → compute statistics → plot medians with 95% CI → export figures. The code is exploratory-style (direct commands) rather than wrapped in functions.

#Representative Code for Lateral root density analysis in Arabidopsis transgenics:

>Graphing_LRD_UAS_driver.R
This script analyzes lateral root density (LRD) phenotypes in Arabidopsis transgenic lines expressing GAL4/UAS-driven constructs (mCherry:RxL21 and HA:RxL21 fusions). It loads a CSV file of root trait measurements, filters for viable primary root lengths (PRL > 40 mm), and removes contaminated samples. The script creates multiple visualizations: initial boxplots assessing data quality, day 13-specific LRD comparisons with genotype-specific colors (using ggnewscale to layer boxplot and jitter point colors independently), and time-course line plots tracking mean LRD across days with standard error bars. For statistical analysis, it performs one-way ANOVA on day 13 data followed by Tukey HSD post-hoc testing and compact letter display (CLD) to identify and annotate significant pairwise differences between genotypes. Outputs are organized by genotype category (Wt, mCherry, HA) and faceted by experimental timepoint.
