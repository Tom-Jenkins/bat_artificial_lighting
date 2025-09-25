# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# Genotype-Environment Association
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(LEA)
library(vegan)
library(qvalue)
library(robust)
library(data.table)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(forcats)

# Load SNP genlight data
load("data/WGS/snps_genlight.RData")

# K = 3 informed by the population structure analysis
K = 3

# Expected false discovery rate
alpha <- 0.05

# ---------- #
# Imputation
# ---------- #

# If true run imputation
runImputation = FALSE
if (runImputation) {
  
  # Perform snmf analysis
  set.seed(123)
  project <- snmf(
    input.file = "data/WGS/snps.lfmm",
    K = K,
    project = "new",
    repetitions = 10,
    entropy = TRUE,
    ploidy = 2
  )
  
  # Select the run with the lowest cross-entropy value
  best = which.min(cross.entropy(project, K = K))
  
  # Impute missing genotypes in all SNPs dataset
  LEA::impute(
    object = project,
    input.file = "data/WGS/snps.lfmm",
    method = "mode",
    K = K,
    run = best
  )
}

# ---------- #
# Data Preparation
# ---------- #

# Read in genotypes in LFMM format
geno <- data.table::fread("data/WGS/snps.lfmm_imputed.lfmm")

# Read in metadata
metadata <- read.csv("outputs/WGS_Rhipposideros_heterozygosity_ind.csv")
head(metadata)

# Subset environmental variables
envdata <- c("nightlight", "urban")
envdata <- dplyr::select(metadata, contains(envdata))
head(envdata)

# ---------- #
# RDA
# ---------- #

# Control for population structure between England and North Wales
popstruct = ifelse(metadata$Site == "NW1" | metadata$Site == "NW2", "NWales", "Remaining")

# Test multicollinearity
lm(nightlight ~ urban, data = envdata) |> summary()

# Run a PCA on nightlight and urban values to avoid multicollinearity
pca1 <- as.data.frame(scores(pca(envdata))$sites)

# Run partial RDA
# https://popgen.nescent.org/2018-03-27_RDA_GEA.html
rda1 <- rda(geno ~ . + Condition(popstruct), data = pca1)
rda1
RsquareAdj(rda1)

# Variance Inflation Factors
vif.cca(rda1)

# Variance explained by each canonical axis
summary(eigenvals(rda1, model = "constrained"))

# Screeplot (same number of axes as the number of predictors in the model)
screeplot(rda1)

# Check each canonical axis for significance (**very long run time ~30 mins)
# signif_axis <- anova.cca(rda1, by = "axis", parallel = 4)
# signif_axis

# Extract SNP loadings on constrained axes
snp_loadings <- scores(rda1, choices = 1, display = "species")

# Histogram of the loadings on each RDA axis
hist(snp_loadings, main="Loadings on RDA1")

# Function from Capblancq and Forester 2021 paper which returns q-values from RDA
# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R
# https://doi.org/10.1111/2041-210X.13722
rdadapt <- function(rda, K)
{
  zscores <- rda$CCA$v[, 1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- robust::covRob(resscale, distance = TRUE, na.action = na.omit, estim = "pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5, df=K)
  reschi2test <- pchisq(resmaha/lambda, K, lower.tail=FALSE)
  qval <- qvalue::qvalue(reschi2test)
  q.values_rdadapt <- qval$qvalues
  return(data.frame(p.values = reschi2test, q.values = q.values_rdadapt))
}
rda_pvals <- rdadapt(rda1, K = 2)

# P-values threshold after Bonferroni correction
thres_env <- 0.01/length(rda_pvals$p.values)

# Identify GEA candidate outlier SNPs
rda_outliers_idx <- which(rda_pvals$p.values < thres_env)
length(rda_outliers_idx)

# # Function where x is the vector of loadings and z is the number of standard deviations to use
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x)    
#   x[x < lims[1] | x > lims[2]]
# }

# # Identify candidate outlier SNPs 
# rda_outliers <- outliers(snp_loadings[,1], 3)

# # Number of candidate outliers
# rda_outliers_idx <- as.integer(str_remove(names(rda_outliers), "V"))
# length(rda_outliers_idx)

# # Find unique outliers identified by both methods
# all_outliers <- unique(c(lfmm_outliers_idx,rda1_outliers_idx))
# length(all_outliers)

# # Find common outliers
# common_outliers <- intersect(lfmm_outliers_idx, rda1_outliers_idx)
# length(common_outliers)

# # Find nightlight outliers found in both LFMM and RDA
# (common_nightlight <- lfmm_outliers |> 
#   filter(nightlight < alpha) |> 
#   pull() |> 
#   intersect(common_outliers)
# )

# All nightlight outliers
# all_nightlight <- lfmm_outliers[which(lfmm_outliers$nightlight < alpha),]

# ---------- #
# Outliers Within Genes
# ---------- #

# Print the contig ID and position of SNPs
snps[,rda_outliers_idx]$chromosome
snps[,rda_outliers_idx]$position

# Read in annotation bed file
annotation <- fread("../Lesser_horseshoe_bat_annotation/geneAnnotation.gtf.gz")
annotation <- annotation[, c(1,3,4,5,9)]
colnames(annotation) <- c("Chromosome","Type","Start_position","End_position","Gene_ID")
head(annotation)

# Function to locate if outlier indexes are within gene regions: returns outlier index or FALSE
locate_snps_within_genes <- function(genlight, annotation, outlier_idx) {
  # Extract chromosome ID and SNP position
  outlier_snp <- genlight[ , outlier_idx]
  outlier_snp_chrom <- str_remove(as.character(outlier_snp$chromosome), "\\.1$")
  outlier_snp_pos <- as.numeric(outlier_snp$position)

  # Filter annotation by chromosome
  outlier_chrom <- filter(annotation, annotation[["Chromosome"]] == outlier_snp_chrom)

  # Find within-gene indices (+- 1000bp)
  within_gene_idx <- outlier_chrom[
    (outlier_snp_pos >= Start_position - 1000) &
    (outlier_snp_pos <= End_position + 1000),
    which = TRUE
  ]

  # Return outlier index if true otherwise return false
  if (length(within_gene_idx) > 0){
    return(outlier_idx)
  } else {
    return(NA)
  }
}

# Test function
# locate_snps_within_genes(snps, annotation, 67280)

# Find the outlier indexes that are within gene regions (remove NAs which are not within gene regions)
idx_within_gene <- sapply(rda_outliers_idx, \(x) locate_snps_within_genes(snps, annotation, x))
idx_within_gene <- idx_within_gene[!is.na(idx_within_gene)]
length(idx_within_gene)

# Get chromosome, position, REF allele and ALt allele for each outlier found within gene regions
idx_within_gene_chrom <- str_remove(as.character(snps[,idx_within_gene]$chromosome), "\\.1$")
idx_within_gene_pos <- snps[,idx_within_gene]$position
idx_within_gene_ref <- snps$other$ref[idx_within_gene]
idx_within_gene_alt <- snps$other$alt[idx_within_gene]

# Function to filter annotation file by the outliers within genes
# Arguments
# ANNOTATION: data.frame containing annotation info
# CHROM: string denoting chromosome ID
# POS: integer denoting SNP position on chromosome
subset_annotation_file <- function(ANNOTATION, CHROM, POS) {
  df <- ANNOTATION |> 
    filter(Chromosome == CHROM) |> 
    filter(POS >= Start_position & POS <= End_position) |> 
    mutate(SNP_position = POS)

  return(df)
}

# Execute function on all outliers and add results to a new data.frame
annotation_outlier <- map2_df(idx_within_gene_chrom, idx_within_gene_pos, ~ subset_annotation_file(annotation, .x, .y))

# Add column for gene symbol
annotation_outlier$Symbol <- str_extract(annotation_outlier$Gene, "(?<=\\.).+?(?=\\.)")
head(annotation_outlier)

# ---------- #
# Enrichr
# ---------- #

# Load package
library(enrichR)

# Select databases
dbs <- c("GO_Biological_Process_2025","MSigDB_Hallmark_2020")

# Get a vector of all unique gene symbols from annotation_outlier data.frame
outlier_symbols <- annotation_outlier$Symbol |> unique() |> sort()
length(outlier_symbols)

# Use all genes annotated in genome as background list
background_genes <- str_extract(annotation$Gene, "(?<=\\.).+?(?=\\.)") |> unique() |> sort()
length(background_genes)

# Run enrichR (returns a list, one for each library in dbs)
enrichr_res <- enrichr(outlier_symbols, dbs, background = background_genes, include_overlap = TRUE)

# Extract significant pathways for each library
enrichr_sig <- enrichr_res |> 
  map_dfr(~ .x, .id = "Library") |> 
  dplyr::select(Library:Adjusted.P.value, Odds.Ratio:Genes) |> 
  filter(Adjusted.P.value < alpha) |> 
  arrange(desc(Combined.Score))

# Visualise
plotEnrich(
  df = enrichr_sig,
  showTerms = 20, numChar = 40,
  y = "Count", orderBy = "Combined.Score"
)

# ---------- #
# Manhattan Plot
# ---------- #

# Create data.frame containing SNP_ID, chromosome, SNP position and P-value
manhattan_df <- data.frame(
  SNP_ID = locNames(snps),
  chromosome = snps$chromosome,
  position = snps$position,
  Pvalue = rda_pvals$p.values,
  RDA_outlier = ifelse(1:length(rda_pvals$p.values) %in% rda_outliers_idx, "Yes", "No")
)

# Read in chromosome lookup file
chromosome_lookup <- read.delim("data/WGS/chromosome_lookup.tsv")[1:29, c("Chromosome.name","GenBank.seq.accession","Seq.length")]

# Add a column for the chromosome numbers
manhattan_df <- left_join(manhattan_df, chromosome_lookup, by = join_by(chromosome == GenBank.seq.accession))

# Convert X chromosome to numeric
manhattan_df <- manhattan_df |> 
  drop_na() |> 
  mutate(Chromosome.name = factor(Chromosome.name, c(1:28, "X")))

# Calculate cumulative position for each SNP
manhattan_df_cum <- manhattan_df |> 
  group_by(Chromosome.name) |> 
  summarise(chr_length = max(position), .groups = "drop") |> 
  mutate(total = cumsum(chr_length) - chr_length) |> 
  left_join(manhattan_df, ., by = "Chromosome.name") |> 
  mutate(basepair_cum = position + total)
head(manhattan_df_cum)

# Get chromosome centre positions for x-axis labels
xaxis_df <- manhattan_df_cum |> 
  group_by(Chromosome.name) |> 
  summarise(centre = mean(basepair_cum), .groups = "drop")
xaxis_df

# Genes in key biological pathways
# retinal_genes <- c("ROBO2","PTPRM","EPHB1")
# chemical_synaptic_trans_genes <- c("USP46","GRID2","NLGN1","MCTP1","CLSTN2","GRIK3","GRIK4","BTBD9","SYN3","SLC8A3","GRM7","NRG3","DLGAP3")
# cell_to_cell_adhesion_genes <- c("KIRREL3","PTPRT","ROBO2","GRID2","NLGN1","PCDHGB4","CRTAM","PTPRM","ROBO1","MPZ","IL1RAPL1","CDH23","CDH12","CNTN4","FAT4")
# ERBB4_genes <- c("NRG3","ERBB4","NRG2")
# synapse_vesicle <- c("NLGN1","SYN3","SYN2")
genes_enrichr <- c("ROBO2","EPHA3","EPHA6","NLGN1","BTBD9","SLC8A3","PAX3","COL4A3","CDH12","GRIK3","DLGAP3","SYN3","EPHB1","CNTN4","DLG5","NRG3","PTPRT","NRG2","PCDHGB4","CRTAM","GRIK4","ATP6V0A4","NTN3","IL1RAPL2","IL1RAPL1","ANOS1")
genes_other <- c("MDGA2","TTC6","SERPINB1","MYLK4","PUM1","MATN1","ZW10","USP2")
all_enrichr_genes <- unique(unlist(strsplit(enrichr_sig$Genes, ";")))
(key_genes_df <- annotation_outlier |> 
  mutate(SNP_ID = paste(paste0(Chromosome, ".1"), SNP_position, sep = "_")) |> 
  left_join(manhattan_df_cum, by = "SNP_ID") |> 
  filter(Symbol %in% c(genes_enrichr,genes_other)) |> 
  # filter(Symbol %in% outlier_symbols) |> 
  dplyr::select(basepair_cum, Symbol) |> 
  distinct(Symbol, .keep_all = TRUE)
)

# Filter key genes data.frame
# key_genes_df <- filter(key_genes_df, basepair_cum > 409034402 & basepair_cum < 459034402)

# Manhattan plot
(man_plot <- ggplot(data = manhattan_df_cum, aes(x=basepair_cum, y=-log10(Pvalue)))+
  geom_point(aes(colour = as.factor(Chromosome.name)), alpha = 0.7, size = 1.5, show.legend = FALSE)+
  scale_colour_manual(
    values = rep(c("#2166ac", "#053061"), length.out = length(unique(manhattan_df_cum$Chromosome.name))),
  )+
  geom_hline(aes(yintercept=-log10(thres_env), linetype = "Significance Threshold"), colour = "red", linewidth = 0.75)+
  geom_vline(
    data = key_genes_df,
    aes(xintercept=basepair_cum),
    linewidth = 0.1, colour = "grey50", linetype = "dashed"
  )+
  ggrepel::geom_label_repel(
    data = key_genes_df,
    aes(x = basepair_cum, y = 24.5, label = paste0("bolditalic('", Symbol, "')")),
    box.padding = 0.4,
    point.padding = 0.1,
    alpha = 0.95,
    max.overlaps = Inf,
    size = 3,
    min.segment.length = 0,
    parse = TRUE
  )+
  scale_linetype_manual(name = "", values = 1)+
  scale_x_continuous(
    label = xaxis_df$Chromosome.name,
    breaks = xaxis_df$centre,
    expand = c(0.01, 0.01)
  )+
  scale_y_continuous(
    limits = c(0, 28),
    expand = c(0, 0)
  )+
  labs(x = "Chromosome", y = expression(-log[10](italic(P))))+
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.border = element_blank(),
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 12),
  )
)

# ---------- #
# Enrichr Plot
# ---------- #

# Create enrichr plotting data.frame
enrichr_sig_df <- enrichr_sig |> 
  mutate(
    # Create column for % DEGs in set
    Hit_Genes = as.numeric(str_extract(Overlap, "^[0-9]+")),
    Total_Genes = as.numeric(str_extract(Overlap, "(?<=/)[0-9]+")),
    Percent_Hit = Hit_Genes / Total_Genes * 100,
    # Reorder data.frame by highest combined score
    Term = fct_reorder(Term, Combined.Score, .desc = FALSE)
  )
enrichr_sig_df

# Dot plot
(enrichr_plt <- ggplot(data = enrichr_sig_df, aes(x = Term, y = Library))+
  geom_point(
    aes(fill = Combined.Score, size = Percent_Hit),
    shape = 21, stroke = 0.2, colour = "black",
  )+
  scale_fill_distiller(
    name = "Combined Score",
    palette = "YlOrRd",
    direction = 1
  )+
  scale_size_continuous(
    name = "Genes in Set",
    breaks = c(10, 20, 50),
    labels = c("10%", "20 %", "50 %"),
    range = c(3, 10)
  )+
  coord_flip()+
  guides(size = guide_legend(reverse = TRUE))+
  scale_y_discrete(labels = c("GO Biological Process"))+
  xlab("Biological Pathway (Term)")+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.background = element_rect(fill = "#F5F5F5"),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, vjust = -0.5),
    axis.text.y = element_text(size = 13, vjust = 0.4, colour = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
  )
)

# ---------- #
# Figure 5
# ---------- #

# Prepare figure using patchwork
library(patchwork)

# Figure composition
layout <- "
A
B
"
figure5 <- wrap_plots(
  A = man_plot+ labs(tag = "A"),
  # free() removes alignment from plot
  B = free(enrichr_plt+ labs(tag = "B")),
  design = layout
)
ggsave("figures/Figure5.png", plot = figure5, width = 11, height = 12, units = "in", dpi = 600)
ggsave("figures/Figure5.pdf", plot = figure5, width = 11, height = 12, units = "in")