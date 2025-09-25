# ========== #
#
# Rhipposideros WGS Analysis 2025
#
# SNP QC and Filtering
#
# ========== #

# In RStudio set working directory to the path where this R script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load packages
library(vcfR)
library(SNPfiltR)
library(adegenet)
library(pegas)
library(stringr)
library(dplyr)
library(tidyr)
library(mapmixture)
library(ggrepel)
library(dartR)
# dartR::gl.install.vanilla.dartR() # install dartR dependencies
library(snpStats) # BiocManager::install("snpStats")
library(qvalue) # BiocManager::install("qvalue")
library(OutFLANK) # devtools::install_github("whitlock/OutFLANK")

# ---------- #
# SNP QC
# ---------- #

# Read in variants
vcf <- read.vcfR("data/WGS/Rhipposideros_0.10missing_4mac.vcf.gz")
vcf

# Check for biallelic SNPs
summary(is.biallelic(vcf))

# Only keep biallelic SNPs
# vcf <- vcf[is.biallelic(vcf) == TRUE]

# Filter SNPs for linkage
vcf <- distance_thin(vcf, min.distance = 1000)

# Extract REF and ALT alleles
ref_alleles <- vcf@fix[, "REF"]
alt_alleles <- vcf@fix[, "ALT"]
snp_ids <- paste(vcf@fix[, "CHROM"], vcf@fix[, "POS"], sep = "_")

# Convert vcf to genlight object
snps <- vcfR2genlight(vcf)
snps

# Add REF and ALT allele information to genlight
match_ids <- match(locNames(snps), snp_ids)
snps@other$ref <- ref_alleles[match_ids]
snps@other$alt <- alt_alleles[match_ids]

# Missing data for individuals
indmiss <- rowSums(is.na(as.matrix(snps))) / adegenet::nLoc(snps)

# Plot using barplot
barplot(indmiss, ylim = c(0,1), ylab = "Complete genotypes (proportion)", las = 2)
indmiss[ which.max(indmiss) ] # maximum missing genotypes

# Add site groups
snps$pop <- as.factor(gsub("[^a-zA-Z]", "", indNames(snps)))
summary(snps$pop)

# Check monomorphic SNPs
gl.report.monomorphs(snps)

# Filter loci that depart from HWE
snps_hwe <- gl.filter.hwe(snps, n.pop.threshold = 10, min_sample_size = 5)

# ---------- #
# Outlier Loci
# ---------- #

# Convert genlight to geno format (one column per SNP)
gl2geno(snps_hwe, outfile = "outflank_snps_hwe", outpath = "data/WGS/")

# OutFLANK
SNPmat <- read.table("data/WGS/outflank_snps_hwe.lfmm")
fst_mat <- MakeDiploidFSTMat(SNPmat, locNames(snps_hwe), snps_hwe$pop)
outflank_results <- OutFLANK(fst_mat, NumberOfSamples = nPop(snps_hwe))
OutFLANKResultsPlotter(outflank_results)
summary(outflank_results$results$OutlierFlag)

# Outliers
snps_hwe_neutral <- snps_hwe[ ,which(!outflank_results$results$OutlierFlag) ]
snps_hwe_neutral

# ---------- #
# Export Genotypes
# ---------- #

# Export all SNPs as .RData file
snps <- snps[order(snps$ind.names)]
save(snps, file = "data/WGS/snps_genlight.RData")

# Export neutral SNPs as .geno and .lfmm file
gl2geno(snps, outfile = "snps", outpath = "data/WGS/")

# Export neutral SNPs as .RData file
snps_hwe_neutral <- snps_hwe_neutral[order(snps_hwe_neutral$ind.names)]
save(snps_hwe_neutral, file = "data/WGS/snps_hwe_neutral_genlight.RData")

# Export neutral SNPs as .geno and .lfmm file
gl2geno(snps_hwe_neutral, outfile = "snps_hwe_neutral", outpath = "data/WGS/")

# Export genlight object as PopCluster input file
# snps_hwe_neutral_gen <- gl2gi(snps_hwe_neutral)
# snps_hwe_neutral_df <- as.data.frame(snps_hwe_neutral_gen)
# snps_hwe_neutral_df$IND <- indNames(snps_hwe_neutral)
# snps_hwe_neutral_df$POP <- factor(snps_hwe_neutral$pop, labels = 1:nPop(snps_hwe_neutral))
# snps_hwe_neutral_df <- select(snps_hwe_neutral_df, IND, POP, everything())
# write.table(snps_hwe_neutral_df, file = "data/WGS/popcluster.txt", col.names = FALSE, row.names = FALSE)
