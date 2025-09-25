# ========== #
#
# Rhipposideros RNA Analysis 2025
#
# Gene Expression Modelling
#
# ========== #

# Import libraries
import os
import re
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import seaborn as sns

# Import gene expression packages
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import make_MA_plot
import gseapy as gp
# from gseapy import barplot, dotplot

# Set working directory to file source location
wd = os.getcwd()
os.chdir(wd)

# Save results
SAVE = True

# Set output directory
if SAVE:
    OUT_DIR = "outputs/"
    os.makedirs(OUT_DIR, exist_ok=True)

# Path to input data
DATA_DIR = "data/RNA/"

# ---------- #
# Data Preprocessing
# ---------- #

# Read in CSV count data
counts_df = pd.read_csv(os.path.join(DATA_DIR, "kallisto_quant_results_Rhipposideros.csv"), index_col=0)
# counts_df

# Extract gene ID from index (e.g. KCNMA1) and add column to DataFrame
geneIDs = []
for string in counts_df.index.tolist():
    gene = re.findall(r"\.(.*?)\.", string)
    geneIDs.append(gene[0])
counts_df["Gene_ID"] = geneIDs
# counts_df

# Group by gene and sum the counts for each column
counts_df = counts_df.groupby("Gene_ID").sum()
counts_df.index.name = ""

# Convert floats to integers
counts_df = counts_df.astype(int)

# Remove library A and rename library B for samples CH172 and CH188
counts_df = counts_df.drop(columns=["CH172A", "CH188A"])
counts_df = counts_df.rename(columns={"CH172B": "CH172", "CH188B": "CH188"})

# Transpose DataFrame so genes are columns and samples are rows
counts_df = counts_df.T.sort_index(ascending=True)

# Read in Excel metadata
metadata = pd.read_excel("../Bat_sampling_database.xlsx", sheet_name="Blood_Sampling")
# metadata

# Subset RNA-seq samples for analysis
metadata = (
    metadata[metadata["RNAseq"] == "Yes"] # chaining using Pandas
    [["Sample_ID", "Treatment"]] # select only these columns
    .sort_values("Sample_ID") # order alphabetically
    .set_index("Sample_ID") # set sample ID column as index
)

# ---------- #
# Filtering
# ---------- #

# Omit genes with less than 10 counts in total
count_threshold = 10
genes_to_keep = counts_df.sum(axis=0) > count_threshold
counts_df_filt = counts_df[counts_df.columns[genes_to_keep]]
counts_df_filt

# Print original and filtered number of genes
print("Pre-Filter: ", len(counts_df.columns))
print("Post-Filter:", len(counts_df_filt.columns))

# ---------- #
# Single Factor Model
# ---------- #

# Model design
design = "~Treatment"

# Create a DeseqDataSet object
dds = DeseqDataSet(
    counts=counts_df_filt,
    metadata=metadata,
    design=design,
    refit_cooks=True,
    n_cpus=4
)

# Run DESeq2 model
dds.deseq2()

# Fit size factors
dds.fit_size_factors()

# Calculate VST counts
dds.vst()

# Save results using pickle.dump
if SAVE:
    with open(os.path.join(OUT_DIR, "Rhipposideros_dds.pkl"), "wb") as file:
        pkl.dump(dds, file)

# # Load previous results using pickle.load
# with open(os.path.join(OUT_DIR, "Rhipposideros_dds.pkl"), "rb") as file:
#     dds = pkl.load(file)

# Print general information about model
dds # model summary
dds.X # count data in numpy array object
dds.obs # metadata factors
dds.var_names # gene names
dds.obsm["design_matrix"] # design factors
dds.varm["dispersions"] # dispersion values
dds.varm["LFC"] # original log fold change values
dds.layers["normed_counts"] # normalised counts
dds.layers["vst_counts"] # VST counts

# Print size factor estimates
np.sort(dds.obsm["size_factors"]) # size factor estimates
size_factors = {
    "sample": dds.obs.index.tolist(),
    "size_factor": dds.obsm["size_factors"],
} # create dictionary then convert to pandas DataFrame
pd.DataFrame(size_factors).sort_values(by="size_factor")


# ---------- #
# Statistical Tests
# ---------- #

# Thresholds
alpha = 0.05
lfc = 1

# Factor to contrast
contrast = ["Treatment", "Light", "Dark"]

# Run Wald tests
results = DeseqStats(dds, contrast=contrast, alpha=alpha)
results.summary()
results_df = results.results_df

# Run LFC shrinkage
results.lfc_shrink(coeff="Treatment[T.Light]")

# Diagnosis plots
dds.plot_dispersions()
make_MA_plot(results_df)

# Highlight upregulated or downregulated genes
upregulated = results_df[
    (results_df["log2FoldChange"] >= lfc) & (results_df["padj"] < alpha)    
]
downregulated = results_df[
    (results_df["log2FoldChange"] <= -lfc) & (results_df["padj"] < alpha)    
]

# Volcano plot
plt.scatter(x=results_df["log2FoldChange"], y=-np.log10(results_df["padj"]), s=5, color="grey", label="Not significant")
plt.scatter(x=upregulated["log2FoldChange"], y=-np.log10(upregulated["padj"]), s=5, color="red", label="Upregulated")
plt.scatter(x=downregulated["log2FoldChange"], y=-np.log10(downregulated["padj"]), s=5, color="blue", label="Downregulated")
plt.xlabel("logFC")
plt.ylabel("-log10 Padjust")
plt.axvline(-lfc, color="grey", linestyle="--")
plt.axvline(lfc, color="grey", linestyle="--")
plt.axhline(-np.log10(alpha), color="grey", linestyle="--")
plt.legend(loc="best")

# ---------- #
# Differentially Expressed Genes (DEGs)
# ---------- #

# Check LFC shrinkage was run
results.shrunk_LFCs

# Subset DEGs from DataFrame
DEGs = results_df.loc[
    lambda df: (df["padj"] < alpha) & (df["log2FoldChange"].abs() > lfc)
].sort_values(by=["padj", "log2FoldChange"], ascending=[True,False])
# DEGs

# Save DEG results
if SAVE:
    DEGs.to_csv(os.path.join(OUT_DIR, "Rhipposideros_DEGs.csv"), index=True)

# ---------- #
# Gene Enrichment Analysis
# ---------- #

# List of DEGs as Entrez gene symbols
DEGs_symbols = DEGs.index.tolist()

# List of background genes (expressed in blood)
background_genes = dds.var_names.tolist()

# # Parse EnrichR human libraries into dictionary (run once then comment out)
# # https://maayanlab.cloud/Enrichr/
# GO_Biological_Process = gp.get_library(name="GO_Biological_Process_2023", organism="Human")
# MSigDB_Hallmark = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")

# # Merge dictionaries into a single dictionary
# enrichr_database = {**GO_Biological_Process, **MSigDB_Hallmark}

# # Save database using pickle.dump
# if SAVE:
#     with open(os.path.join(OUT_DIR, "enrichr_database.pkl"), "wb") as file:
#         pkl.dump(enrichr_database, file)

# Load database using pickle.load
with open(os.path.join(OUT_DIR, "enrichr_database.pkl"), "rb") as file:
    enrichr_database = pkl.load(file)

# # Run enrichment analysis offline
# enrich = gp.enrich(
#     gene_list=DEGs_symbols,
#     background=background_genes,
#     gene_sets=enrichr_database,
#     outdir=None,
#     verbose=True,
# )

# Run enrichment analysis through web server
gene_sets = ["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
enrich = gp.enrichr(
    gene_list=DEGs_symbols,
    background=background_genes,
    gene_sets=gene_sets,
    organism="Human",
    outdir=None,
    verbose=True,
)

# Filter results by alpha threshold and sort by combined score
enrich_filt = enrich.results.loc[
    lambda df: df["Adjusted P-value"] < alpha
].sort_values(by="Combined Score", ascending=False)

# Add percentage of genes in gene set to DataFrame
terms = enrich_filt["Term"].tolist()
percent_terms = []
for term in terms:
    all_genes_in_set = enrichr_database[term]
    DEGs_in_set = enrich_filt[enrich_filt["Term"] == term]["Genes"].tolist()[0].split(";")
    common_genes = len(set(all_genes_in_set) & set(DEGs_in_set)) # count how many genes are in both lists
    percentage = round((common_genes / len(all_genes_in_set) * 100), 1) # as a percentage of set
    percent_terms.append(percentage)
enrich_filt["Percentage DEGs in Gene Set"] = percent_terms

# Save gene enrichment results
if SAVE:
    enrich_filt.to_csv(os.path.join(OUT_DIR, "Rhipposideros_enrichr.csv"), index=False)

