# ========== #
#
# Rhipposideros RNA Analysis 2025
#
# Gene Expression Visualisation
#
# ========== #

# Import libraries
import os
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec
import seaborn as sns
from sklearn.decomposition import PCA

# Set working directory to file source location
wd = os.getcwd()
os.chdir(wd)

# Path to input data
INPUT_DIR = "outputs/"

# ---------- #
# DEGs: Scatter Plot of Counts
# ---------- #

# Read in DEGs
DEGs = pd.read_csv(os.path.join(INPUT_DIR, "Rhipposideros_DEGs.csv"), index_col=0)

# The only downregulated gene
downregulated_gene = "PTGS2"

# Read in dds object
with open(os.path.join(INPUT_DIR, "Rhipposideros_dds.pkl"), "rb") as file:
    dds = pkl.load(file)

# Function to get counts for a gene
def get_gene_counts(deseq, gene=None, counts="normalised"):
    """
    Returns a numpy array of normalised counts
    """
    # Get column index for gene
    gene_idx = deseq.var_names.get_loc(gene)
    
    # Count data to return
    counts_type = "normed_counts" if counts == "normalised" else "vst_counts"
    
    # Extract normalised counts from numpy array
    norm_counts = dds.layers[counts_type][:, gene_idx]
    
    # Return gene counts
    return(norm_counts.round(1))
        
    
# Test function
# get_gene_counts(dds, "HBB")

# DataFrame for plotting
plt_counts = pd.DataFrame(
    {
     "Sample": dds.obs.index.tolist(),
    "Treatment": dds.obs["Treatment"].tolist(),
    }
)

# Gene to plot
gene = "DDB1"
gene = "BSG"
gene = "SOD1"
gene = "PDK2"
gene = "SIRT2"
gene = "PARK7"
gene = "CAT"
gene = "GADD45A"
gene = "PRDX2"
gene = "GPX1"

# Add gene counts to DataFrame
plt_counts[gene] = get_gene_counts(dds, gene, counts="normalised")
# plt_counts

# Palette
palette = {"Light": "#FFA500", "Dark": "black"}

# Plot jitter per gene
sns.set_theme(style="darkgrid")
plt.figure(figsize=(6, 5))
sns.stripplot(x="Treatment", y=gene, data=plt_counts, s=8, hue="Treatment", palette=palette, edgecolor="white", linewidth=0.4) 
plt.title(f"Gene expression: {gene}")
plt.xlabel("")
plt.ylabel("Normalised counts\n")
plt.tight_layout()
# plt.savefig(f"{gene}.jpeg", dpi=300)

# List of genes to plot altogether
# genes_to_plot = DEGs.index.tolist()[0:30]
genes_to_plot = ["DDB1","BSG","SOD1","PDK2","SIRT2","PARK7","CAT","GADD45A","PRDX2","GPX1"]

# Extract normalised counts for these genes
count_array = [get_gene_counts(dds, gene, "vst_counts") for gene in genes_to_plot]

# Create a DataFrame of counts and rename columns and rows
count_df = pd.DataFrame(count_array).T # transpose DataFrame
count_df.columns = genes_to_plot
count_df.index = dds.obs.index.tolist()
count_df["Treatment"] = dds.obs["Treatment"].tolist() # add column for treatment

# Convert DataFrame to long format using melt()
count_df_long = count_df.reset_index(names="Sample").melt(
    id_vars=["Sample","Treatment"], var_name="Gene", value_name="Counts"
)

# Plot jitter for all genes in list
sns.set_theme(style="darkgrid")
plt.figure(figsize=(12, 6))
sns.stripplot(data=count_df_long, x="Gene", y="Counts", s=6, alpha=0.9, hue="Treatment", palette=palette)
plt.xlabel("")
plt.xticks(rotation=90)
plt.ylabel("VST Counts")
plt.title("Gene Expression Profiles Between Light and Dark Group")
plt.tight_layout()
# plt.savefig("DEGs_counts.jpeg", dpi=300)

# ---------- #
# DEGs: Sample Positions
# ---------- #

# Get normalised counts for all genes in DataFrame
all_gene_counts = [get_gene_counts(dds, gene, "normalised") for gene in DEGs.index.tolist()]
all_gene_counts = pd.DataFrame(all_gene_counts).T # transpose DataFrame
all_gene_counts.columns = DEGs.index.tolist()
all_gene_counts.index = dds.obs.index.tolist()

# Find index of maximum value and save the row index in list
sample_max_count = []
for col in all_gene_counts.columns:
    max_idx = all_gene_counts[col].idxmax()
    sample_max_count.append(max_idx)
    # print(f"Column '{col}': max value at index '{max_idx}'")

# Length of new list
len(sample_max_count)

# Print table where each count represents the sample with the highest expression 
pd.DataFrame(sample_max_count, columns=["Sample"]).value_counts()

# Subset individuals exposed to artificial light
light = all_gene_counts.loc[dds.obs[dds.obs["Treatment"] == "Light"].index]

# Remove PTGS2 downregulated gene from the light expression matrix
light = light.drop(columns=["PTGS2"])

# Calculate the mean of expression in Light samples for each gene
light_mean = light.mean()

# For each Light individual, mark as upregulated if expression > mean of Light group
upregulated = light.gt(light_mean)

# Count in how many light samples each gene is upregulated
n_samples_upregulated = upregulated.sum(axis=0)

# Summarise how many genes are upregulated in 10 Light individuals
n_samples_upregulated.value_counts().sort_index()
upregulated.sum(axis=1)
    
# ---------- #
# DEGs: Scatter Plot of PCA
# ---------- #

## Extract VST counts for all samples and DEGs

# Create a dictionary mapping gene names to their index in dds.var_names
gene_to_index = {gene: i for i, gene in enumerate(dds.var_names)}

# Get indices in the original order of DEGs.index
DEG_idx = [gene_to_index[gene] for gene in DEGs.index if gene in gene_to_index]

# Get the normalised counts for DEGs (subset the numpy array)
DEG_vst_counts = dds.layers["vst_counts"][:, DEG_idx]
DEG_norm_counts = dds.layers["normed_counts"][:, DEG_idx]

# Principal components analysis
pca = PCA(n_components=3)
DEG_pca = pca.fit_transform(DEG_norm_counts)

# Get percentage of variance explained
percent_variance = (pca.explained_variance_ratio_ * 100).round(1)
percent_variance

# Add samples to DataFrame
DEG_pca_pd = pd.DataFrame(DEG_pca)
DEG_pca_pd.index = dds.obs.index
DEG_pca_pd.columns = ["PC1", "PC2", "PC3"]
DEG_pca_pd["Treatment"] = dds.obs["Treatment"]

# # Map colours to Treatment list
# pca_colour_map = dds.obs["Treatment"].map(palette).tolist()

# # Plot the first two PCs
# plt.scatter(
#     x=DEG_pca_pd[0], y=DEG_pca_pd[1],
#     c=pca_colour_map,
#     marker="o",
#     s=50,
#     edgecolor="white",
#     linewidth=0.4,
# )
# plt.xlabel(f"PC1 ({percent_variance[0]}%)")
# plt.ylabel(f"PC2 ({percent_variance[1]}%)")
# plt.title("PCA of DEGs using normalised counts")
# plt.tight_layout()
# plt.savefig("outputs/RNA_DEGs_PCA.jpeg", dpi=300)

# Plot PCA using Seaborn
sns.set_theme(style="darkgrid")
plt.figure(figsize=(8, 5))
sns.scatterplot(
    data=DEG_pca_pd,
    x="PC1",
    y="PC2",
    hue="Treatment",          # Color by treatment group
    palette=palette,          # Use your predefined palette
    s=50,                     # Point size
    edgecolor="white",        # White outline
    linewidth=0.4
)

# Labels and title
plt.xlabel(f"PC1 ({percent_variance[0]}%)")
plt.ylabel(f"PC2 ({percent_variance[1]}%)")
plt.title("PCA of DEGs using normalised counts")
plt.legend(title="Treatment", bbox_to_anchor=(0.85, 1), loc="upper left")
plt.tight_layout()

# Save plot
plt.savefig("outputs/RNA_DEGs_PCA.jpeg", dpi=300)

# ---------- #
# Gene Enrichment
# ---------- #

# Read in enrichment results
enrichr = pd.read_csv(os.path.join(INPUT_DIR, "Rhipposideros_enrichr.csv"))

# Convert DataFrame to dictionary
enrichr = pd.DataFrame.to_dict(enrichr, orient="list")

# Gene set order on x axis
x_order = ["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
x_order_labs = [s.replace("_", " ").replace(" 2023", "").replace(" 2020", "") for s in x_order]

# Convert gene set column to indices based on the order
gene_set_to_idx = {gene: idx for idx, gene in enumerate(x_order)}
x_vals = [gene_set_to_idx[gene] for gene in enrichr["Gene_set"]]

# Extract y values, color aes, size aes
y_vals = enrichr["Term"]
color_by = enrichr["Combined Score"]
size_by = enrichr["Percentage DEGs in Gene Set"]
size_by_scale = 4

# Colour palette
cmap = cm.YlOrRd

# Create figure and axes
fig, ax = plt.subplots(figsize=(10, 6))

# Create a scatter plot using Matplotlib's scatter()
scatter = ax.scatter(
    x_vals,
    y_vals,
    c=color_by, # Hue is the combined score
    s=np.array(size_by)*size_by_scale, # Size is percentage DEGs in Gene Set
    cmap=cmap,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.9,
)

# Set x-ticks with gene sets on x-axis
ax.set_xticks(np.arange(len(x_order)))
ax.set_xticklabels(x_order_labs)

# Adjust margin between gene sets
x_margin = 0.55
ax.set_xlim(-x_margin, 1 + x_margin)

# Colorbar legend
cbar = plt.colorbar(scatter, ax=ax, fraction=0.05, aspect=17, anchor=(0.0, 0.1))
cbar.set_label("\nCombined Score", rotation=270, labelpad=35)

# Size legend values and scaling
size_legend_values = [20, 10, 5]
size_legend_scaled = [s * size_by_scale for s in size_legend_values]

# Create legend handles
legend_handles = [
    plt.scatter([], [], s=s, edgecolor="black", color="white", label=f"{v} %")
    for s, v in zip(size_legend_scaled, size_legend_values)
]

# Add legend handles to plot
ax.legend(
    handles=legend_handles,
    title="Genes in Set", title_fontsize=12,
    loc="upper right", bbox_to_anchor=(1.40, 0.95),
    fontsize=12, handletextpad=0.1,
)

# Adjust the y-axis labels
ax.set_yticks(np.unique(y_vals)) # Make sure all Terms are displayed
ax.set_yticklabels(np.unique(y_vals))
ax.tick_params(axis="y", labelsize=10, labelcolor="black")

# Customise plot
ax.set_xlabel("")
ax.set_ylabel("Biological Pathway (Term)\n")
ax.set_title("Biological Pathways Statistically Impacted by Artificial Light")
ax.invert_yaxis()
plt.tight_layout()
plt.subplots_adjust(right=0.90) # Adjust right margin (default is 1.0)
# plt.savefig("DEGs_enrichr.jpeg", dpi=300)

# ---------- #
# Figure 2
# ---------- #

# DataFrame for plotting
plt_counts = pd.DataFrame(
    {
     "Sample": dds.obs.index.tolist(),
    "Treatment": dds.obs["Treatment"].tolist(),
    }
)

# Add gene counts to DataFrame
plt_counts["DDB1"] = get_gene_counts(dds, "DDB1", counts="normalised")
plt_counts["GADD45A"] = get_gene_counts(dds, "GADD45A", counts="normalised")
plt_counts["BSG"] = get_gene_counts(dds, "BSG", counts="normalised")
plt_counts["H2BC11"] = get_gene_counts(dds, "H2BC11", counts="normalised")
plt_counts["MACROD1"] = get_gene_counts(dds, "MACROD1", counts="normalised")
plt_counts["MACROH2A1"] = get_gene_counts(dds, "MACROH2A1", counts="normalised")
plt_counts["UBE2B"] = get_gene_counts(dds, "UBE2B", counts="normalised")
plt_counts["C1D"] = get_gene_counts(dds, "C1D", counts="normalised")
plt_counts["UBE2B"] = get_gene_counts(dds, "UBE2B", counts="normalised")
plt_counts["IFI27"] = get_gene_counts(dds, "IFI27", counts="normalised")
plt_counts["CAT"] = get_gene_counts(dds, "CAT", counts="normalised")
plt_counts["SOD1"] = get_gene_counts(dds, "SOD1", counts="normalised")
plt_counts["PARK7"] = get_gene_counts(dds, "PARK7", counts="normalised")
plt_counts["PDK2"] = get_gene_counts(dds, "PDK2", counts="normalised")
plt_counts["PRDX2"] = get_gene_counts(dds, "PRDX2", counts="normalised")
plt_counts["SIRT2"] = get_gene_counts(dds, "SIRT2", counts="normalised")
plt_counts

# Figure with custom grid specification
fig = plt.figure(figsize=(8,14))
sns.set_theme(style="darkgrid")

# Create a gridspec for the top 4 plots (two rows, two columns)
gs_top = gridspec.GridSpec(
    nrows=2, ncols=3,
    figure=fig,
    # height_ratios=[1, 1],
    wspace=0.1,
    left=0, right=1, bottom=0.55, top=0.95
)

# Create a separate gridspec for the bottom plot (one row, one column)
gs_bottom = gridspec.GridSpec(
    nrows=1, ncols=1,
    figure=fig,
    # height_ratios=[2],
    left=0.60, right=0.90, bottom=0.10, top=0.50
)

# Row1 Col1
ax00 = fig.add_subplot(gs_top[0, 0])
sns.stripplot(
    ax=ax00,
    x="Treatment", y="DDB1", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax00.set_xlabel("")
ax00.set_xticklabels("")
ax00.set_yticklabels("")
ax00.set_ylabel("Normalised counts")
ax00.set_title("DDB1", fontdict={"style": "italic"})

# Row1 Col2
ax01 = fig.add_subplot(gs_top[0, 1])
sns.stripplot(
    ax=ax01,
    x="Treatment", y="BSG", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax01.set_xlabel("")
ax01.set_xticklabels("")
ax01.set_yticklabels("")
ax01.set_ylabel("")
ax01.set_title("BSG", fontdict={"style": "italic"})

# Row1 Col3
ax02 = fig.add_subplot(gs_top[0, 2])
sns.stripplot(
    ax=ax02,
    x="Treatment", y="IFI27", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax02.set_xlabel("")
ax02.set_xticklabels("")
ax02.set_yticklabels("")
ax02.set_ylabel("")
ax02.set_title("IFI27", fontdict={"style": "italic"})

# Row2 Col1
ax10 = fig.add_subplot(gs_top[1, 0])
sns.stripplot(
    ax=ax10,
    x="Treatment", y="SOD1", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax10.set_xlabel("")
# ax10.set_xticklabels("")
ax10.set_yticklabels("")
ax10.set_ylabel("Normalised counts")
ax10.set_title("SOD1", fontdict={"style": "italic"})

# Row2 Col2
ax11 = fig.add_subplot(gs_top[1, 1])
sns.stripplot(
    ax=ax11,
    x="Treatment", y="PARK7", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax11.set_xlabel("")
# ax11.set_xticklabels("")
ax11.set_yticklabels("")
ax11.set_ylabel("")
ax11.set_title("PARK7", fontdict={"style": "italic"})

# Row2 Col3
ax12 = fig.add_subplot(gs_top[1, 2])
sns.stripplot(
    ax=ax12,
    x="Treatment", y="H2BC11", data=plt_counts,
    s=8, hue="Treatment",
    palette=palette,
    edgecolor="white",
    linewidth=0.4
) 
ax12.set_xlabel("")
# ax12.set_xticklabels("")
ax12.set_yticklabels("")
ax12.set_ylabel("")
ax12.set_title("H2BC11", fontdict={"style": "italic"})

# Row2 Col1-2
ax3 = fig.add_subplot(gs_bottom[0, 0])
scatter = ax3.scatter(
    x=x_vals,
    y=y_vals,
    c=color_by, # Hue is the combined score
    s=np.array(size_by)*size_by_scale, # Size is percentage DEGs in Gene Set
    cmap=cmap,
    edgecolor="black",
    linewidth=0.5,
    alpha=0.9,
)

# Gene set order on x axis
x_order = ["GO_Biological_Process_2023", "MSigDB_Hallmark_2020"]
x_order_labs = [s.replace("GO_", "GO ").replace("_Process_2023", "\nProcess").replace("_Hallmark_2020", "\nHallmark") for s in x_order]
x_order_labs

# Set x-ticks with gene sets on x-axis
ax3.set_xticks(np.arange(len(x_order)))
ax3.set_xticklabels(x_order_labs)

# Adjust margin between gene sets
ax3.set_xlim(-0.75, 1.75)

# Colorbar legend
cbar = plt.colorbar(scatter, ax=ax3, fraction=0.05, aspect=17, anchor=(0.0, 0.5))
cbar.set_label("\nCombined Score", rotation=270, labelpad=35)

# Size legend values and scaling
size_legend_values = [20, 10, 5]
size_legend_scaled = [s * size_by_scale for s in size_legend_values]

# Create legend handles
legend_handles = [
    plt.scatter([], [], s=s, edgecolor="black", color="white", label=f"{v} %")
    for s, v in zip(size_legend_scaled, size_legend_values)
]

# Add legend handles to plot
ax3.legend(
    handles=legend_handles,
    title="Genes in Set", title_fontsize=12,
    loc="upper right", bbox_to_anchor=(1.60, 1.01),
    fontsize=12, handletextpad=0.1,
)

# Adjust the y-axis labels
ax3.set_yticks(np.unique(y_vals)) # Make sure all Terms are displayed
ax3.set_yticklabels(np.unique(y_vals))
ax3.tick_params(axis="y", labelsize=10, labelcolor="black")

# Customise plot
ax3.set_xlabel("")
ax3.set_ylabel("Biological Pathway (Term)\n")
# ax3.set_title("Biological Pathways Statistically Impacted by Artificial Light")
ax3.invert_yaxis()

# Add figure tags to the top-left of each subplot
# ax00.text(x=-0.3, y=1.10, s="A", color="black", fontsize=15, fontweight="bold", ha="left", va="top", transform=ax00.transAxes)
ax3.text(x=-2.4, y=2.20, s="A", color="black", fontsize=15, fontweight="bold", ha="left", va="top", transform=ax3.transAxes)
ax3.text(x=-2.4, y=1.05, s="B", color="black", fontsize=15, fontweight="bold", ha="left", va="top", transform=ax3.transAxes)

# Save figure
plt.savefig("figures/Figure2.png", dpi=600, bbox_inches="tight")
plt.savefig("figures/Figure2.pdf", bbox_inches="tight")
