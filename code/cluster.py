
# IMPORT LIBRARIES
import pandas as pd

from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np
import hdbscan
import scanpy as sc
import seaborn as sns

# DOWNLOAD EXPRESSION MATRIX
expr = pd.read_csv('../data/expression_annotated_corrected.csv')
print("Shape of expression matrix:", expr.shape)

# SCALING
scaler = StandardScaler()
expr_scale = scaler.fit_transform(expr)

# CLUSTERING
tags_name = np.array(expr.columns)
expr_scale = pd.DataFrame(expr_scale, columns = tags_name)
adata = sc.AnnData(X = expr_scale)

sc.pp.neighbors(adata, n_pcs = 15)
sc.tl.leiden(adata, key_added = "Clusters")

sc.tl.umap(adata)  # Compute UMAP
sc.pl.umap(adata, color = "Clusters", legend_loc = "on data")
plt.savefig("umap_unlabeled", dpi = 300, bbox_inches = "tight")
plt.close()

def umap_stacked(adata):
    # Define the number of plots per row
    n_cols = 4  # Adjust based on preference
    n_plots = len(tags_name)
    n_rows = (n_plots // n_cols) + (n_plots % n_cols > 0)  # Calculate number of rows needed

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))  # Adjust figure size

    # Flatten axes for easy iteration
    axes = axes.flatten()

    for i, tag in enumerate(tags_name):
        sc.pl.umap(
            adata,
            color = tag,
            vmin = 0,
            vmax = "p99",
            sort_order = False,
            frameon = False,
            cmap = "Reds",
            show = False,  # Prevent immediate display
            ax = axes[i]   # Assign to correct subplot
        )

    # Hide unused subplots (if any)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig("../images/umap_stacked.png", dpi=300, bbox_inches="tight")  # Save the figure
    plt.close()

def genes_dot(adata, n_genes = 5):
    sc.tl.rank_genes_groups(adata, groupby = "groups", method = "wilcoxon")
    sc.pl.rank_genes_groups_dotplot(adata, n_genes)
    plt.savefig("../images/ranked_genes_dotplot", dpi = 300, bbox_inches = "tight")
    plt.close()

markers_to_cells = {
    "17": "Cytotoxic T cells (CD8a)",
    "13": "B cells (CD20)",
    "16": "Macrophages (CD68)",
    "19": "Epithelial cells (E/P Cadherin)",
    "6": "Breast cancer cells (PanCK, GATA3)",
    "10": "Potentially breast cancer (PanCK, GATA3)",
    "15": "Potentially breast cancer (PanCK, GATA3, CAIX)",
    "0": "Some T cells (CD3)",
    "12": "MMP9 High cells",
    "11": "Potentional fibroblasts (fibronectin)",
    "2": "Potentially fibroblasts (fibronectin)",
    "9": "Potentially fibroblasts (fibronectin)",
    "4": "Potentially fibroblasts (fibronectin)",
    "14": "Potentially malignant cells (CD34, Podoplanin)"
}

adata.obs["manual_celltype_annotation"] = adata.obs.groups.map(markers_to_cells)
sc.pl.umap(adata, color=["manual_celltype_annotation"])