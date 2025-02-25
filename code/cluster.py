
# IMPORT LIBRARIES
import numpy as np
import pandas as pd
import scanpy as sc
import plotly.express as px
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from sklearn.preprocessing import StandardScaler

# DOWNLOAD EXPRESSION MATRIX, CELL COORDINATES
expr = pd.read_csv('../data/expression_annotated_corrected.csv')
coords = pd.read_csv('../data/cell_coordinates.csv', index_col = 0)
print("Shape of expression matrix:", expr.shape)

# VARIABLES
marker_to_cell = {
    "17": "Cytotoxic T cells (CD8a)",
    "13": "B cells (CD20)",
    "16": "Macrophages (CD68)",
    "19": "Epithelial cells (E/P Cadherin)",
    "6": "Breast cancer cells (PanCK, GATA3)",
    "10": "Potentially breast cancer (PanCK, GATA3)",
    "15": "Potentially breast cancer (PanCK, GATA3, CAIX)",
    "0": "Some T cells (CD3)",
    "12": "MMP9 High cells",
    "11": "Potentially fibroblasts (fibronectin)",
    "2": "Potentially fibroblasts (fibronectin)",
    "9": "Potentially fibroblasts (fibronectin)",
    "4": "Potentially fibroblasts (fibronectin)",
    "14": "Potentially malignant cells (CD34, Podoplanin)"
}

# SCALING
scaler = StandardScaler()
expr_scale = scaler.fit_transform(expr)

# CLUSTERING 
tags_name = np.array(expr.columns)
expr_scale = pd.DataFrame(expr_scale, columns = tags_name)
adata = sc.AnnData(X = expr_scale)

sc.pp.neighbors(adata, n_pcs = 15)
sc.tl.leiden(adata, key_added = "Clusters")
print(f"Number of clusters: {adata.obs['Clusters'].nunique()}")

sc.tl.umap(adata)  # Compute UMAP
sc.pl.umap(adata, color = "Clusters", legend_loc = "on data", return_fig = True)
plt.savefig("../images/umap_unlabeled", dpi = 300, bbox_inches = "tight")
plt.close()

# ANNOTATION
def umap_stacked(adata):
    """
    Creates stacked image with plot corresponding
    to the expression level of certain marker
    @adata: scanpy object
    """
    n_cols = 4  # Define the number of plots per row
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
    plt.savefig("../images/umap_stacked.png", dpi = 300, bbox_inches = "tight")  # Save the figure
    plt.close()

def genes_dot(adata, n_genes = 5):
    """
    Creates top n_genes dotplot
    @adata: scanpy object
    @n_genes: top n genes to visualize
    """
    sc.tl.rank_genes_groups(adata, groupby = "Clusters", method = "wilcoxon")
    fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes = n_genes, return_fig = True)
    fig.savefig("../images/ranked_genes_dotplot.png", dpi = 300, bbox_inches = "tight")
    plt.close()

def map_cells(adata, marker_to_cell):
    """
    Performs manual annotation based on cells cluster location
    Functions `genes_dot` and `umap_stacked` should help in annotation
    @adata: scanpy object
    @marker_to_cell(dict): mapping from cell cluster # to prospective cell type
    """
    adata.obs["manual_celltype_annotation"] = adata.obs["Clusters"].map(marker_to_cell)
    sc.pl.umap(adata, color = ["manual_celltype_annotation"], return_fig = True)
    plt.savefig("../images/cell_types.png", dpi = 300, bbox_inches = "tight")
    plt.close()

# CONCORDANCE WITH 3D DATA
def spatial_cells(adata, coords):
    """
    Creates spatial DataFrame object with cell coordinates and cell clusters
    For `spatial_plot` function
    @adata: scanpy object
    """
    clusters = np.array(adata.obs["Clusters"]).astype(str)
    cells = adata.obs["manual_celltype_annotation"].astype('category')
    cells = cells.cat.add_categories(['Unknown'])

    cells.fillna('Unknown', inplace = True)
    cells = np.array(cells)

    coords["Clusters"] = clusters
    coords["Cells"] = cells
    coords['z'] = coords['z'] * 2 # every slice is 2 micrometers

    # Create a DataFrame with your coordinates and cluster information
    cell_df = pd.DataFrame({
        "Area": coords["area"],
        'X': coords['x'],
        'Y': coords['y'],
        'Z': coords['z'],
        'Clusters': coords['Clusters'],
        'Cells': coords['Cells']
    })
    cell_df.to_csv('../data/cell_coordinates.csv') # overwrites original file without clusters

    # Create a window for visualization in `spatial_plot` function
    scene = dict(
        xaxis = dict(range = [0, 900]),
        yaxis = dict(range = [0, 900]),
        zaxis = dict(range = [0, 200])
        )
    return cell_df, scene

def spatial_plot(adata, coords, accent, title, output):
    """
    Visualization either cell clusters or cell annotation in 3D
    @adata: scanpy object
    @coords(pd.DataFrame): contains cell coordinates
    @accent(str): either "Clusters" or "Cells"
    @title(str): description for cell clusters or types
    @output(str): path to the output image
    """
    cell_df, scene = spatial_cells(adata, coords)
    fig = px.scatter_3d(cell_df, x = 'X', y = 'Y', z = 'Z', color = accent, 
                        title = title, opacity = 0.7)
    fig.update_traces(marker = dict(size = 1))

    fig.update_layout(scene = scene)
    fig.write_image(output, scale = 2, engine = "kaleido")
    
# Run the functions you want
umap_stacked(adata)
genes_dot(adata, n_genes = 5)
map_cells(adata, marker_to_cell)
spatial_plot(adata, coords, 'Clusters', '3D Scatter - Cell Clusters', "../images/cell_clusters_3d.png")
spatial_plot(adata, coords, 'Cells', '3D Scatter - Cell Types', "../images/cell_types_3d.png")