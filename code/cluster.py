
# IMPORT LIBRARIES
import numpy as np
import pandas as pd
import scanpy as sc
sc.settings.n_jobs = -1
import plotly.express as px
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from sklearn.preprocessing import StandardScaler

# VARIABLES
marker_to_cell = {
    "5": "Cytotoxic T cells (CD8a, CD3, CD45)",
    "18": "B cells (CD20)",
    "19": "Macrophages (CD68)",
    "12": "Basal layer epithelium (CK5, SMA, CK14)",
    "14": "Endothelium (vWF+CD31, Vimentin)",
    "8": "Carcinoma (CK19, Her2, CK7, panCK)",
    "4": "Potentially carcinoma (CK19, Her2, CK7, panCK)",
    "21": "Potentially neutrophiles (MPO, pS6, CD44)",
    "22": "Potentially neutrophiles (MPO, pS6, CD44)",
    "9": "Fibroblasts (CollagenI)",
    "16": "Potentially plasma cells (CD138, cRARP+cCasp3)",
} # complete according to umaps and dotplot

rem_list = ['Hoechst0', 'Hoechst1', 'Hoechst2', 'Hoechst3', 'Hoechst4', 
                     'Hoechst5', 'Hoechst6', 'Hoechst7', 'Hoechst8', 'Hoechst9']

def umap_image(adata, sup_text = ""):
    """
    Create umap images for 
    """
    sc.tl.umap(adata)
    sc.pl.umap(adata, color = "Clusters", legend_loc = "on data", return_fig = True)
    adata.write(f"../data/adata{sup_text}.h5ad")
    plt.savefig(f"../images/umap_unlabeled{sup_text}", dpi = 300, bbox_inches = "tight")
    plt.close()

# PREPROCESSING
def preprocess(expr, rem_list):
    """
    Performs scaling and cleaning of expression data
    @expr(pd.DataFrame): table with expressios
    @rem_list(list): list of variables to remove from adata
    """
    # SCALING
    scaler = StandardScaler()
    expr_scale = scaler.fit_transform(expr)

    # ADATA CONVERTION
    adata = sc.AnnData(X = expr_scale)
    adata.var_names = expr.columns
    print("Adata convertion complete")# Downloading the spatial centered data of specific cell type

    # ERASING USELESS VARIABLES
    adata = adata[:, ~adata.var_names.isin(rem_list)].copy()
    print(f"Removed {len(rem_list)} Hoechst markers. New shape: {adata.shape}")
    return adata

# CLUSTERING 
def clustering(adata):
    """
    Leiden clustering of adata based on expr data
    If dim(expr)[0] > 1000000 better use clustering_big
    @adata(scanpy object): object with expressions
    """
    # LEIDEN CLUSTERING 
    sc.pp.pca(adata, n_comps = 15, use_highly_variable = False) 
    sc.pp.neighbors(adata, use_rep = 'X_pca', n_pcs = 15)
    print("Neighboring complete")
    sc.tl.leiden(adata, key_added = "Clusters", flavor = "igraph")
    print(f"Number of clusters: {adata.obs['Clusters'].nunique()}")

    umap_image(adata) # visualization of clustering

def clustering_big(adata, sub_size = 100000):
    """
    Leiden clustering of big adata based on expr data
    @adata(scanpy object): object with expressions
    @sub_size(int): n_cells in sample to perform clusterization on
    """
    # VARIABLES
    mean_distances = np.zeros(adata.n_obs)

    # OBTAINING DENSITIES
    sc.pp.neighbors(adata, n_neighbors = 15, use_rep = 'X')
    distances = adata.obsp['distances']
    for i in range(adata.n_obs):
        neighbor_distances = distances[i].data
        mean_distances[i] = np.mean(np.sort(neighbor_distances)[1:15])
        if i % 100000 == 0:
            print(f"Percent neighbors computed {i / adata.n_obs}")
    # Unnormalized probabilities (reversed mean distance)
    prob_abs = 1 / (mean_distances + 1e-8)
    # Normalized probabilities
    probs = prob_abs / np.sum(prob_abs)

    # SAMPLING
    sample_ind = np.random.choice(adata.n_obs,
        size = sub_size, replace = False, p = probs)
    adata_sample = adata[sample_ind].copy()
    print("Subsample adata created")

    # CLUSTERIZATION ON SUBSAMPLE
    sc.pp.pca(adata_sample, n_comps = 15, use_highly_variable = False)  # Force PCA on all 30 markers
    sc.pp.neighbors(adata_sample, use_rep = 'X_pca')
    sc.tl.leiden(adata_sample, key_added = "Clusters")

    umap_image(adata_sample, "_sample") # visualization of clustering on sample data

    # CLUSTER TRANSFERING
    sc.tl.ingest(adata, adata_sample, obs = 'Clusters', embedding_method = 'pca')
    umap_image(adata, "_complete") # visualization of clustering on whole dataset

# ANNOTATION
def umap_stacked(adata):
    """
    Creates stacked image with plot corresponding
    to the expression level of certain marker
    @adata: scanpy object
    """
    n_cols = 4  # Define the number of plots per row
    n_plots = len(adata.var_names)
    n_rows = (n_plots // n_cols) + (n_plots % n_cols > 0)  # Calculate number of rows needed

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows))  # Adjust figure size

    # Flatten axes for easy iteration
    axes = axes.flatten()

    for i, tag in enumerate(adata.var_names):
        sc.pl.umap(adata, color = tag, vmin = 0, vmax = "p99", sort_order = False,
                   frameon = False, cmap = "#c46754",
                   show = False,  # Prevent immediate display
                   ax = axes[i])   # Assign to correct subplot

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
    adata.write("../data/adata_umap.h5ad")
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

    coords["cluster"] = clusters
    coords["cell"] = cells

    # Create a DataFrame with your coordinates and cluster information
    cell_df = pd.DataFrame({
        "area": coords["area"],
        'x': coords['x'],
        'y': coords['y'],
        'z': coords['z'],
        'cluster': coords['cluster'],
        'cell': coords['cell']
    })
    cell_df.to_csv('../data/cell_coordinates.csv', index = False) # overwrites original file without clusters

    # Create a window for visualization in `spatial_plot` function
    scene = dict(
        xaxis = dict(range = [0, 1400]),
        yaxis = dict(range = [0, 1400]),
        zaxis = dict(range = [0, 1400])
        )
    return cell_df, scene

def spatial_plot(adata, coords, accent, title, output, type = "-1"):
    """
    Visualization either cell clusters or cell annotation in 3D
    @adata: scanpy object
    @coords(pd.DataFrame): contains cell coordinates
    @accent(str): either "Clusters" or "Cells"
    @title(str): description for cell clusters or types
    @output(str): path to the output image
    """
    cell_df, scene = spatial_cells(adata, coords)
    if type != "-1":
        cell_df = cell_df[cell_df["cluster"] == type]
    fig = px.scatter_3d(cell_df, x = 'x', y = 'y', z = 'z', color = accent, 
                        title = title, opacity = 0.7)
    fig.update_traces(marker = dict(size = 1))

    fig.update_layout(scene = scene)
    fig.write_image(output, scale = 2, engine = "kaleido")

if __name__ == "__main__":
    # Comment after annotation
    expr = pd.read_csv('../data/expression_annotated_corrected.csv')
    coords = pd.read_csv('../data/cell_coordinates.csv')
    print("Shape of expression matrix:", expr.shape)
    adata = preprocess(expr, rem_list)
    clustering(adata)
    adata = sc.read("../data/adata_umap.h5ad")
    umap_stacked(adata)
    genes_dot(adata, n_genes = 5)

    # Uncomment after annotation
    # adata = sc.read("../data/adata_umap.h5ad")
    # map_cells(adata, marker_to_cell)
    # spatial_plot(adata, coords, 'cluster', '3D Scatter - Cell Clusters', "../images/cell_clusters_3d.png")
    # spatial_plot(adata, coords, 'cell', '3D Scatter - Cell Types', "../images/cell_types_3d.png")