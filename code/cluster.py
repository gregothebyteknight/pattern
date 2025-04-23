
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
    "6": "Basal or myoepithelial (CK5, SMA, some Vimentin)",
    "19": "B (CD20)",
    "8": "Endothelial with mesenchymal features (vWF+CD31, Vimentin, some SMA, Collagen)",
    "5": "Smooth muscle (SMA, minimal Vimentin and Collagen)",
    "20": "ECM-producing stromal (CollagenI, lacking Vimentin, SMA)",
    "21": "Predominanly Cytotoxic T (CD45+, 60% CD8a)",
    "17": "Apoptotic macrophages (cPARP+Casp3, CD68, some CD45)",
    "24": "T (CD3)",
    "7": "Proliferating (Ki-67)",
    "4": "Mitotic (Phospho-H3)",
    "11": "Tumor epithelial (some Her2, panCK, CK19, CK7, CK8or18, some EorP Cadherin)",
    "3": "Tumor epithelial (some Her2, pS6, EorP Cadherin, CAIX)"
} # complete according to tsne stack and dotplot

rem_list = ['Hoechst0', 'Hoechst1', 'Hoechst2', 'Hoechst3', 'Hoechst4', 
                     'Hoechst5', 'Hoechst6', 'Hoechst7', 'Hoechst8', 'Hoechst9']

def tsne_image(adata):
    """
    Create tsne images for adat object
    @adata(scanpy object): object with expressions
    """
    sc.tl.tsne(adata)
    sc.pl.tsne(adata, color = "Clusters", legend_loc = "on data", return_fig = True)
    plt.savefig("../images/tsne_unlabeled", dpi = 300, bbox_inches = "tight")
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
    n_comps = 50 if adata.n_vars >= 50 else adata.n_vars - 1
    sc.pp.pca(adata, n_comps = n_comps, use_highly_variable = False, svd_solver = "auto")
    sc.pl.pca_variance_ratio(adata, log = True, show = False)
    plt.savefig("../images/pca_variance.png")
    plt.close()

    sc.pp.neighbors(adata, use_rep = 'X_pca', n_pcs = 16, n_neighbors = 30)
    print("Neighboring complete")
    sc.tl.leiden(adata, key_added = "Clusters", flavor = "igraph", resolution = 1.3)
    print(f"Number of clusters: {adata.obs['Clusters'].nunique()}")

    tsne_image(adata) # visualization of clustering

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

    tsne_image(adata_sample, "_sample") # visualization of clustering on sample data

    # CLUSTER TRANSFERING
    sc.tl.ingest(adata, adata_sample, obs = 'Clusters', embedding_method = 'pca')
    tsne_image(adata, "_complete") # visualization of clustering on whole dataset

# ANNOTATION
def tsne_stacked(adata):
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
        sc.pl.tsne(adata, color = tag, vmin = 0, vmax = "p99", sort_order = False,
                   frameon = False, cmap = "cividis",
                   show = False,  # Prevent immediate display
                   ax = axes[i])   # Assign to correct subplot

    # Hide unused subplots (if any)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig("../images/tsne_stacked.png", dpi = 300, bbox_inches = "tight")  # Save the figure
    plt.close()

def annotate(adata, n_genes = 3):
    """
    Helps with annotation by producing dot and violin plot
    @adata: scanpy object
    @n_genes: top n genes to visualize
    """
    sc.tl.rank_genes_groups(adata, groupby = "Clusters", method = "wilcoxon")

    fig = sc.pl.rank_genes_groups_dotplot(adata, n_genes = n_genes, return_fig = True)
    fig.savefig("../images/ranked_genes_dotplot.png", dpi = 300, bbox_inches = "tight")
    plt.close()

    fig = sc.pl.rank_genes_groups_stacked_violin(adata, groupby = 'Clusters', n_genes = n_genes, return_fig = True)
    fig.savefig("../images/ranked_genes_violin.png", dpi = 300, bbox_inches = "tight")
    plt.close()

def map_cells(adata, marker_to_cell):
    """
    Performs manual annotation based on cells cluster location
    Functions `genes_dot` and `tsne_stacked` should help in annotation
    @adata: scanpy object
    @marker_to_cell(dict): mapping from cell cluster # to prospective cell type
    """
    adata.obs["manual_celltype_annotation"] = adata.obs["Clusters"].map(marker_to_cell)
    sc.pl.tsne(adata, color = ["manual_celltype_annotation"], return_fig = True)
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
        xaxis = dict(range = [coords['x'].min(), coords['x'].max()]),
        yaxis = dict(range = [coords['y'].min(), coords['y'].max()]),
        zaxis = dict(range = [coords['z'].min() - 30, coords['z'].max() + 30])
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
    # expr = pd.read_csv('../data/expression_annotated_corrected.csv')

    # print("Shape of expression matrix:", expr.shape)
    # adata = preprocess(expr, rem_list)
    # clustering(adata)
 
    # tsne_stacked(adata)
    # annotate(adata, n_genes = 5)
    # adata.write("../data/adata.h5ad")

    # Uncomment after annotation
    adata = sc.read("../data/adata.h5ad")
    coords = pd.read_csv('../data/cell_coordinates.csv')
    
    map_cells(adata, marker_to_cell)
    spatial_plot(adata, coords, 'cluster', '3D Scatter - Cell Clusters', "../images/cell_clusters_3d.png")
    spatial_plot(adata, coords, 'cell', '3D Scatter - Cell Types', "../images/cell_types_3d.png")