
# IMPORT MODULES
import pandas as pd
import argparse

coords = pd.read_csv("../data/cell_loc_complete_cluster.csv", index_col = 0)

parser = argparse.ArgumentParser(
                    prog = 'choose.py',
                    description = 'Choose the cluster of cells ')

parser.add_argument('--clusters', '-c', type = int, nargs = '+', default = [-1])
parser.add_argument('--sample_size', '-s', type = int, default = None)
args = parser.parse_args()
cluster_type = args.clusters
n_sample = args.sample_size

if cluster_type[0] == -1:
    cluster_type = coords["cluster"].unique()

def select_cell(coords, cluster_type, n_sample):
    # Selecting cell type
    cell_coords = coords[coords['cluster'].isin(cluster_type)]
    print("Number of cells in cell type before:", cell_coords.shape[0])
    # Subsampling
    if n_sample is not None:
        if cell_coords.shape[0] >= n_sample:
            cell_coords = cell_coords.sample(n = n_sample)
        else: print("Less than 50,000 cells available; using all cells.")
    # Centering
    cell_coords_center = cell_coords.iloc[:, 1:4] - cell_coords.iloc[:, 1:4].mean()
    cell_coords_center["cluster"] = cell_coords["cluster"]
    cell_coords_center["cell"] = cell_coords["cell"]

    print("Selected cell types:", *cell_coords_center["cell"].unique())
    print("Mean by X coordinate:", cell_coords_center.iloc[:, 0].mean())
    print("Mean by Y coordinate:", cell_coords_center.iloc[:, 1].mean())
    print("Mean by Z coordinate:", cell_coords_center.iloc[:, 2].mean())
    print("Number of cells in cell type after:", cell_coords_center.shape[0])
    return cell_coords_center

cell_coords_center = select_cell(coords, cluster_type, n_sample)
cell_coords_center.to_csv('../data/selected_cell_coordinates.csv')