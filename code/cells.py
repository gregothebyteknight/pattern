
# IMPORT MODULES
import numpy as np
import pandas as pd

coords = pd.read_csv("../data/cell_coordinates.csv", index_col = 0)
cell_type = [6] # some cell types are arranged between several clusters

def select_cell(coords, cell_type):
    cell_coords = coords[coords['Clusters'].isin(cell_type)]
    cell_coords_center = cell_coords.iloc[:,1:4] - cell_coords.iloc[:,1:4].mean()
    cell_coords_center["cluster"] = cell_coords["Clusters"]
    cell_coords_center["cell"] = cell_coords["Cells"]

    print("Mean by X coordinate:", cell_coords_center.iloc[:, 0].mean())
    print("Mean by Y coordinate:", cell_coords_center.iloc[:, 1].mean())
    print("Mean by Z coordinate:", cell_coords_center.iloc[:, 2].mean())
    return cell_coords_center

cell_coords_center = select_cell(coords, cell_type)
cell_coords_center.to_csv('../data/selected_cell_coordinates.csv')