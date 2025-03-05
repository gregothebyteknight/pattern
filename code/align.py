
# IMPORT LIBRARIES
import os
import torch
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from STalign import STalign
from pathlib import Path as path
from contextlib import redirect_stdout # to supress output

# INITIALIZE VARIABLES
directory = path('../data/cell_crc/')
csv_list = sorted([f for f in directory.iterdir() if f.is_file()])
dfs = [pd.read_csv(csv_list[0])] # Add reference slice (CRC1-002)
dfs[0]['z'] = 7.5 # set start z value

def coord_align(df_s, df_t):
    """
    Shift coordinates of sample dataset towards target
    using STalign approach
    """
    x_s = np.array(df_s['x']) # s is for source
    y_s = np.array(df_s['y']) # source - subject to alignment
    x_t = np.array(df_t['x']) # t is for target
    y_t = np.array(df_t['y']) # target - reference

    # RASTERIZE
    with open(os.devnull, 'w') as fnull, redirect_stdout(fnull):
        x_s_rast, y_s_rast, s_tensor, _ = STalign.rasterize(x_s, y_s ,dx = 30, blur = 1.5)
        x_t_rast, y_t_rast, t_tensor, _ = STalign.rasterize(x_t, y_t, dx = 30, blur = 1.5)

    # SOLVE MAPPING
    if torch.cuda.is_available():
        device = 'cuda' # or cuda:i if many
    else:
        device = 'cpu'
    torch.set_default_device(device) # for point mapping

    params = {'niter': 10000, 'device': device, 'epV': 50}

    out = STalign.LDDMM([y_s_rast, x_s_rast], s_tensor, [y_t_rast, x_t_rast], t_tensor, **params)
    plt.close('all')

    a = out['A']
    v = out['v']
    xv = out['xv']

    # APPLY MAPPING
    points = STalign.transform_points_source_to_target(xv, v, a, np.stack([y_s, x_s], 1))

    # switch tensor from cuda to cpu for plotting with numpy
    if points.is_cuda:
        points = points.cpu()

    x_lddmm = points[:, 1]
    y_lddmm = points[:, 0]

    df_s['x'] = x_lddmm
    df_s['y'] = y_lddmm

    return df_s

# Create a 5x5 grid of subplots with an appropriate figure size
fig, axes = plt.subplots(nrows = 5, ncols = 5, figsize = (15, 15))
ax = axes[0, 0]
sc = ax.scatter(x = dfs[0]['x'], y = dfs[0]['y'], s = 0.5)
ax.set_title(csv_list[0].stem, fontsize = 8)
ax.set_xticks([])
ax.set_yticks([])

# Loop over the CSV files (only the first 25 if there are more)
for i, csv_path in enumerate(csv_list[1:], start = 1):
    # Determine the row and column indices
    ax = axes[i // 5, i % 5]
    
    # Read CSV file into a DataFrame
    df_t = dfs[i - 1] # reference slice
    df_s = coord_align(pd.read_csv(csv_path), df_t)
    n_slice = int(csv_path.stem.split('-')[-1])
    
    if i <= 19:
        df_s['z'] = n_slice * 5 - 2.5
    else:
        df_s['z'] = (85 * 5) + (n_slice - 85) * 4 - 2
    
    dfs.append(df_s)
    
    # Scatter plot using 'AREA' column for the color
    sc = ax.scatter(x = df_s['x'], y = df_s['y'], s = 0.5)
    ax.set_title(csv_path.stem, fontsize = 8)
    ax.set_xticks([])
    ax.set_yticks([])
    print('/\/\/\/\/\/\/\/\/\/\/\/\/\/\/')

plt.tight_layout()
plt.savefig("../images/slice_align.png", dpi = 300, bbox_inches = "tight")
plt.close()

df_all = pd.concat(dfs, ignore_index = True)
df_all.to_csv('../data/cell_coordinates.csv', index = False)
df_all.to_csv('../data/cell_intensities.csv', index = False)
