
# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore")

import pandas as pd
import matplotlib.pyplot as plt

from pathlib import Path
from module import coord_align, scatter, map_slice, affine

# INITIALIZE VARIABLES
directory = Path('../data/cell_crc/')
csv_list = sorted([f for f in directory.iterdir() if f.is_file()])
dfs = [pd.read_csv(csv_list[0])] # Add reference slice (CRC1-002)
dfs[0]['z'] = 7.5 # set start z value

# Loop over the CSV files (only the first 25 if there are more)
for i, csv_path in enumerate(csv_list[1:], start = 1):
    df_t = dfs[i - 1] # reference slice
    df_s = affine(0, pd.read_csv(csv_path), df_t)
    df_s = coord_align(df_s, df_t)

    n_slice = int(csv_path.stem.split('-')[-1])
    df_s['z'] = n_slice * 5 - 2.5 if i <= 19 else (85 * 5) + (n_slice - 85) * 4 - 2
    dfs.append(df_s)
    print(f'{i}. slice is complete')

plt.close()

df_all = pd.concat(dfs, ignore_index = True)
df_all['x'] = df_all['x'].round(2)
df_all['y'] = df_all['y'].round(2)
df_all.iloc[:, -4:].to_csv('../data/cell_coordinates.csv', index = False)
df_all.iloc[:, :-4].to_csv('../data/cell_intensities.csv', index = False)

# Create a 5x5 grid of subplots with an appropriate figure size
fig, axes = plt.subplots(nrows = 5, ncols = 5, figsize = (15, 15), 
                         sharex = 'col', sharey = 'row')
layers = df_all['z'].unique()

for i, layer in enumerate(layers):
    title = "WD-76845-00" + str(map_slice.get(i, ''))
    scatter(df_all[df_all['z'] == layer], axes[i // 5, i % 5], title)

fig.supxlabel('X coordinate')
fig.supylabel('Y coordinate')
plt.tight_layout()
plt.savefig("../images/slice_align.png", dpi = 300, bbox_inches = "tight")
plt.close()