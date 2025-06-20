# `extract.py`
**Purpose:** Load 3D masks and multi‑channel TIFF stacks, extract cell features, and save to CSV

## Overview & Usage
### Input
+ `full_mask_final_segmentation.tif` (3D TIFF mask; shape: [z, y, x])
+ A set of channel directories under `../data/init/`, each containing 2D TIFF slices named in acquisition order

### Output
+ `cell_coordinates.csv` (CSV with columns: area, x, y, z)
+ `cell_intensities.csv` (CSV with mean intensity per channel; headers from folder names)

## Pipeline Steps & Functions
### Setup & Imports
Load libraries and suppress warnings

### Mask Loading & Reordering
`imread` imports 3D mask; `np.swapaxes` reorders to [x, y, z]

### Channel Discovery & Preallocation
Identify subdirectories as channels; define image array of shape [x, y, z, n_channels]

### Reshape image `reshape_image(channels)`
Reads each per-channel TIFF stack into the image array and returns multi-channel volume

### Feature Extraction
Compute cell properties:

+ mean_intensity across channels via `measure.regionprops_table(masks, intensity_image = image)`
+ area and centroid (scaled for z-spacing) via `measure.regionprops_table(masks)`

### Save Results
Export Pandas DataFrames to CSV

## How to Run
`python extract.py`

+ Ensure relative paths are correct, and that TIFF filenames are zero-padded for sorted() accuracy
+ Adjust voxel-size comment or metadata reading if z-spacing differs from 2 µm. (Parquet/HDF5) for large CSVs. versions
+ Adjust TIFF dimensions with the `x_dim, y_dim, z_dim` 
---
# `correct.r` 
**Purpose**: Compensate IMC intensities and filter out technical channels

## Overview & Usage
### Input
+ `cell_intensities.csv` (CSV with mean intensities per cell × channel)
+ `tags_to_channels.csv` (mapping of metal tags to cleaned channel names). Will be loaded to `comp_mat`
+ `compensation_matrix.csv` (raw spillover matrix; rows/columns = metal tags)

### Output
+ `heat_unc.png` & `heat_cor.png` (correlation heatmaps before / after compensation)
+ `expression_annotated_corrected.csv` (compensated intensities, filtered)

## Pipeline Steps & Functions
### Setup & Imports
`suppressPackageStartupMessages()` hides startup noise; loads `pheatmap`

### Data Loading & Renaming
Read expression and panel tables; relabel `expr` columns using `panel$dollar_sign$clean_target`

### Uncompensated Heatmap
Plot correlation of raw `expr` with Ward clustering to `heat_unc.png`

### Compensation Matrix Preparation
Load `comp_mat`; filter to channels present in data; expand with zeros for missing tags; reorder to full channel list; set diagonal = 1

### Apply Compensation
Invert `comp_mat` via `solve()`; right-multiply raw matrix to get `expr_corr`; relabel columns

### Post-Compensation Heatmap
Plot correlation of `expr_corr` to `heat_cor.png`

### Optional Cleaning
Remove non-detection channels (like Ir193, Ir191, etc.) from `expr_corr_filt`

### Save Results
Write `expr_corr_filt` to `expression_annotated_corrected.csv`

## How to Run
`Rscript correct.r`

+ Ensure working directory contains `../data/` and `../images/ subfolders`
+ Verify column names in panel match `expr` metal tags exactly

---

# `cluster.py`
**Purpose:** cells clusterization and annotation

## Overview & Usage
### Inputs:
+ expression_annotated_corrected.csv (raw or scaled expression matrix)
+ cell_coordinates.csv (area, x, y, z)

### Outputs:
+ PNG figures in ../images/ (t‑SNE plots, PCA variance, dot/violin plots, stacked t‑SNE, 3D scatter of clusters and cell types)
+ Updated adata.h5ad with clustering and manual annotations
+ Overwritten cell_coordinates.csv augmented with cluster and cell labels

## Pipeline Steps & Key Functions
### Setup & Imports
NumPy, Pandas, Scanpy (`n_jobs = -1`), Plotly Express, Matplotlib, StandardScaler

+ Warnings suppressed to reduce console clutter

### Marker Configuration (Variables)
+ `marker_to_cell`: maps Leiden cluster IDs (“6”, “19”, etc.) → human‑readable cell‑type names
+ `rem_list`: list of Hoechst channels to drop before clustering

### t-SNE Utility
Runs t‑SNE, colors by “Clusters,” saves high‑res PNG

### Preprocessing
Scales data, constructs AnnData, removes unwanted channels

### Clustering
#### Standard (`clustering`)
+ PCA (up to 50 PCs), neighbor graph (16 PCs, 30 neighbors), Leiden (resolution 1.3)
+ Saves `pca_variance.png` and invokes `tsne_image()`

#### Large dataset (`clustering_big`):
+ Computes per‑cell mean neighbor distance → sampling probabilities
+ Subsamples ~100 k cells, clusters sample, then ingests labels into full AnnData
+ Generates t‑SNE for both sample and full sets

### Annotation Helpers
Complete annotation (cluster number: cell type name) should be put in market_to_cell

+ `tsne_stacked()`: grid of marker‑based t‑SNE plots (`tsne_stacked.png`)
+ `annotate()`: ranked‑genes dotplot + stacked violin plots saved as PNGs
+ `map_cells()`: maps clusters → cell types via `marker_to_cell`, plots colored t‑SNE

### Spatial Integration & 3D Plotting
+ `spatial_cells()`: merges cluster/cell labels into coordinate table; returns DataFrame + Plotly scene config
+ `spatial_plot()`: filters (by cluster, if requested), then makes a 3D scatter with Plotly and saves via Kaleido

### Main Block
+ Offers both raw‑from‑CSV pipeline (commented) and quick-load of precomputed adata.h5ad
+ Executes map_cells() and two spatial_plot() calls for clusters and cell types

## How to Run
`python cluster.py`

+ If starting from scratch, uncomment the preprocessing + clustering block and comment out the adata.h5ad loader
+ Ensure your working directory contains the required CSVs and that Plotly’s Kaleido engine is installed

---

# `choose.py`
## Overview & Usage
### Input
+ `cell_coordinates.csv` (with at least x, y, z, cluster, cell columns). Parsed with `--path_to_coords` or `-p` flag

#### Parameters
+ `--clusters`, `-c`: list of integer cluster IDs to select (default -1 → all clusters)
+ `--sample_size`, `-s`: number of cells to randomly sample per selection (optional)

### Output:
- `cell_coordinates_selected.csv` in the same folder, with centered coordinates and retained annotations

## Pipeline Steps & Key Logic
### Argument Parsing
+ Builds CLI with defaults; constructs out_path by appending _selected to the filename.

### Data Loading & Default Cluster Handling
+ Reads coordinate table; -1 flag means “select all clusters”

### Selection function (`select_cell()`)
+ Filtering, subsampling, and mean-centering of the first three columns
+ Preserves cluster (and optionally cell) annotations

### Execution & Saving
+ Saves the centered subset to `<original>_selected.csv`

## How to Run
```bash
python choose.py \
  --path_to_coords ../data/cell_coordinates.csv \
  --clusters 3 7 11 \
  --sample_size 50000
```
+ Omit `--clusters` to select all, omit `--sample_size` for no subsampling

---

# `prespace.r`
**Purpose**: cumulative pcf plot generation, angle analysis, compute Clark-Evans index for further analysis

## Overview & Usage
### Input
+ CSV of centered cell coordinates (columns: x, y, z, cluster, cell) passed as first command line argument
+ Optional second and third args define a valid range of cell counts per slice

### Output
+ In‐memory tables: `ce_tbl` (Clark–Evans indices) and `pcf_naive` (max PCF per slice plus true 3D)
+ PNG figure from `cumul_pcf()` (PCF curves)
+ Logs (`../logs/logs_angle_<date>.txt`) with session info
+ CSVs for `ce_tbl` and `pcf_naive`
+ Optional slice‐range analysis via `angle_analysis()`

### Dependencies
+ `module.r`
    - `pc_for_slice` - cell cluster rotation and dissection process
    - `angle_analysis` - select dissection within a range of cell counts
    - `ce_idx` - Clark-Evans index computation
+ `visual.r`
	- `cumul_pcf` - cumulative PCF visualization
    
## Pipeline Steps & Key Logic
### Setup & Imports
Defines working directory, loads custom modules, suppresses startup messages

### Argument Parsing
Reads the path, loads coords; sets `n_range` for filtering valid slices

### True 3D PCF & CE
Computes 3D point‐pattern object and its isotropic pair‐correlation estimate

### Slice‐by‐Slice PCF Computation
Repeatedly rotates the 3D cloud, projects to 2D slices via `pc_for_slice()`, and records PCF and Clark–Evans indices for each

### Aggregate True & Naïve Metrics
Appends the true‐3D values for direct comparison

### Plotting & Logging
Generates cumulative PCF figure; saves a log with package versions

### Optional Outputs & Analysis
+ CSV writes for `ce_tbl` and `pcf_naive`
+ If `n_range` is set, runs `angle_analysis()` on slices within the valid range
+ Prints “Optimal max_r” based on slice with most cells

## How to Run
```bash
Rscript angle.r ../data/cell_coordinates_selected.csv min_cells max_cells
```
+ The second/third arguments are optional; omit to process all slices without range filtering

---

# `space.r`
**Purpose**: mean PCF 2D vs 3D visualization, save PCF values for further analysis

## Overview & Usage
### Inputs
+ Path to centered cell coordinates CSV (x, y, z, cluster, cell)
+ Maximum radius `r_max` for PCF estimation

### Outputs
+ PCF visualization via `mean_pcf()` (in `visual.r`).
+ Log file: `../logs/logs_slice_<date>.txt`
+ CSVs for slice PCF (`pcf_plane.csv`), space PCF (`pcf_space.csv`), and counts (`pcf_num.csv`)

### Dependencies
+ `module.r`
    - `pc_for_slice` - cell cluster rotation and dissection process
+ `visual.r`
	- `mean_pcf` - mean PCF visualization

## Pipeline Steps & Key Logic
### Setup & Imports
Defines project root, loads custom modules and required packages

### Argument Parsing & Variables
Reads coordinates, constructs an equally spaced radius vector (`r_grid`), and extracts `cell_type`

### True 3D PCF
Builds a 3D point‐pattern object and computes its isotropic pair correlation function

### Slice‐wise PCF Computation
Randomly rotates the 3D cloud to generate 2D slices via `pc_for_slice()`, storing PCF curves (`g`) and slice cell counts

### Visualization
Calls the plotting routine in `visual.r` to compare the true 3D PCF against slice distributions

### Archiving & Data Preparation
Logs session info, prepares data frames for downstream CSV export, and ensures target directories exist

## How to Run
```bash
Rscript slice.r path/to/cell_coordinates.csv 50
```
+ Replace 50 with your desired maximum radius (use max_r from `prespace` or based on the cumulative PCF plot)
+ To export CSVs, uncomment the `write.csv()` lines at the end

---

# `reveal.r`
**Purpose**: visualize the difference of spatial statistics between 2D and 3D as a box plot

## Overview & Usage
### Inputs
+ Folders under `../data/pcf/{mode}/{dataset}/{cell}` containing: `pcf_plane.csv`, `pcf_space.csv`, and `pcf_num.csv` (or `ce_tbl.csv/pcf_naive.csv` for "naive" mode).

#### Parameters (hard‑coded at bottom):
+ mode = "pcf" or "naive"
+ type = "gamma" or "oscillation" (or "ce"/"naive" if mode == "naive")
+ sim = FALSE (skip simulation datasets)
+ g = TRUE (apply transformed parameters for PCF fits)

### Outputs
+ `temp/tbl_all.csv`: aggregated fit parameters for every (`data`, `cell`, `dim`)
+ `temp/pval_tbl.csv`: per‐cell two‐sided p‐values (and z‐scores) comparing space vs. plane
+ Final visualization via `plot_space()`, saved according to `labels$\dollar_sign$file` 

## Pipeline Steps & Key Logic
### Setup & Imports
Establish project root, load plotting & fitting modules, and data‐manipulation tools

#### Configuration Objects
+ `fallback`: default fit parameters for “gamma” and “oscillation.”
+ `labs`: nested lists of titles, file names, and y‐axis labels for each mode/type

### Fit Functions

#### `pcf_mode()`
Reads `pcf_{dim}.csv` and `pcf_num.csv`, loops over slices, fits the chosen model (fit_func), computes $R^2$.

If `g = TRUE`, transforms gamma‐fit parameters to "amplitude" and "range"

Returns a tibble of (`data`, `cell`, `dim`, `par_1`, `par_2`, `r2`, `n`)

### `naive_mode()`
Reads either ce_tbl.csv or pcf_naive.csv, filters out small slices, and returns (`data`, `cell`, `dim`, `par_1`, `par_2`)

### Data Aggregation: `reveal_all()`
Iterates over `../data/{mode}` subdirectories, skips simulation runs if `sim = FALSE`, and calls the appropriate mode function for both "plane" and "space" dims

Filters out non‐finite parameters before saving

### P‑Value Computation
#### `reveal_pval()`
Mid‐p test directly comparing lists of plane‐values vs. single space‐value per (data, cell)

#### `reveal_pval_z()`
Converts plane‐values to z‐scores relative to plane median & SD; computes mid‐p via z; adjusts for multiple testing (Bonferroni, Holm, BH, BY)

### Visualization: `plot_space()`
Uses labs metadata to generate a comparative plot of plane vs. space parameters, annotated with p‑values

## How to Run
```bash
Rscript analysis_pipeline.r
```
+ Adjust the hard‑coded flags (`mode`, `type`, `sim`, `g`, `par`) at the bottom of the script to switch between modes, fitting options, and which parameter to test

---

# `evaluate.r`
**Purpose**: visualize computed p-values of spatial difference as a bar plot

## Overview & Usage
### Input
+ `../temp/pval_tbl.csv` containing columns: `data`, `cell`, `p.BH` (adjusted p‑value), and `n` (number of slices)

### Output
+ Bar plot saved as PDF (`../images/FINALS/Pvalue/pvalue.pdf`), comparing $–\log_{10}FDR‑BH$ across cell types and datasets

## Pipeline Steps & Key Logic
### Setup & Imports
Ensures working directory is project root; loads plotting and data‑manipulation packages; checks for Cairo capability

### Data Preparation
Reads p‑value table, computes $–\log_{10}$ of adjusted p‑values (`p_log`), and casts cell as factor

### Color & Label Configuration
Defines a manual color palette and corresponding legend labels for each dataset

### Plot Construction
Horizontal bar chart (`geom_col`), grouped by data, with a dashed threshold line at $–\log_{10}0.05 \approx 1.3$

### Export to PDF via Cairo
Creates output directory, opens a Cairo PDF device, renders the plot, and closes the device

### How to Run
```bash
Rscript plot_pvalues.r
```
+ Ensure `pval_tbl.csv` exists under `../temp/`
+ Cairo support is required for PDF output with embedded fonts

---

# `heat.r`
**Purpose**: visualize short- and long-term differences in a single heat plot

## Overview & Usage
### Input
+ Per-test p‑value tables in `../temp/{stat_type}/pval_tbl.csv`, each with columns including `data`, `cell`, and `p`

### Output
+ Combined CSV `../temp/pval_all.csv` with $–\log_{10}p$ values per test, and a clustered heatmap PNG at `../images/FINALS/heatmap/pval_heat.png`

## Pipeline Steps & Key Logic
### Setup & Imports
Sets project root; loads `pheatmap` and defines the temp directory

### Read Base Table
Imports the initial “space vs. plane” p‑value summary for each (`data`, `cell`)

### Aggregate $–\log_{10}p$ Across Tests
Iterates over each subdirectory in `../temp/`, reads its p‑values, replaces zeros with a small floor, computes $–\log_{10}p$, and appends as a new column

### Prepare Heatmap Matrix
Constructs a numeric matrix (`mat`) of size (samples × tests), labels rows by dataset|cell, and defines a row annotation for dataset grouping

### Plot Heatmap
Generates a clustered heatmap of $–\log_{10}p$ values, with rows annotated by dataset, and saves it as a PNG

## How to Run
```bash
Rscript plot_pval_heatmap.r
```
+ Ensure that each subdirectory of `../temp/` contains its own `pval_tbl.csv`
+ Confirm that `../temp/pval_tbl.csv` (the master comparison table) exists at the start
