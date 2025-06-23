# `align.py`
**Purpose**: align multiplexed imaging slices into 4D dataframe 

## Overview & Usage
### Inputs
+ A directory of per‐slice CSVs in `../data/cell_crc/`, each with columns `x`, `y`, `area`, `cluster`, `cell` (or similar)

### Outputs
+ `../data/cell_coordinates_1.csv`: concatenated, aligned XYZ coordinates per cell
+ `../data/expression_annotated_corrected_1.csv`: concatenated expression data across all slices
+ `../images/slice_align_1.png`: 5×5 grid of scatter plots by Z‐layer

### Dependencies
+ `module.py`
    - `affine` - affine transformations
    - `coord_align` - affine and non-affine transformations with STalign
    - `scatter` - maps plot to the panel
    - `map_slice` - dictionary with mapping order: slice number

## Pipeline Steps & Key Logic
### Setup & Imports
Suppresses warnings, loads core libraries and custom alignment/plotting functions

### Initialize Reference Slice
Reads the first CSV as the reference (`CRC1‑002`), assigns a fixed Z

### Iterative Alignment Loop
Aligns each subsequent slice to its predecessor, computes its Z position, and appends to the list

### Concatenate & Export Tables
Merges all aligned slices, rounds XY, and writes two CSVs: one for coordinates/annotations, one for expression

### Visualization: 5×5 Grid of Slices
Creates a 5×5 panel of 2D scatter plots for each unique Z‐layer, saving a high‐res PNG

## How to Run
```bash
python slice_alignment.py
```
Ensure `../data/cell_crc/` contains your slice CSVs and that `module.py` with `coord_align`, `affine`, and `scatter` is available

# Gamma fitting `Mimic.r`
**Purpose**: Fit either a Gamma‐shaped or oscillatory model to empirical pair‐correlation data (PCF) via nonlinear least squares

## Overview & Usage
### Inputs
+ Vectors `x` (distance) and `y` (PCF values)
+ `type = "gamma"` or `"oscillation"`
+ `fallback`: named list of default parameters for each model
+ `show_plot`: if TRUE, display a diagnostic plot

### Output
+ Named vector of fitted parameters (`height`, `shape`, `scale` for Gamma; `amp`, `decay`, `rate` for oscillation)

## Pipeline Steps & Key Logic
```r
fit_func(x, y, type, fallback, show_plot)
```

### Model Library
Defines two parametric forms: a shifted Gamma PDF and a damped sine wave

### Data Prep
```r
df <- data.frame(x, y) |> filter(is.finite(x), is.finite(y))
stop("No data") if empty
```

### Initial Guess:
+ **Gamma**: uses weighted mean / variance of `y − 1` to derive `shape_init` and `scale_init`; `height = max(y − 1)`
+ **Oscillation**: sets `amp`, `decay = 1`, and `rate` from the location of the first peak (or `fallback` if invalid)

### Fallback Check:
If fit fails or returns non‐finite parameters, retries with `fallback[[type]]`

### Plotting (Optional)
Overlays data points and fitted curve, annotates parameter values

## How to Call
```r
source("model_fitting.r")

# Define PCF data vectors x and y, and a fallback list:
fallback <- list(
  gamma = c(height = 15, shape = 3, scale = 4),
  oscillation = c(amp = 15, decay = 1, rate = 4)
)

# Fit a Gamma model
pars_gamma <- fit_func(x, y, type = "gamma", fallback = fallback, show_plot = TRUE)

# Fit an oscillatory model
pars_osc <- fit_func(x, y, type = "oscillation", fallback = fallback)
```

# Additional module in Python `Module.py`
**Purpose**: Provide functions for aligning 2D cell-coordinate slices via affine transforms or STalign diffeomorphic mapping, and for plotting aligned points

## Overview & Usage
### Dependencies:
+ `STalign` (for rasterization and LDDMM mapping)
+ `torch`, `numpy`, `scipy.stats.trim_mean`
+ `matplotlib` for plotting

### Key Outputs
+ Aligned DataFrame with updated x, y columns
+ Scatter plots on user-provided Matplotlib axes

## `affine(deg, df_s, df_t = None)`

### Inputs
+ `deg` (int): rotation angle in degrees (clockwise)
+ `df_s` (pd.DataFrame): source coordinates (x, y)
+ `df_t` (pd.DataFrame, optional): target coordinates for scaling & translation

### Process:
Rotate points by −deg via a 2×2 rotation matrix

If df_t given, scale source to match target standard deviations, then translate using 10% trimmed means to align centroids

### Output
Returns df_s with transformed x, y

## `coord_align(df_s, df_t)` 
### Inputs
Two DataFrames of x, y coordinates (source and target)

### Process
Rasterize both point sets into blurred count‐tensors via STalign.rasterize (silencing stdout)

Perform LDDMM mapping (STalign.LDDMM) to compute velocity fields v and transformation A

Apply transform to source points (STalign.transform_points_source_to_target) and update df_s

### Output
Returns df_s with diffeomorphically aligned x, y

## `scatter(df, ax, title)`
### Inputs
+ `df` (pd.DataFrame): must include x, y, area
+ `ax` (matplotlib Axes): subplot in which to draw
+ `title` (str): subplot title

### Process
Plots points colored by area at size 0.5, and sets the title with small font

### Output
None (draws on ax)

## How to Use
```py
from alignment_utils import affine, coord_align, scatter
import pandas as pd

# Load two slices
df_src = pd.read_csv("slice1.csv")
df_tgt = pd.read_csv("slice2.csv")

# Affine align (rotate + scale + translate)
df_aligned = affine(deg=10, df_s=df_src, df_t=df_tgt)

# Or perform diffeomorphic alignment
df_lddmm = coord_align(df_src, df_tgt)

# Plot onto existing Matplotlib axes
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
scatter(df_lddmm, ax, title="Aligned Slice")
plt.show()
```

# Additional module in R `module.r`
**Purpose**: Provide core spatial‐statistics utilities: defining point‐pattern windows, 3D rotations & slicing, per‐slice PCF & Clark–Evans computations, and selecting an optimal slice

## Overview & Usage
### Inputs
+ `cell_mat`: an n×3n×3 numeric matrix of 3D coordinates
+ Rotation angles (a, b, g), slice thickness, and valid angle pairs

### Outputs
+ `pc_for_slice()`: list with pcf object, num_cells, and ce index for a given orientation
+ `angle_analysis()`: prints slice size and saves a PNG of the “best” slice by angle
+ `ce_idx()`: numeric Clark–Evans index for 2D or 3D point sets

## Function Breakdown
### Window Builder: `pp_box(mat)`
Inspects mat columns; returns either a 2D owin() or 3D box3() spanning the data ranges

### Rotation Matrix: `rotation(a,b,g)`
Yaw (a), pitch (b), roll (g) → 3×3 rotation matrix applied as R xRx

### Slice Selector: `slice(a,b,c,slice_size,mut_frame)`
Given plane normal (a,b,c)(a,b,c) and half‐thickness, returns logical mask of points within the slab

### Per‐Slice Analysis: pc_for_slice(a,b,g,cell_mat,r_grid=NULL)
Rotates cell_mat by R(a,b,g)R(a,b,g), filters a central 10-unit slice, and:
+ If ≤ 1 point, returns NA-filled PCF of length 513.
+ Else, builds 2D ppp object, computes isotropic `Kest()` and its `pcf.fv()`, and calls `ce_idx()`

**Returns**: list(pcf = …, num_cells = …, ce = …)

### Optimal Slice Plotter: `angle_analysis(cell_mat, valid_pairs, type)`
From valid_pairs (roll, pinch), selects the pair whose sum is closest to the mean sum

Produces and saves a PNG of the 10-unit thick slice in that orientation, reporting its cell count

### Clark–Evans Index: `ce_idx(df)`
+ **2D**: wraps clarkevans() on a ppp object (needs ≥4 points)
+ **3D**: approximates expected nearest‐neighbor distance from convex‐hull volume & surface area, then computes observed mean nearest-neighbor via `FNN::get.knnx()`, and returns robs/rexprobs​/rexp​

## How to Use
```r
source("spatial_utils.R")

# Compute PCF for a rotated slice
res <- pc_for_slice(a=0, b=pi/4, g=pi/6, cell_mat, r_grid=seq(0,10,length=513))

# Perform slice‐angle analysis and save the best slice plot
angle_analysis(cell_mat, valid_pairs=your_pairs, type="your_cell_type")

# Directly compute Clark–Evans for arbitrary points
ce_value <- ce_idx(data.frame(X=...,Y=...,Z=...))
```

# Automation `multiply.sh`
**Purpose**: Automate running angle.r and slice.r across 30 simulation folders, extracting an optimal radius from each run and scaling it before the slice analysis

## Overview & Usage
### Inputs
+ Base directory: `../../spheres/data/sim_<i>/`, each containing one CSV of cell coordinates
+ R scripts: `angle.r` (computes optimal max_r) and `slice.r` (requires r_max)

### Outputs
+ Console logs for each simulation, including parsed and scaled max_r
+ angle.r and slice.r outputs stored as defined in those scripts

## Pipeline Steps & Key Logic
### Setup & Loop
Iterates simulations 1–30, skipping any missing directories

### Input File Detection
Picks the first CSV in each folder or skips if none found.

### Run `angle.r` & Parse `max_r`
Captures angle.r stdout, greps the “Optimal max_r” line, and extracts the numeric value

### Scale & Validate
Multiplies by 1.1 and ensures a valid result before proceeding

### Run slice.r with Scaled Radius
Executes slice analysis using the scaled radius

## How to Run
```bash
chmod +x multiply.sh
./multiply.sh
```
+ Ensure the script has execute permissions and that `angle.r` and `slice.r` paths are correct
+ Adjust `BASE_DIR` or the simulation range (`seq 1 30`) as needed