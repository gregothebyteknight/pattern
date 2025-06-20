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

#