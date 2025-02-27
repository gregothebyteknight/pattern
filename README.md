# Objective
Question of validity of using sections data as a representation of real tissue statistical properties. This validity is discussed and it is believed that 2D projections of the tissues do not provide with the accurate spatial distribution patterns

The **aim** of this pipeline is to compare spatial metrics both on true spatial data and on 2D sections, also to elaborate the method how to overcome the bias derived from two-dimensional sampling

# Theoretical Part
## Pair-Correlation Function
### Quick Introductions
Pair-correlation function $g(r)$ measures the probability of finding a pair of objects separated by a distance r relative to what would be expected for a completely random (Poisson) distribution
#### Derivation
Given that space is isotropic let suggest that probability of observing a point around location x is $$p(x)=\lambda(x)dx$$ where $\lambda(x)$ is intensity, or density of the point process (like Poisson process) and $dV$ - small region of volume

Let $p(x,y)$ be a probability of observing 1 point in x and 1 point in y, then $$p(x,y)=\lambda(x)dx\lambda(y)dy\bullet g(x,y)$$ where $g(x,y)=\frac{p(x,y)}{\lambda^2dxdy}$ - **pair-correlation function** (in simple case)

#### Intuition
Layers of spheres get more diffuse, so for large distances, the probability of finding two spheres with a given separation is essentially constant > a more dense system has more spheres, this it's more likely to find two of them with a given distance

The PC function accounts for these factors by normalizing by the density; this at large values of r it goes to 1, uniform probability

#### Calculation in 3D Case
1. Pick a value of `dr`
2. Loop over all values of r
    1. Count all particles that are a distance between `r` and `r+dr` away from the particle you're considering
    2. Divide your total count by N > average local density
    3. Divide this number by $4\pi r^2dr$ (volume of the spherical shell)
    4. Divide this by the particle number density $\rho$ - ensures that $g(r)=1$

For 2D volume correction will be just $2\pi r dr$

## Calculations
**Pipeline**: `reshape` > `express` > `cluster` > `choose` > `slices` > `main` (dependent of `slices`)
- **`reshape.py`** - takes stacks of images corresponding to the certain cell markers and cell masks ~ clusterization. This script creates cell datasets with coordinates `cell_coordinates.csv` and with intensities of certain cell markers - `cell_intensities.csv`
- **`express.r`** - creates expression matrices, both raw `expression_annotated.csv` and corrected `expression_annotated_corrected.csv` with compensation matrix `compensationMatrix.csv`
- **`cluster.py`** - performs clusterization and manual annotation of expression matrix `expression_annotated_corrected.csv`. Makes various plots for analysis
- **`choose.py`** - selects cell cluster obtain in the previous script > `cell_coordinates_clusters.csv`
- **`module.r`** - contains function to perform random sections based on the 3d dataset `cell_coordinates`. Compute **Pair-correlation function** for each of the slices
- **`slice.r`** - visualize both pcf of 3d initial dataset and of its 2d slices
