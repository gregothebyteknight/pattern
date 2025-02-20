# Objective
The aim of this pipeline is to ...

# Theoretical Part
## Pair-Correlation Function
### Quick Introductions
Pair-correlation function $g(r)$ measures the probability of finding a pair of objects separated by a distance r relative to what would be expected for a completely random (Poisson) distribution
#### Derivation
Given that space is isotropic let suggest that probability of observing a point around location x is $$p(x)=\lambda(x)dx$$where $\lambda(x)$ is intensity, or density of the point process (like Poisson process) and $dV$ - small region of volume

Let $p(x,y)$ be a probability of observing 1 point in x and 1 point in y, then $$p(x,y)=\lambda(x)dx\lambda(y)dy\bullet g(x,y)$$where $g(x,y)=\frac{p(x,y)}{\lambda^2dxdy}$ - **pair-correlation function** (in simple case)

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

### `clusters.ipynb`


### `sections.rmd` 
**Input**: cell_specific matrix with coordinates in X,Y,Z
- dataset is previously centered by columns

**Output**: set of sections produced by applying consecutive section on the initial dataset in variable direction. Difference of orientation is achieved by performing matrix multiplication of initial dataset with the rotation matrices (roll, pinch and yaw)

WORKFLOW
1. Converting Dataframe to matrix with rows as cells and columns as coordinates (n_cells * 3 coordinates)
2. 