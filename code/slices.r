# Installing the libraries
library(tidyverse)
library(spatstat)
library(magrittr)

# Define the ranges for each coordinate
pp_box <- function(mat) {
  "Define the range of the point pattern object
  mat(matrix): the 2D or 3D point pattern matrix
  "
  x_range <- c(min(mat[, 1]), max(mat[, 1]))
  y_range <- c(min(mat[, 2]), max(mat[, 2]))

  if (dim(mat)[2] == 3) {
    z_range <- c(min(mat[, 3]), max(mat[, 3]))
    # Create the scale frame for point pattern (either 3D or 2D)
    spatstat.geom::box3(xrange = x_range, yrange = y_range, zrange = z_range)
  } else {
    spatstat.geom::owin(xrange = x_range, yrange = y_range)
  }
}

# Declare the rotation matrices
rotation <- function(a, b, g) {
  "Compute the rotation matrix
   a(alpha of yaw): angle of rotation about the z-axis
   b(beta of pinch): angle of rotation about the y-axis
   g(gamma of roll): angle of rotation about the x-axis"
  matrix(c(
           cos(a) * cos(b), # first row
           cos(a) * sin(b) * sin(g) - sin(a) * cos(g),
           cos(a) * sin(b) * cos(g) + sin(a) * sin(g),
           sin(a) * cos(b), # second row
           sin(a) * sin(b) * sin(g) + cos(a) * cos(g),
           sin(a) * sin(b) * cos(g) - cos(a) * sin(g),
           -sin(b), # third row
           cos(b) * sin(g),
           cos(b) * cos(g)),
  nrow = 3, byrow = TRUE)
}

# Declare the slice function
slice <- function(a, b, c, slice_size, mut_frame) {
  "Compute the slice of the 3D point pattern.
   a, b, c: the normal vector of the slice (d = 0)
   slice_size: the depth of the slice
   mut_frame: the 3D point dataframe undergo transformation"
  value <- a * mut_frame[, 1] + b * mut_frame[, 2] + c * mut_frame[, 3]
  # Return a logical vector: TRUE if value is within the slice, FALSE otherwise
  value >= -slice_size / 2 & value <= slice_size / 2
}

# Define the function to compute the dissections
pc_for_slice <- function(a, b, g, cell_mat) {
  "Compute the dissections of the 3D point pattern
   a(alpha of yaw): angle of rotation about the z-axis
   b(beta of pinch): angle of rotation about the y-axis
   g(gamma of roll): angle of rotation about the x-axis
   cell_mat: the 3D point matrix"
  mut_mat <- t(rotation(a, b, g) %*% t(cell_mat))
  colnames(mut_mat) <- c("X", "Y", "Z")
  # Convert mutation_matrix to dataframe object to apply slice function
  mut_frame <- as.data.frame(mut_mat)
  # Create a slice of the 3D point pattern
  slices <- mut_frame %>%
    filter(slice(0, 0, 1, 10, as.matrix(mut_frame[, 1:3])))
  print(sprintf("Number of cells in slice: %s", dim(slices)[1]))
  # Create a 2D point pattern object
  pp_obj <- spatstat.geom::ppp(slices$X, slices$Y,
                               pp_box(as.matrix(slices[, 1:2])))
  # Compute the 2D K-Ripley function
  k <- spatstat.explore::Kest(pp_obj)
  pcf <- spatstat.explore::pcf.fv(k, method = "b", spar = 0.5)

  list(pcf = pcf, num_cells = dim(slices)[1]) # return pcf and number of cells
}