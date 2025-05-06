# Installing the libraries
suppressPackageStartupMessages({
  library(tidyverse) # data manipulation
  library(spatstat) # spatial statistics
  library(dplyr) # data manipulation with pipes
  library(FNN) # fast nearest neighbor search
  library(geometry) # volume computation
})
# Define the ranges for each coordinate
pp_box <- function(mat) {
  "Define the range of the point pattern object
  @mat(matrix): the 2D or 3D point pattern matrix
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
   @a(alpha of yaw): angle of rotation about the z-axis
   @b(beta of pinch): angle of rotation about the y-axis
   @g(gamma of roll): angle of rotation about the x-axis"
  ca <- cos(a) # precompute cos and sin values
  sa <- sin(a)
  cb <- cos(b)
  sb <- sin(b)
  cg <- cos(g)
  sg <- sin(g)

  matrix(c(
    ca * cb, ca * sb * sg - sa * cg, ca * sb * cg + sa * sg,
    sa * cb, sa * sb * sg + ca * cg, sa * sb * cg - ca * sg,
    -sb, cb * sg, cb * cg
  ), nrow = 3, byrow = TRUE)
}

# Declare the slice function
slice <- function(a, b, c, slice_size, mut_frame) {
  "Compute the slice of the 3D point pattern.
   @a, b, c: the normal vector of the slice (d = 0)
   @slice_size: the thickness of the slice
   @mut_frame: the 3D point dataframe which undergo transformation"
  value <- a * mut_frame[, 1] + b * mut_frame[, 2] + c * mut_frame[, 3]
  # Return a logical vector: TRUE if value is within the slice, FALSE otherwise
  value >= -slice_size / 2 & value <= slice_size / 2
}

# Define the function to compute the dissections
pc_for_slice <- function(a, b, g, cell_mat, r_grid = NULL) {
  "Compute the dissections of the 3D point pattern
   @a(alpha of yaw): angle of rotation about the z-axis
   @b(beta of pinch): angle of rotation about the y-axis
   @g(gamma of roll): angle of rotation about the x-axis
   @cell_mat: the 3D point matrix"
  mut_mat <- t(rotation(a, b, g) %*% t(cell_mat))
  colnames(mut_mat) <- c("X", "Y", "Z")
  # Convert mutation_matrix to dataframe object to apply slice function
  mut_frame <- as.data.frame(mut_mat)
  # Create a slice of the 3D point pattern
  slices <- mut_frame %>%
    filter(slice(0, 0, 1, 10, as.matrix(mut_frame[, 1:3])))
  print(sprintf("Number of cells in slice: %s", dim(slices)[1]))
  if (dim(slices)[1] <= 1) {
    return(list(pcf = data.frame(r = runif(n = 513, min = 0, max = 10),
                                 pcf = rep(NA, 513)), num_cells = NA))
  }
  # Create a 2D point pattern object
  pp_obj <- spatstat.geom::ppp(slices$X, slices$Y,
                               pp_box(as.matrix(slices[, 1:2])))
  # Compute the 2D PCF function
  k <- spatstat.explore::Kest(pp_obj, r = r_grid, correction = "isotropic")
  pcf <- spatstat.explore::pcf.fv(k, method = "b", spar = 0.5)

  # Compute 2D Clark and Evans index and return the results
  list(pcf = pcf, num_cells = dim(slices)[1], ce = ce_idx(slices))
}

angle_analysis <- function(cell_mat, valid_pairs, cell_type) {
  "Image slice of point pattern corresponds to the 
  selected n_cell_range
  @cell_mat: the 3D pint matrix"
  valid_pairs <- as.data.frame(valid_pairs)
  colnames(valid_pairs) <- c("roll", "pinch")
  valid_pairs$angle_sum <- valid_pairs$roll + valid_pairs$pinch

  # Identify the row with the minimal sum
  mean_sum <- mean(valid_pairs$angle_sum)
  closest_index <- which.min(abs(valid_pairs$angle_sum - mean_sum))
  closest_pair <- valid_pairs[closest_index, c("roll", "pinch")]

  pinch_value <- as.numeric(closest_pair["pinch"])
  roll_value <- as.numeric(closest_pair["roll"])

  mut_mat <- t(rotation(0, pinch_value,
                        roll_value) %*% t(cell_mat))
  colnames(mut_mat) <- c("X", "Y", "Z")
  # Convert mutation_matrix to dataframe object to apply slice function
  mut_frame <- as.data.frame(mut_mat)
  # Create a slice of the 3D point pattern
  slices <- mut_frame %>%
    filter(slice(0, 0, 1, 10, as.matrix(mut_frame[, 1:3])))
  print(sprintf("Number of cells in selected slice: %s", dim(slices)[1]))

  png(filename = sprintf("../images/angle_%s.png", cell_type))
  plot(slices[, "X"], slices[, "Y"], lwd = 2,
       xlab = "X axis", ylab = "Y axis",
       main = "Slice of cell with the closest angle sum")
}

ce_idx <- function(df) {
  "Compute the Clark and Evans index for a set of points in 2D or 3D space
  @df(data.frame): A data frame containing the coordinates of the points"
  # R EXPECTED
  if (nrow(df) < 4) return(NA) # not enough points to compute the index
  if (ncol(df) == 2) {
    pts <- spatstat.geom::ppp(df[, 1], df[, 2],
                              window = spatstat.geom::owin(range(df[, 1]),
                                                           range(df[, 2])))
    ce_2d <- spatstat.explore::clarkevans(pts, correction = "Donnelly")
    return(as.numeric(ce_2d))

  } else if (ncol(df) == 3) {
    # Compute the convex geometry of the points
    ch <- geometry::convhulln(df, options = "FA")
    surface <- ch$area
    volume <- ch$vol

    r_exp <- gamma(4 / 3) * (4 / 3 * pi * (nrow(df) / volume))^(-1 / 3) *
      (1 + 0.75 * (surface / (nrow(df)^(2 / 3) * volume^(2 / 3))))
  } else {
    stop("Data frame must have 2 or 3 columns.")
  }

  # R OBSERVED
  nn <- FNN::get.knnx(df, df, k = 2) # find the nearest neighbor
  r_obs <- mean(nn$nn.dist[, 2]) # mean distance to the nearest neighbor

  # CLARK AND EVANS INDEX
  round(r_obs / r_exp, digits = 3)
}