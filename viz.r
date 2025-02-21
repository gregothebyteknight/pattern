
source("sections.r")
library(spatstat)

pp_box_3 <- function(mat) {
  x_range <- c(min(mat[, 1]), max(mat[, 1]))
  y_range <- c(min(mat[, 2]), max(mat[, 2]))
  z_range <- c(min(mat[, 3]), max(mat[, 3]))

  # Create the 3D box (an object of class "box3")
  spatstat.geom::box3(xrange = x_range, yrange = y_range, zrange = z_range)
}

# Downloading the spatial centered data of specific cell type
coords <- read.csv("cell_type_coordinates.csv",
                   row.names = 1) # data for cancer cells

# Convertion the coordinates from the data to matrix
cell_mat <- as.matrix(coords[, 1:3])

true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box_3(cell_mat))

# Define a vector of angles
angles <- seq(0, 2 * pi, length.out = 10)
print(roll_angles)
print(cos(angles))

# Create an empty list to store the pcf objects
pcf_list <- vector("list", length(angles))

# Loop over the roll angles
for (i in seq_along(angles)) {
  angle <- angles[i]
  # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
  pcf_list[[i]] <- pc_for_slice(a = 0, b = angle, g = 0,
                                cell_mat = cell_mat)
}

# Plot the first pcf curve
pcf_true <- spatstat.explore::pcf3est(true_pattern)


plot(pcf_true,
     main = "Pair Correlation Functions for Varying Yaw Angles",
     xlab = "r (Distance)", ylab = "pcf", col = "red")

# Overlay the rest of the pcf curves
for (i in seq_along(pcf_list)) {
  points(pcf_list[[i]]$r, pcf_list[[i]]$pcf, type = "l")
}

png(filename = "./images/pc_3d_to_2d.png")