
# INSTALLING PACKAGES
setwd(this.path::here())
source("module.r")

library(spatstat) # for spatial statistics
library(scales) # for col_numeric
library(fields) # for image.plot (colorbar)
library(RColorBrewer) # for color palettes

# VARIABLES
# Downloading the spatial centered data of specific cell type
coords <- read.csv("../data/sci_embryo/selected_cell_coordinates.csv")
thr <- 0.01 / log10(nrow(coords)) # threshold for the PCF difference

# COMPUTING PCF FOR SLICES
# Creating a point pattern object
cell_mat <- as.matrix(coords[, 1:3])

true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))

# Define a vector of angles
n_iter <- 10
rolls <- runif(n = n_iter, min = 0, max = 2 * pi)
pinches <- runif(n = n_iter, min = 0, max = 2 * pi)
# yaws <- runif(n = n_iter, min = 0, max = 2 * pi) just in case

# Create an empty lists to store the pcf objects and number of cells
pcf_list <- vector("list", length = n_iter^2)
num_cells_list <- numeric(length = n_iter^2)
valid_pairs <- matrix(nrow = 0, ncol = 2)  # column 1 = roll, column 2 = pinch
n_cell_range <- c(0, 100000) # set according to the pcf plot

# Loop over the roll angles
index <- 1
for (roll in rolls) {
  for (pinch in pinches) {
    # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
    result <- pc_for_slice(a = 0, b = pinch, g = roll,
                           cell_mat = cell_mat)
    pcf_list[[index]] <- result$pcf
    num_cells_list[index] <- result$num_cells
    if (n_cell_range[1] != -1) {
      if ((n_cell_range[1] <= result$num_cells) &&
            (result$num_cells <= n_cell_range[2])) {
        valid_pairs <- rbind(valid_pairs, c(roll, pinch))
      }
    } # activated when n_cell_range is set
    index <- index + 1
  }
}

# COMPUTING TRUE SPATIAL PCF
pcf_true <- spatstat.explore::pcf3est(true_pattern)

# PLOTTING THE PCF CURVES
# Plot the true PCF
png(filename = "../images/pc_total.png")
plot(pcf_true, main = "Pair Correlation Functions for Varying Angles",
     xlab = "r (Distance)", ylab = "pcf", col = "#84a98c")

# Prepare continuous color panel
col_fun <- scales::col_numeric(palette = viridis::cividis(100),
                               domain = range(num_cells_list))
curve_colors <- col_fun(num_cells_list)

# Overlay the rest of the pcf curves
for (i in seq_along(pcf_list)) {
  points(pcf_list[[i]]$r, pcf_list[[i]]$pcf, type = "l", col = curve_colors[i])
}

# Adjust zlim if num_cells_list has no variation
z_range <- range(num_cells_list)
if (z_range[1] == z_range[2]) {
  z_range <- z_range + c(-1, 1)  # add a small margin
}

# Saving the plot and session info
fields::image.plot(legend.only = TRUE,
                   zlim = z_range,
                   col = viridis::cividis(100),
                   legend.lab = "Number of cells",
                   legend.line = 2,
                   legend.mar = 3)
dev.off()

dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: main.r"), capture.output(sessionInfo())),
           file.path("../logs", paste0("logs_main_", Sys.Date(), ".txt")))

# Slice analysis
angle_analysis(cell_mat, valid_pairs) # only when n_cell_range is set

# Print optimal r max for slice.r script
max_index <- which.max(num_cells_list)
print(sprintf("Optimal max_r: %s", max(pcf_list[[max_index]]$r)))