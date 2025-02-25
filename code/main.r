
source("./slices.r")
setwd("./code")
getwd()
library(spatstat)
library(viridis)    # for color palettes
library(scales)     # for col_numeric
library(fields)     # for image.plot (colorbar)
library(RColorBrewer)

# Downloading the spatial centered data of specific cell type
coords <- read.csv("../data/selected_cell_coordinates.csv",
                   row.names = 1) # data for cancer cells

# Creating a point pattern object
cell_mat <- as.matrix(coords[, 1:3])

true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))

# Define a vector of angles
angles <- seq(0, 2 * pi, length.out = 10)

# Create an empty lists to store the pcf objects and number of cells
pcf_list <- vector("list", length(angles * angles))
num_cells_list <- numeric(length(angles * angles))

# Loop over the roll angles
index <- 1
for (roll in seq_along(angles)){
  for (pinch in seq_along(angles)){
    # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
    result <- pc_for_slice(a = 0, b = angles[pinch], g = angles[roll],
                           cell_mat = cell_mat)
    pcf_list[[index]] <- result$pcf
    num_cells_list[index] <- result$num_cells
    max_pcf <- max(max_pcf, max(pcf_list[[i]]$pcf)) # for plot's axis setting
    index <- index + 1
  }
}

# Plot the first pcf curve
pcf_true <- spatstat.explore::pcf3est(true_pattern)
max_pcf <- max(pcf_true$pcf)

png(filename = "../images/pc_plot.png")
plot(pcf_true,
     main = "Pair Correlation Functions for Varying Angles",
     xlab = "r (Distance)", ylab = "pcf", col = "#84a98c")

# Prepare continuous color panel
col_fun <- scales::col_numeric(palette = viridis::cividis(100),
                               domain = range(num_cells_list))
curve_colors <- col_fun(num_cells_list)

# Overlay the rest of the pcf curves
for (i in seq_along(pcf_list)) {
  points(pcf_list[[i]]$r, pcf_list[[i]]$pcf, type = "l", col = curve_colors[i])
}

# Saving the plot and session info
fields::image.plot(legend.only = TRUE,
                   zlim = range(num_cells_list),
                   col = viridis::cividis(100),
                   legend.lab = "Number of cells",
                   legend.line = 2,
                   legend.mar = 3)
dev.off()

dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: main.r"),
             capture.output(sessionInfo())),
           file.path("../logs",
                     paste0("logs_main_", Sys.Date(), ".txt")))
