
# INSTALLING LIBRARIES & ARGUMENTS PARSING
library(Hmisc) # for weighted variance
library(spatstat) # spatial statistics

setwd(this.path::here())
source("./module.r")

# PARSING ARGUMENTS
args <- commandArgs(TRUE) # parsing the command line arguments
if (length(args) < 2) {
  stop("Requires 2 arguments: r_max and path to the file")
}

r_grid <- seq(0, as.numeric(args[2]), length.out = 513)

# VARIABLES' DECLARATION
# Downloading the spatial centered data of specific cell type
coords <- read.csv(args[1]) # using the path argument from command line

# COMPUTING PCF FOR SLICES
# Creating a point pattern object
cell_mat <- as.matrix(coords[, 1:3])
cell_type <- coords[1, "cell"] # reading the cell type

true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))

# Define a vector of angles
n_iter <- 30
rolls <- runif(n = n_iter, min = 0, max = 2 * pi)
pinches <- runif(n = n_iter, min = 0, max = 2 * pi)
# yaws <- runif(n = n_iter, min = 0, max = 2 * pi) have no effect

# Create an empty lists to store the pcf objects and number of cells
pcf_mat <- matrix(ncol = n_iter^2, nrow = 513) # locations of dots
num_cells_list <- numeric(length = length(rolls) * length(pinches))

# Loop over the roll angles
index <- 1
for (roll in rolls) {
  for (pinch in pinches) {
    # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
    result <- pc_for_slice(a = 0, b = pinch, g = roll, cell_mat, r_grid)
    pcf_mat[, index] <- result$pcf$pcf # 513 pcf values
    num_cells_list[index] <- result$num_cells
    index <- index + 1
  }
}

# COMPUTING TRUE SPATIAL PCF
pcf_true <- spatstat.explore::pcf3est(true_pattern)

# COMPUTING MEAN AND SD OF ALL 2D PCF
# Compute the mean and SD across all curves at each common r value
weighted_mean_pcf <- apply(pcf_mat, 1, function(x) {
  weighted.mean(x, weights = num_cells_list, na.rm = TRUE)
})
weighted_sd_pcf <- apply(pcf_mat, 1, function(x) {
  sqrt(Hmisc::wtd.var(x, weights = num_cells_list, na.rm = TRUE))
})

# VISUALIZATION
# Plot the true PCF, then overlay the mean PCF with ±1 standard deviation
png(filename = sprintf("../images/pc_mean_%s.png", cell_type),
    width = 800, height = 600)

y_max <- max(max(pcf_true$iso), max(weighted_mean_pcf + weighted_sd_pcf))
plot(pcf_true, col = "#84a98c", lwd = 2,
     xlab = "r (Distance)", ylab = "pcf", main = "Mean PCF with ±1 SD")

# Addition of dashed lines for mean ± SD
lines(r_grid, weighted_mean_pcf + weighted_sd_pcf, col = "gray", lty = 2)
lines(r_grid, weighted_mean_pcf - weighted_sd_pcf, col = "gray", lty = 2)

# Addition of a shaded region between mean + SD and mean - SD
polygon(c(r_grid, rev(r_grid)),
        c(weighted_mean_pcf + weighted_sd_pcf,
          rev(weighted_mean_pcf - weighted_sd_pcf)),
        col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)

# Additiom of the true PCF for reference
lines(r_grid, weighted_mean_pcf, type = "l", lwd = 2, col = "#7373c9")

dev.off()

dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: slice.r"), capture.output(sessionInfo())),
           file.path("../logs", paste0("logs_slice_", Sys.Date(), ".txt")))