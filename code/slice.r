
source("./module.r")
getwd() # wd need to be set explicitly, in terminal :(
library(spatstat)
library(viridis)    # for color palettes
library(scales)     # for col_numeric
library(fields)     # for image.plot (colorbar)
library(RColorBrewer)
library(Hmisc)

# Downloading the spatial centered data of specific cell type
coords <- read.csv("../data/selected_cell_coordinates.csv")

# COMPUTING PCF FOR SLICES
# Creating a point pattern object
cell_mat <- as.matrix(coords[, 1:3])

true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))

# Define a vector of angles
rolls <- runif(n = 30, min = 0, max = 2 * pi)
pinches <- runif(n = 30, min = 0, max = 2 * pi)
# yaws <- runif(n = 30, min = 0, max = 2 * pi) have no effect

# Create an empty lists to store the pcf objects and number of cells
pcf_list <- vector("list", length = length(rolls) * length(pinches))
num_cells_list <- numeric(length = length(rolls) * length(pinches))

# Loop over the roll angles
index <- 1
for (roll in rolls){
  for (pinch in pinches){
    # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
    result <- pc_for_slice(a = 0, b = pinch, g = roll,
                           cell_mat = cell_mat)
    pcf_list[[index]] <- result$pcf
    num_cells_list[index] <- result$num_cells
    index <- index + 1
  }
}

# COMPUTING TRUE SPATIAL PCF
pcf_true <- spatstat.explore::pcf3est(true_pattern)

# COMPUTING MEAN AND SD OF ALL 2D PCF
# Define a shared r grid (using the r values from true PCF)
r_grid <- pcf_true$r

# Function to interpolate a PCF object onto the r grid
interp_pcf <- function(pcf_list, r_grid) {
  approx(x = pcf_list$r, y = pcf_list$pcf, xout = r_grid, rule = 2)$y
}

# Interpolate each PCF in your list onto r_grid
pcf_matrix <- sapply(pcf_list, interp_pcf, r_grid = r_grid)

# Compute the mean and SD across all curves at each common r value
weighted_mean_pcf <- apply(pcf_matrix, 1, function(x) {
  weighted.mean(x, weights = num_cells_list, na.rm = TRUE)
})
weighted_sd_pcf <- apply(pcf_matrix, 1, function(x) {
  sqrt(wtd.var(x, weights = num_cells_list, na.rm = TRUE))
})

# VISUALIZATION
# Plot the true PCF, then overlay the mean PCF with ±1 standard deviation
png(filename = "../images/pc_mean_sd.png", width = 800, height = 600)

y_max <- max(max(pcf_true$iso), max(weighted_mean_pcf + weighted_sd_pcf))
plot(pcf_true, col = "#84a98c", lwd = 2, ylim = c(0, y_max),
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
writeLines(c(paste("Executed Script: slice.r"),
             capture.output(sessionInfo())),
           file.path("../logs",
                     paste0("logs_slice_", Sys.Date(), ".txt")))