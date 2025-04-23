
# INSTALLING LIBRARIES & ARGUMENTS PARSING
library(spatstat) # spatial statistics
library(tibble) # dataframes

setwd(this.path::here())
source("./module.r")
source("./visual.r")

# PARSING ARGUMENTS
args <- commandArgs(TRUE) # parsing the command line arguments
if (length(args) < 2) {
  stop("Requires 2 arguments: r_max and path to the file")
}
r_grid <- seq(0, as.numeric(args[2]), length.out = 513)

# VARIABLES' DECLARATION
# Downloading the spatial centered data of specific cell type
coords <- read.csv(args[1]) # using the path argument from command line
cell_type <- coords[1, "cell"] # reading the cell type
dataset_name <- basename(dirname(args[1])) # reading the dataset name
cell_mat <- as.matrix(coords[, 1:3])

# COMPUTING TRUE SPATIAL PCF
true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))
pcf_true <- spatstat.explore::pcf3est(true_pattern)

# COMPUTING PCF FOR SLICES
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

# VISUALIZATION
mean_pcf(pcf_true, r_grid, pcf_mat, cell_type, num_cells_list) # module visual.r

# ARCHIVISATION
dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: slice.r"), capture.output(sessionInfo())),
           file.path("../logs", paste0("logs_slice_", Sys.Date(), ".txt")))

# Final edits
pcf_plain <- data.frame(r = r_grid, pcf_mat, check.names = FALSE)
colnames(pcf_plain)[-1] <- paste0("pcf_", seq_len(ncol(pcf_mat)))

pcf_space <- tibble(r = pcf_true$r, pcf = pcf_true$iso)

utils::write.csv(pcf_plain, sprintf("../data/pcf_temp/pcf_plain_%s.csv",
                                    cell_type), row.names = FALSE)
utils::write.csv(pcf_space, sprintf("../data/pcf_temp/pcf_space_%s.csv",
                                    cell_type), row.names = FALSE)