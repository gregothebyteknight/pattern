
# SETTING THE ENVIRONMENT
setwd(this.path::here())
source("./module.r")
source("./visual.r")

suppressPackageStartupMessages({
  library(spatstat.geom) # geometric objects (ppp, pp3, box, owin)
  library(spatstat.explore) # exploratory analysis (Kest, pcf3est)
  library(tibble) # dataframes
})

# PARSING ARGUMENTS
args <- commandArgs(TRUE) # parsing the command line arguments
if (length(args) < 2) {
  stop("Requires 2 arguments: r_max and path to the file")
}

# VARIABLES' DECLARATION
# Downloading the spatial centered data of specific cell type
coords <- read.csv(args[1]) # using the path argument from command line
r_grid <- seq(0, as.numeric(args[2]), length.out = 513)
cell_type <- coords[1, "cell"] # reading the cell type
cell_mat <- data.matrix(coords[, c("x", "y", "z")])

# COMPUTING TRUE SPATIAL PCF
true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))
pcf_true <- spatstat.explore::pcf3est(true_pattern, correction = "isotropic")

# COMPUTING PCF FOR SLICES
# Define a vector of angles
n_iter <- 900
rolls <- runif(n = n_iter, min = 0, max = 2 * pi)
pinches <- runif(n = n_iter, min = 0, max = 2 * pi)
# yaws <- runif(n = n_iter, min = 0, max = 2 * pi) have no effect

# Create an empty lists to store the pcf objects and number of cells
pcf_mat <- matrix(ncol = n_iter, nrow = 513) # locations of dots
num_cells_list <- numeric(n_iter)

# Loop over the roll angles
for (i in 1:n_iter) {
  # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
  result <- pc_for_slice(a = 0, b = pinches[i], g = rolls[i], cell_mat, r_grid)
  g <- result$pcf$pcf
  length(g) <- length(r_grid) # expanding the vector to 513
  pcf_mat[, i] <- g # 513 pcf values
  num_cells_list[i] <- result$num_cells
}

# VISUALIZATION
mean_pcf(pcf_true, r_grid, pcf_mat, cell_type, num_cells_list) # module visual.r

# ARCHIVISATION
dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: slice.r"), capture.output(sessionInfo())),
           file.path("../logs", sprintf("logs_slice_%s.txt", Sys.Date())))

# Final edits
pcf_plane <- data.frame(r = r_grid, pcf_mat, check.names = FALSE)
colnames(pcf_plane)[-1] <- paste0("pcf_", seq_len(ncol(pcf_mat)))

pcf_space <- tibble(r = pcf_true$r, pcf = pcf_true$iso)

pcf_num <- tibble(dim = "plane", num = num_cells_list)
pcf_num <- bind_rows(pcf_num, tibble(dim = "space", num = nrow(cell_mat)))

path <- sprintf("../data/pcf/%s/%s", basename(dirname(args[1])), cell_type)
dir.create(path, showWarnings = FALSE, recursive = TRUE)
# write.csv(pcf_plane, file.path(path, "pcf_plane.csv"), row.names = FALSE)
# write.csv(pcf_space, file.path(path, "pcf_space.csv"), row.names = FALSE)
# write.csv(pcf_num, file.path(path, "pcf_num.csv"), row.names = FALSE)