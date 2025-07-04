
# SETTING THE ENVIRONMENT
setwd(this.path::here())
source("module.r")
source("visual.r")

suppressPackageStartupMessages({
  library(spatstat.geom) # geometric objects (ppp, pp3, box, owin)
  library(spatstat.explore) # exploratory analysis (Kest, pcf3est)
  library(tibble) # dataframes
  library(dplyr) # data manipulation
})

# PARSING ARGUMENTS
args <- commandArgs(TRUE)
if (length(args) < 1) stop("Requires 1 argument: path to the file")

# VARIABLES
# Downloading the spatial centered data of specific cell type
coords <- read.csv(args[1]) # parsing the path argument
cell_type <- coords[1, "cell"] # reading the cell type
cell_mat <- data.matrix(coords[, c("x", "y", "z")]) # 3D point pattern matrix

ce_tbl <- tibble(dim = character(), num = integer(), ce = double())
pcf_naive <- tibble(dim = character(), pcf = double(),
                    r = double(), num = integer())

# COMPUTING TRUE SPATIAL PCF
true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))
pcf_true <- spatstat.explore::pcf3est(true_pattern, correction = "isotropic")

# COMPUTING PCF FOR SLICES
# Define a vector of angles
n_iter <- 900 # number of iterations
rolls <- runif(n = n_iter, min = 0, max = 2 * pi)
pinches <- runif(n = n_iter, min = 0, max = 2 * pi)
# yaws <- runif(n = n_iter, min = 0, max = 2 * pi) just in case

# Create an empty lists to store the pcf objects and number of cells
pcf_list <- vector("list", length = n_iter)
num_cells_list <- numeric(length = n_iter)
valid_pairs <- matrix(nrow = 0, ncol = 2)  # column 1 = roll, column 2 = pinch
if (length(args) == 3) n_range <- as.numeric(args[2:3]) else n_range <- c(-1, 1)
print(sprintf("Selected range of number of cells: %s", n_range))

# Loop over the roll angles
for (i in 1:n_iter) {
  # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
  res <- pc_for_slice(a = 0, b = pinches[i], g = rolls[i], cell_mat = cell_mat)
  pcf_list[[i]] <- res$pcf
  num_cells_list[i] <- res$num_cells

  if (!is.na(res$num_cells)  && res$num_cells > 1) {
    ce_tbl <- bind_rows(ce_tbl, tibble(dim = "plane",
                                       num = res$num_cells, ce = res$ce))
    pcf_naive <- bind_rows(pcf_naive, tibble(dim = "plane",
                                             pcf = max(res$pcf$pcf,
                                                       na.rm = TRUE),
                                             r = max(res$pcf$r, na.rm = TRUE),
                                             num = res$num_cells))
    if ((n_range[1] <= res$num_cells) &&
          (res$num_cells <= n_range[2])) {
      valid_pairs <- rbind(valid_pairs, c(rolls[i], pinches[i]))
    }
  } # activated when n_range is set
}

# Add space data to the ce table
ce_tbl <- bind_rows(ce_tbl, tibble(dim = "space", num = nrow(coords),
                                   ce = ce_idx(coords[, c("x", "y", "z")]))
)
pcf_naive <- bind_rows(pcf_naive, tibble(dim = "space", pcf = max(pcf_true$iso,
                                                                  na.rm = TRUE),
                                         r = max(pcf_true$r, na.rm = TRUE),
                                         num = nrow(coords)))

# PLOTTING THE PCF CURVES
cumul_pcf(pcf_true, pcf_list, cell_type, num_cells_list, n_range)

# SAVING THE FILES
dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: angle.r"), capture.output(sessionInfo())),
           file.path("../logs", sprintf("logs_angle_%s.txt", Sys.Date())))
path <- sprintf("../data/naive/%s/%s", basename(dirname(args[1])), cell_type)
dir.create(path, showWarnings = FALSE, recursive = TRUE)
# write.csv(ce_tbl, file.path(path, "ce_tbl.csv"), row.names = FALSE)
# write.csv(pcf_naive, file.path(path, "pcf_naive.csv"), row.names = FALSE)

# SLICE ANALYSIS
if (n_range[1] != -1) angle_analysis(cell_mat, valid_pairs, cell_type)

# Print optimal r max for slice.r script
max_index <- which.max(num_cells_list)
print(sprintf("Optimal max_r: %s", round(max(pcf_list[[max_index]]$r), 2)))