
# INSTALLING PACKAGES
setwd(this.path::here())
source("module.r")
source("visual.r")

library(spatstat) # for spatial statistics

# PARSING ARGUMENTS
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Requires 1 argument: path to the file")
}

# VARIABLES
# Downloading the spatial centered data of specific cell type
coords <- read.csv(args[1]) # parsing the path argument
cell_type <- coords[1, "cell"] # reading the cell type
data_name <- basename(dirname(args[1])) # technology name
ce_tbl <- tibble(data = character(), cell = character(),
  dim = character(), num = integer(), ce = double()
) # empty table for storing the ce values
cell_mat <- as.matrix(coords[, 1:3]) # 3D point pattern matrix

# COMPUTING TRUE SPATIAL PCF
true_pattern <- spatstat.geom::pp3(cell_mat[, 1], cell_mat[, 2],
                                   cell_mat[, 3], pp_box(cell_mat))
pcf_true <- spatstat.explore::pcf3est(true_pattern)

# COMPUTING PCF FOR SLICES
# Define a vector of angles
n_iter <- 10
rolls <- runif(n = n_iter, min = 0, max = 2 * pi)
pinches <- runif(n = n_iter, min = 0, max = 2 * pi)
# yaws <- runif(n = n_iter, min = 0, max = 2 * pi) just in case

# Create an empty lists to store the pcf objects and number of cells
pcf_list <- vector("list", length = n_iter^2)
num_cells_list <- numeric(length = n_iter^2)
valid_pairs <- matrix(nrow = 0, ncol = 2)  # column 1 = roll, column 2 = pinch
n_cell_range <- c(-1, 100000) # set according to the pcf plot

# Loop over the roll angles
index <- 1
for (roll in rolls) {
  for (pinch in pinches) {
    # Call the pc_for_slice function to compute pcf for slice of 3D dataframe
    result <- pc_for_slice(a = 0, b = pinch, g = roll, cell_mat)
    pcf_list[[index]] <- result$pcf
    num_cells_list[index] <- result$num_cells

    if (!is.na(result$num_cells)  && result$num_cells > 1) {
      ce_tbl <- bind_rows(ce_tbl,
        tibble(data = data_name, cell = cell_type,
          dim = "plane", num = result$num_cells, ce = result$ce
        )
      )

      if ((n_cell_range[1] <= result$num_cells) &&
            (result$num_cells <= n_cell_range[2])) {
        valid_pairs <- rbind(valid_pairs, c(roll, pinch))
      }
    } # activated when n_cell_range is set
    index <- index + 1
  }
}

# Add space data to the ce table
ce_tbl <- bind_rows(ce_tbl,
  tibble(data = data_name, cell = cell_type, dim = "space",
    num = nrow(coords), ce = ce_idx(coords[, c("x", "y", "z")])
  )
)

# PLOTTING THE PCF CURVES
cumul_pcf(pcf_true, pcf_list, cell_type, num_cells_list)

# SAVING THE FILES

dir.create("../logs", showWarnings = FALSE)
writeLines(c(paste("Executed Script: angle.r"), capture.output(sessionInfo())),
           file.path("../logs", paste0("logs_angle_", Sys.Date(), ".txt")))
utils::write.csv(ce_tbl, sprintf("../data/ce_temp/ce_tbl_%s.csv",
                                 cell_type), row.names = FALSE)

# SLICE ANALYSIS
if (n_cell_range[1] != -1) angle_analysis(cell_mat, valid_pairs, cell_type)

# Print optimal r max for slice.r script
max_index <- which.max(num_cells_list)
print(sprintf("Optimal max_r: %s", max(pcf_list[[max_index]]$r)))