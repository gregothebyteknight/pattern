# LOADING LIBRARIES
suppressPackageStartupMessages({
  library(pheatmap)
})

# Downloading Data (expression matrix, table metal - marker)
expr <- read.delim("../data/cell_intensities.csv", sep = ",")
panel <- read.delim("../data/model_201709_panel.csv", sep = ",")

chan_full <- colnames(expr)
colnames(expr) <- panel$clean_target[match(chan_full, panel$Metal.Tag)]
# Plot Heatmap for uncompensated expression matrix
pheatmap::pheatmap(cor(expr), clustering_method = "ward",
                   filename = "../images/heat_unc.png")

# FILTERING COMPENSATION MATRIX
comp_mat <- read.delim("../data/compensationMatrix.csv", sep = ",",
                       row.names = 1)
chan_true <- chan_full[chan_full %in% colnames(comp_mat)] # in comp matrix
print(sprintf("Dim of Compensation Matrix before filter %s", dim(comp_mat)))
chan_to_add <- chan_full[!chan_full %in% colnames(comp_mat)] # !in comp matrix
comp_mat <- comp_mat[chan_true, chan_true] # select part of mat with chan_true

# EXPAND COMPENSATION MATRIX WITH ZEROS
comp_mat <- cbind(comp_mat, matrix(0, ncol = length(chan_to_add),
                                   nrow = nrow(comp_mat)))
comp_mat <- rbind(as.matrix(comp_mat), matrix(0, nrow = length(chan_to_add),
                                              ncol = ncol(comp_mat)))
rownames(comp_mat) <- c(chan_true, chan_to_add)
colnames(comp_mat) <- c(chan_true, chan_to_add)
comp_mat <- comp_mat[chan_full, chan_full] # full-scale comp matrix
diag(comp_mat) <- 1
print(sprintf("Dim of Compensation Matrix after filter %s", dim(comp_mat)))

# APPLY COMPENSATION
comp_mat_inv <- solve(comp_mat) # find inverse of comp matrix
expr_corr <- as.matrix(expr) %*% comp_mat_inv

colnames(expr_corr) <- panel$clean_target[match(chan_full, panel$Metal.Tag)]
pheatmap::pheatmap(cor(expr_corr), filename = "../images/heat_cor.png")

# ADDITIONAL CLEANING (OPTIONAL)
col_to_remove <- c("Ir193", "Ir191", "Iridium", "Ruthenium", "NaN")
expr_corr_filt <- expr_corr[, !colnames(expr_corr) %in% col_to_remove]
print(sprintf("N_col of Filtered Corrected Mat %s", dim(expr_corr_filt)))

# SAVING EXPRESSION FILES
write.table(expr_corr_filt, "../data/expression_annotated_corrected.csv",
            sep = ",", row.names = FALSE)