
# Load libraries
setwd(this.path::here())
library(pheatmap)

# P value heatmap
suppressPackageStartupMessages({
  path <- "../temp/"
})

pval_all <- read.csv("../temp/pval_tbl.csv")
pval_all <- pval_all |>
  dplyr::select(data, cell)

reveal_all <- function(pval_all) {
  # DATA AGGREGATION
  for (type in list.dirs(path = "../temp/", recursive = FALSE)) {
    print(sprintf("Statistics type: %s", basename(type)))
    pval_tbl <- read.csv(file.path(type, "pval_tbl.csv"))
    pval_tbl$p[pval_tbl$p == 0] <- 1e-3
    type_name <- basename(type) # now need only the dataset name, not path
    print(-log(pval_tbl$p, base = 10))
    pval_all[[type_name]] <- -log(pval_tbl$p, base = 10)
  }
  pval_all # return the data frame with all the data)
}

pval_all <- reveal_all(pval_all)
write.csv(pval_all, file = "../temp/pval_all.csv", row.names = FALSE)

pval_all <- read.csv("../temp/pval_all.csv", stringsAsFactors = FALSE)

# 2) Build your matrix: rows = unique rows, cols = tests
mat <- as.matrix(pval_all[, -(1:2)])    # drop data + cell cols
rownames(mat) <- paste(pval_all$data, pval_all$cell, sep = "|")

# 3) Build annotation_row correctly against mat’s rownames
annotation_row <- data.frame(Dataset = factor(pval_all$data),
                             row.names = rownames(mat))

# 4) Plot
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_row,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "-log10(p) across cell × test",
  filename = "../images/FINALS/heatmap/pval_heat.png"
)
