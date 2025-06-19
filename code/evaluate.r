
# SETTING UP THE ENVIRONMENT
setwd(this.path::here())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

if (!capabilities("cairo")) {
  stop("Your R was compiled without Cairo support.")
}

ce_table <- read.csv("../temp/pval_tbl.csv")

ce_table <- ce_table |>
  mutate(p_log = -log(p.BH, base = 10), cell = factor(cell)) |>
  select(data, cell, p_log, n)

deep_cols <- c(
  "nat_blood" = "#4C72B0",
  "nat_mainher2" = "#55A868",
  "nat_lymph" = "#C44E52",
  "nat_secher2" = "#8172B2",
  "sci_brain" = "#CCB974",
  "sci_embryo" = "#64B5CD"
)

pretty_labels <- c(
  "LVIBlood...",
  "MainHer2...",
  "LVILymph...",
  "SecondHer2...",
  "Hippocampal...",
  "Embryo..."
)


ce_plot <- ggplot(ce_table, aes(y = cell, x = p_log, fill = data)) +
  geom_col(width = 1, position = position_dodge2(width = 1,
                                                 preserve = "single",
                                                 padding = 0)) +
  # geom_text(aes(label = paste0("Number of slices = ", n)),
  #           position = position_dodge2(width = 0.9),
  #           hjust = -0.1, size = 3, check_overlap = TRUE) +
  scale_y_discrete(expand = expansion(add = 0)) +
  geom_vline(aes(xintercept = 1.3), linetype = "dashed", color = "grey") +
  labs(title = "Comparison of p-values (FDR-BH) between 2D and 3D values",
       y = "Cell Type", x = "-log10(adj p-value)", fill  = "Dataset name") +
  scale_fill_manual(
    name = "Dataset name", # legend title
    values = deep_cols, # your deep palette
    labels = pretty_labels # in the same order as deep_cols
  ) +
  # scale_y_discrete(label = function(x) stringr::str_trunc(x, 15)) +
  xlim(0, 2.5) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        text = element_text(family = "Arial"))

out_dir <- "../images/FINALS/Pvalue"

# 2) Call svglite() WITHOUT a 'family' argument:
outfile_pdf <- file.path(out_dir, "pvalue.pdf")
cairo_pdf(
  filename = outfile_pdf, width = 9, height = 7, family = "Arial"
)

print(ce_plot)
dev.off()

cat("PDF written to:", normalizePath(outfile_pdf), "\n")