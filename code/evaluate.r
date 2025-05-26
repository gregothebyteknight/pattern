
# SETTING UP THE ENVIRONMENT
setwd(this.path::here())
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

ce_table <- read.csv("../temp/pval_tbl.csv")

ce_table <- ce_table |>
  mutate(p_log = -log(p.BH, base = 10)) |>
  select(data, cell, p_log)

ce_plot <- ggplot(ce_table, aes(y = cell, x = p_log, fill = data)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  geom_vline(aes(xintercept = 1.3), linetype = "dashed") +
  labs(title = "Comparison of p-values (FDR-BH) between 2D and 3D values",
       y = "Cell Type", x = "-log10(adj p-value)", fill  = "Dataset name") +
  scale_fill_brewer(palette = "Set2") +
  xlim(0, 3.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_discrete(breaks = c("nat_blood", "nat_mainher2", "nat_lymph",
                                 "nat_secher2", "sci_brain", "sci_embryo"),
                      labels = c("LVIBloodBreastCancerModel",
                                 "MainHer2BreastCancerModel",
                                 "LVILymphBreastCancerModel",
                                 "SecondHer2BreastCancerModel",
                                 "Mouse hippocampal slice",
                                 "Intact E6.5-7.0 mouse embryo"))

ggsave("../images/FINALS/Pvalue/compare_bar.png", plot = ce_plot, bg = "white",
       width = 10, height = 6)