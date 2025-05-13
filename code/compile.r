


library(ggplot2)
library(dplyr)
ce <- read.csv("../temp/tbl_all.csv")

ce <- ce |>
  group_by(data, cell, dim) |>
  summarise(
    mu = mean(par_2, na.rm = TRUE),
    .groups = "drop"
  )

ce_table <- read.csv("../temp/pval_tbl.csv")

pretty_data_names <- c(
  nat_blood  = "LVIBloodBreastCancerModel",
  nat_mainher2 = "MainHer2BreastCancerModel",
  nat_lymph = "LVILymphBreastCancerModel",
  nat_secher2 = "SecondHer2BreastCancerModel",
  sci_brain = "Mouse hippocampal slice",
  sci_embryo = "Intact E6.5-7.0 mouse embryo"
)

ce_table <- ce_table |>
  mutate(p_log = -log(p, base = 10)) |>
  select(data, cell, p_log, p)

plot <- ggplot(ce, aes(y = cell, x = mu, fill = dim)) +
  geom_col(width = 1, position = "dodge2") +
  geom_vline(aes(xintercept = 1), color = "grey",
             linetype = "dashed") +
  # Add text labels showing the 'dim' value on each bar
  geom_text(aes(label = dim),
            position = position_dodge2(width = 1),
            size = 3) +
  facet_wrap(~ data, ncol = 2, scales = "free_y",
             labeller = labeller(data = as_labeller(pretty_data_names))) +
  labs(title = "Mean CE index by Cell Type, Data & Dimension",
       y = "Cell Type", x = "Mean of Clark–Evans index (μ)",
       fill  = "Dimention") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")

ce_plot <- ggplot(ce_table, aes(y = cell, x = p_log, fill = data)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  geom_vline(aes(xintercept = 1.3),
             linetype = "dashed") +
  labs(title = "Comparison of p-values between 2D and 3D values",
       y = "Cell Type", x = "-log10(p-value)",
       fill  = "Dataset name") +
  scale_fill_brewer(palette = "Set2") +
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

ggsave("../temp/gamma_scale_compare_bar.png", plot = ce_plot, bg = "white",
       width = 10, height = 6)