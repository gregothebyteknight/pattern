
# INSTALLING LIBRARIES & ARGUMENTS PARSING
suppressPackageStartupMessages({
  library(Hmisc) # for weighted variance
  library(scales) # for col_numeric
  library(fields) # for image.plot (colorbar)
  library(ggplot2) # visualization
  library(dplyr) # for data manipulation
  library(ggpubr) # for stat_compare_means
  library(RColorBrewer) # for color palettes
})
setwd(this.path::here())

mean_pcf <- function(pcf_true, r_grid, pcf_mat, cell_type, num_cells_list) {
  "visualize the mean pcf with ±1 SD
   @pcf_true: the true pcf vector
   @r_grid: the grid of r values
   @pcf_mat: the matrix of pcf for each slice
   @cell_type: the cell type
   @num_cells_list: the list of number of cells in each slice"
  # Compute the mean and SD across all curves at each common r value
  mean_pcf_w <- apply(pcf_mat, 1, function(x) {
    weighted.mean(x, weights = num_cells_list, na.rm = TRUE)
  })
  sd_pcf_w <- apply(pcf_mat, 1, function(x) {
    sqrt(Hmisc::wtd.var(x, weights = num_cells_list, na.rm = TRUE))
  })
  # Plot the true PCF, then overlay the mean PCF with ±1 standard deviation
  png(filename = sprintf("../images/pc_mean_%s.png", cell_type))

  plot(pcf_true, col = "#84a98c", lwd = 2,
       xlab = "r (Distance)", ylab = "pcf", main = "Mean PCF with ±1 SD")

  # Addition of dashed lines for mean ± SD
  lines(r_grid, mean_pcf_w + sd_pcf_w, col = "gray", lty = 2)
  lines(r_grid, mean_pcf_w - sd_pcf_w, col = "gray", lty = 2)

  # Addition of a shaded region between mean + SD and mean - SD
  polygon(c(r_grid, rev(r_grid)), c(mean_pcf_w + sd_pcf_w,
                                    rev(mean_pcf_w - sd_pcf_w)),
          col = rgb(0.1, 0.1, 0.8, 0.2), border = NA)

  # Addition of the true PCF for reference
  lines(r_grid, mean_pcf_w, type = "l", lwd = 2, col = "#7373c9")
  dev.off()
}

cumul_pcf <- function(pcf_true, pcf_list, cell_type, num_cells_list, n_range) {
  "visualize the true PCF against all the pcf for plain slices
   @pcf_true: the true pcf vector
   @pcf_list: the list of pcf for each slice
   @cell_type: the cell type
   @num_cells_list: the list of number of cells in each slice
   @n_range: arbitrary (preselected) range of number of cells"
  # Plot the true PCF
  png(filename = sprintf("../images/pc_total_%s_%s.png", cell_type, n_range[1]))
  plot(pcf_true, main = "Pair Correlation Functions for Varying Angles",
       xlab = "r (Distance)", ylab = "pcf", col = "#84a98c")

  # Prepare continuous color panel
  col_fun <- scales::col_numeric(palette = viridis::cividis(100),
                                 domain = range(num_cells_list, na.rm = TRUE))
  curve_colors <- col_fun(num_cells_list)

  # Overlay the rest of the pcf curves
  # Only plot slices whose num_cells lie in n_range (or all if n_range[1]==-1)
  if (n_range[1] == -1) {
    valid_idx <- seq_along(pcf_list)
  } else {
    valid_idx <- which(num_cells_list >= n_range[1] &
                         num_cells_list <= n_range[2])
  }
  for (i in valid_idx) {
    points(pcf_list[[i]]$r, pcf_list[[i]]$pcf,
           type = "l", col = curve_colors[i])
  }

  # Adjust zlim if num_cells_list has no variation
  z_range <- range(num_cells_list, na.rm = TRUE)
  if (z_range[1] == z_range[2]) {
    z_range <- z_range + c(-1, 1)  # add a small margin
  }

  # Saving the plot and session info
  fields::image.plot(legend.only = TRUE, zlim = z_range,
                     col = viridis::cividis(100),
                     legend.lab = "Number of cells",
                     legend.line = 3, legend.mar = 10, legend.shrink = 0.4)
  dev.off()
}

plot_space <- function(tbl_all, pval_tbl = NULL, param, labels, lim = NA) {
  "creates boxplot to compare difference between spatial
   treats of 2D slices and full 3D
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, num, ce"
  deep_cols <- c(
    "#4C72B0",  # deep[1]: blue
    "#55A868",  # deep[2]: green
    "#C44E52",  # deep[3]: red
    "#8172B2",  # deep[4]: purple
    "#CCB974",  # deep[5]: mustard
    "#64B5CD"   # deep[6]: cyan
  )
  dodge_width <- 0.8

  p <- ggplot(tbl_all, aes_string(x = "cell", y = param, fill = "dim")) +
    geom_boxplot(alpha = 0.7, position = position_dodge(width = dodge_width)) +
    geom_point(inherit.aes = FALSE,
      data = subset(tbl_all, dim == "space"),
      aes(x = cell, y = !!sym(param)),
      color    = "#55A868",
      position = position_nudge(x = 0.2),
      size = 2
    ) +
    facet_wrap(~ data, ncol = 2, scales = "free_y") +
    geom_hline(yintercept = 1, col = "grey", linetype = "dashed") +
    scale_fill_manual(
      name   = "Dimension",
      values = deep_cols[1:2],
      labels = c("2D slice", "3D full")
    ) +
    labs(x = "Cell type", y = labels$y_label,
         fill = "Dimention",
         title = labels$title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, vjust = 0.5)) +
    scale_x_discrete(label = function(x) stringr::str_trunc(x, 15)) +
    scale_y_log10(limits = c(NA, 1e4)) +
    coord_flip()

  # Saving the plot
  cairo_pdf(filename = sprintf("../images/FINALS/Boxplots/%s_box.pdf",
                               labels$file),
            width = 8, height = 6, family = "Arial")
  print(p)
  dev.off()
}