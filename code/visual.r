
# INSTALLING LIBRARIES & ARGUMENTS PARSING
library(Hmisc) # for weighted variance
library(scales) # for col_numeric
library(fields) # for image.plot (colorbar)
library(ggplot2) # visualization
library(RColorBrewer) # for color palettes

setwd(this.path::here())

mean_pcf <- function(pcf_true, r_grid, pcf_mat, cell_type, num_cells_list) {
  # Compute the mean and SD across all curves at each common r value
  mean_pcf_w <- apply(pcf_mat, 1, function(x) {
    weighted.mean(x, weights = num_cells_list, na.rm = TRUE)
  })
  sd_pcf_w <- apply(pcf_mat, 1, function(x) {
    sqrt(Hmisc::wtd.var(x, weights = num_cells_list, na.rm = TRUE))
  })
  # Plot the true PCF, then overlay the mean PCF with ±1 standard deviation
  png(filename = sprintf("../images/pc_mean_%s.png", cell_type),
      width = 800, height = 600)

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

cumul_pcf <- function(pcf_true, pcf_list, cell_type, num_cells_list) {
  # Plot the true PCF
  png(filename = sprintf("../images/pc_total_%s.png", cell_type))
  plot(pcf_true, main = "Pair Correlation Functions for Varying Angles",
       xlab = "r (Distance)", ylab = "pcf", col = "#84a98c")

  # Prepare continuous color panel
  col_fun <- scales::col_numeric(palette = viridis::cividis(100),
                                 domain = range(num_cells_list, na.rm = TRUE))
  curve_colors <- col_fun(num_cells_list)

  # Overlay the rest of the pcf curves
  for (i in seq_along(pcf_list)) {
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
                     legend.line = 2, legend.mar = 3)
  dev.off()
}

plot_ce <- function(ce_all) {
  png(filename = sprintf("../images/pc_total_%s.png"))

  ggplot(ce_all, aes(x = ce_all$dimension, y = ce_all$ce,
                     fill = ce_all$cell_type)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = ce_all$dataset),
                width = 0.2, size = 1.5, alpha = 0.8) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    labs(title = "Clark–Evans Index: 2D Slices vs. Full 3D",
         x = "Dimension", y = "CE Index",
         fill = "Cell Type", color = "Dataset") +
    theme_minimal()
}