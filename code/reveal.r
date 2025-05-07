
# SETTING UP THE ENVIRONMENT
setwd(this.path::here())
source("visual.r")
source("mimic.r")

suppressPackageStartupMessages(library(dplyr)) # for data manipulation

# Variables
fallback <- list(gamma = c(height = 15, shape = 3, scale = 4),
                 oscillation = c(amp = 15, decay = 1, rate = 4))
same <- "2D Slices vs. Whole 3D Tissue"
labs <- list(
  naive = list(naive = list(y_label = "PCF value (%s)", file = "naive_%s",
                            title = "PCF function (%s): 2D Slices vs. Full 3D"),
               ce = list(y_label = "Clark-Evans Index", file = "ce",
                         title = sprintf("Clark–Evans Index: %s", same))),
  pcf = list(gamma = list(y_label = "PCF Fit with Gamma: %s",
                          file = "pcf_gamma_%s",
                          title = sprintf("PCF Fit on Gamma: %s", same)),
             oscillation = list(y_label = "PCF Fit with Oscillation: %s",
                                file = "pcf_oscillation_%s",
                                title = sprintf("PCF Fit with Ocsillation: %s",
                                                same)))
)

# Declare functions
pcf_mode <- function(cult, dataset, dim, tbl_all, type, fallback) {
  "compute the pcf fit for pcf values
   @cult: the path to the cell type
   @dataset: the name of the dataset
   @dim: the dimension of the pcf (plane or space)
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, par_1, par_2 (2nd and 3rd fit parameters)
   @type: the type of fit (gamma or oscillation)
   @fallback: the fallback parameters for the fit"
  pcf_tbl <- read.csv(file.path(cult, sprintf("pcf_%s.csv", dim)))
  for (i in 2:ncol(pcf_tbl)) {
    if (all(pcf_tbl[, i] <= 1, na.rm = TRUE)) {
      message("Skipping slice ", i, " in ", cult, ": no y > 1")
      next
    }
    params <- fit_func(pcf_tbl$r, pcf_tbl[, i], fallback = fallback, # nolint
                       type = type, show_plot = FALSE)
    tbl_all <- bind_rows(tbl_all, tibble(data = dataset, cell = basename(cult),
                                         dim = dim, par_1 = params[2],
                                         par_2 = params[3])) # returns tbl_all
  }
  tbl_all # return the data frame
}

naive_mode <- function(cult, dataset, tbl_all, type) {
  "concatenate ce values from either plane or space ce values tables
   @cult: the path to the cell type
   @dataset: the name of the dataset
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, 
   par_1 (num of cells or pcf values), par_2 (ce values or r values)
   @type: the type of table (ce or naive)"
  if (type == "ce") name <- "ce_tbl.csv" else name <- "pcf_naive.csv"
  if (type == "ce") par_1 <- "num" else par_1 <- "pcf"
  if (type == "ce") par_2 <- "ce" else par_2 <- "r"

  naive <- read.csv(file.path(cult, name))
  tbl_all <- bind_rows(tbl_all, tibble(data = dataset, cell = basename(cult),
                                       par_1 = naive[[par_1]], dim = naive$dim,
                                       par_2 = naive[[par_2]]))
  tbl_all # return the data frame
}

reveal_all <- function(mode, type, fallback, sim = FALSE) {
  "concatenate spatial functions values from either plane or space
   @mode: the mode of the analysis (pcf or ce)
   @type: the type of fit (gamma or oscillation)"
  tbl_all <- tibble(data = character(), cell = character(), dim = character(),
                    par_1 = double(), par_2 = double()) # initialize dataframe
  # DATA AGGREGATION
  for (dataset in list.dirs(path = sprintf("../data/%s", mode),
                            recursive = FALSE)) {
    if (sim && !grepl("sim_[0-9]+$", basename(dataset))) next
    if (!sim && grepl("sim_[0-9]+$", basename(dataset))) next

    for (cult in list.dirs(path = dataset, recursive = FALSE)) {
      print(sprintf("Cell type:%s", basename(cult)))
      dataset <- basename(dataset) # now need only the dataset name, not path
      if (mode == "naive") {
        tbl_all <- naive_mode(cult, dataset, tbl_all, type)
      } else {
        tbl_all <- pcf_mode(cult, dataset, "plane", tbl_all, type, fallback)
        tbl_all <- pcf_mode(cult, dataset, "space", tbl_all, type, fallback)
      }
    }
  }
  write.csv(tbl_all, file = "../data/tbl_all.csv", row.names = FALSE)
  tbl_all # return the data frame with all the data)
}

reveal_pval <- function(tbl_all, par_name) {
  "compute the p-values for the difference between
   spatial treats of 2D slices and full 3D
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, par_1, par_2"
  cell <- dim <- plane <- space <- p <- NULL # nonsense to avoid check notes
  space_vals <- tbl_all |>
    filter(dim == "space", na.rm = TRUE) |>
    distinct(cell, .keep_all = TRUE) |> # one row per cell
    rename(space = all_of(par_name)) |>
    select(cell, space)

  plane_vals <- tbl_all |>
    filter(dim == "plane", na.rm = TRUE) |>
    rename(plane = all_of(par_name)) |>
    select(cell, plane)

  diff_tbl <- plane_vals |>
    left_join(space_vals, by = "cell", relationship = "many-to-one") |>
    mutate(diff = plane - space)

  pval_tbl <- diff_tbl |>
    group_by(cell) |>
    summarise(p = t.test(diff)$p.value) |>
    mutate(p.label = paste0("p = ", signif(p, 3)),
           group1 = "space", group2  = "plane",
           y.position = max(diff_tbl$plane, na.rm = TRUE) * 20)

  # distribution of diff values
  diff_ggplot <- diff_tbl |>
    ggplot2::ggplot(ggplot2::aes(diff, fill = cell)) +
    ggplot2::geom_histogram() +
    ggplot2::geom_vline(aes(xintercept = 0), color = "#d86666",
                        linetype = "dashed") +
    labs(title = "Distribution of differences",
         x = "Difference (plane - space)",
         y = "Count") +
    ggplot2::facet_wrap(~cell, scales = "free_y") +
    ggplot2::theme_minimal()

  ggplot2::ggsave("../temp/diff_histogram.png",
                  plot = diff_ggplot, bg = "white")

  value_ggplot <- plane_vals |>
    ggplot2::ggplot(ggplot2::aes(plane, fill = cell)) +
    ggplot2::geom_histogram() +
    ggplot2::geom_vline(aes(xintercept = 0), color = "#d86666",
                        linetype = "dashed") +
    labs(title = "Distribution of differences",
         x = "Difference (plane - space)",
         y = "Count") +
    ggplot2::facet_wrap(~cell, scales = "free_y") +
    ggplot2::theme_minimal()

  ggplot2::ggsave("../temp/value_histogram.png",
                  plot = value_ggplot, bg = "white")


  pval_tbl # return the data frame with p-values
}

# MAIN FUNCTIONS
sim <- FALSE
mode <- "naive" # choose mode: "naive" or "pcf"
type <- "ce" # choose type if mode == "pcf": "gamma" or "oscillation"
# and "ce" or "naive" if mode == "naive"
par <- "par_2" # choose parameter to compare

tbl_all <- reveal_all(mode, type, fallback, sim) # get final table
write.csv(tbl_all, file = "../temp/tbl_all.csv", row.names = FALSE)

tbl_all <- read.csv("../temp/tbl_all.csv") # load the table
pval_tbl <- reveal_pval(tbl_all, par) # get p-values
write.csv(pval_tbl, file = "../temp/pval_tbl.csv", row.names = FALSE)

# VISUALIZATION
lbl <- labs[[mode]][[type]]
if (mode == "pcf") {
  if (type == "gamma") {
    if (par == "par_1") f_par <- "shape" else f_par <- "scale"
  } else if (type == "oscillation") {
    if (par == "par_1") f_par <- "decay" else f_par <- "rate"
  }
  lbl$y_label <- sprintf(lbl$y_label, f_par)
  lbl$file <- sprintf(lbl$file, f_par)
} else if (mode == "naive") {
  if (type == "naive") {
    if (par == "par_1") f_par <- "max(pcf)" else f_par <- "max(r)"
    lbl$y_label <- sprintf(lbl$y_label, f_par)
    lbl$title <- sprintf(lbl$title, f_par)
    lbl$file <- sprintf(lbl$file, f_par)
  }
}

plot_space(tbl_all, pval_tbl, par, labels = lbl)