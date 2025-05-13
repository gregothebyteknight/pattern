
# SETTING UP THE ENVIRONMENT
setwd(this.path::here())
source("visual.r")
source("mimic.r")

suppressPackageStartupMessages({
  library(dplyr) # for data manipulation
  library(purrr) # for mapping
})

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
  "compute two‐sided mid‐p on parameter values:
   for each (data, cell) group
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, num, ce
   @par_name: the name of the parameter (par_1 or par_2)"
  cell <- plane_vals <- space_val <- NULL
  big_r <- big_t <- f <- p <- NULL
  tbl_clean <- tbl_all |>
    mutate(par = as.numeric(.data[[par_name]])) |>
    filter(!is.na(par))

  paired <- tbl_clean |>
    group_by(data, cell) |>
    summarise(plane_vals = list(par[dim == "plane"]),
              space_val = mean(par[dim == "space"]),
              .groups = "drop")

  pval_tbl <- paired |>
    rowwise() |>
    mutate(n = length(plane_vals),
           big_r = sum(plane_vals <= space_val),
           big_t = sum(plane_vals == space_val),
           # mid‐p fraction
           f = (big_r - big_t) / n + 0.5 * big_t / n,
           # two‐sided p
           p = 2 * pmin(f, 1 - f), p.label = paste0("p = ", signif(p, 3)),
           y.position = 2.5, group1 = "space", group2 = "plane") |>
    ungroup() |>
    select(data, cell, p, y.position, group1, group2, p.label) # nolint

  pval_tbl
}

reveal_pval_z <- function(tbl_all, par_name) {
  "compute two‐sided mid‐p on z‐scores:
   z = (space_md - plane_md) / plane_sd
   for each (data, cell) group
   @tbl_all: the data frame containing the ce values,
   with the columns: data, cell, dim, num, ce
   @par_name: the name of the parameter (par_1 or par_2)"
  cell <- plane_vals <- space_md <- plane_md <- f <- p <- NULL
  z_space <- plane_sd <- big_r <- big_t <- NULL
  # Clean and cast
  tbl_clean <- tbl_all |>
    mutate(par = as.numeric(.data[[par_name]])) |>
    filter(!is.na(par))

  # Summarise into one row per (data, cell)
  paired <- tbl_clean |>
    group_by(data, cell) |>
    summarise(# list of raw 2D vals
              plane_vals = list(par[dim == "plane"]),
              # plane median & sd
              plane_md = median(par[dim == "plane"]),
              plane_sd = sd(par[dim == "plane"]),
              # single 3D summary
              space_md = median(par[dim == "space"]),
              .groups = "drop")

  # 3) Rowwise compute z‐scores & mid‐p
  pval_tbl <- paired |>
    # drop any zero‐sd or missing
    filter(!is.na(plane_sd), plane_sd > 0, !is.na(space_md)) |>
    rowwise() |>
    mutate(z_space = (space_md - plane_md) / plane_sd,
           # number of plane observations
           n = length(unlist(plane_vals)),
           # mid-p counts
           big_r = sum(((unlist(plane_vals) - plane_md) / plane_sd) <= z_space),
           big_t = sum(((unlist(plane_vals) - plane_md) / plane_sd) == z_space),
           # canonical mid-p fraction and two-sided p
           f = (big_r - 0.5 * big_t) / (n + 1),
           p = 2 * min(f, 1 - f), p.label = paste0("p = ", signif(p, 3)),
           group1 = "space", group2 = "plane",
           #  y.position = pmax(space_md, plane_md + 2.5 * plane_sd))
           y.position = 1e4) |>
    ungroup() |>
    # select exactly your desired columns
    select(data, cell, p, y.position, group1, group2, p.label) # nolint

  pval_tbl
}

# MAIN FUNCTIONS
sim <- FALSE
mode <- "pcf" # choose mode: "naive" or "pcf"
type <- "gamma" # choose type if mode == "pcf": "gamma" or "oscillation"
# and "ce" or "naive" if mode == "naive"
par <- "par_2" # choose parameter to compare

tbl_all <- reveal_all(mode, type, fallback, sim) # get final table
write.csv(tbl_all, file = "../temp/tbl_all.csv", row.names = FALSE)

tbl_all <- read.csv("../temp/tbl_all.csv") # load the table
pval_tbl <- reveal_pval_z(tbl_all, par) # get p-values
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

plot_space(tbl_all, pval_tbl, par, labels = lbl, lim = NA)