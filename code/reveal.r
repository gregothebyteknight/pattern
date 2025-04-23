

# INSTALL LIBRARIES
library(sf) # for spatial data
library(FNN) # fast nearest neighbor search
library(geometry)

setwd(this.path::here())

source("visual.r")

suppressPackageStartupMessages({
  library(Hmisc) # for weighted variance
  library(scales) # for col_numeric
  library(fields) # for image.plot (colorbar)
  library(ggplot2) # visualization
  library(dplyr) # for data manipulation
  library(RColorBrewer) # for color palettes
})

setwd(this.path::here())

ce_all <- tibble(data = character(), cell = character(),
                 dim = character(), num = integer(), ce = double())

for (dataset in list.dirs(path = "../data/ce_temp", recursive = FALSE)) {
  for (cult in list.dirs(path = dataset, recursive = FALSE)) {
    ce_temp <- read.csv(file.path(cult, "ce_tbl.csv"))
    ce_all <- bind_rows(ce_all, tibble(dim = ce_temp$dim, num = ce_temp$num,
                                       ce = ce_temp$ce, cell = basename(cult),
                                       data = basename(dataset)))
  }
}

plot_ce(ce_all)