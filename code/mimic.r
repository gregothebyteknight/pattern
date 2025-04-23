
# SETUP THE ENVIRONMENT
suppressPackageStartupMessages({
  library(minpack.lm)
  library(dplyr)
})

setwd(this.path::here())

"
Optimal (maybe) model parameters:
  - gamma: height = 15, shape = 3, scale = 4
  - oscillation: amp = 15, decay = 1, rate = 4
"

model_lib <- list(
  gamma =  function(x, pars) {
    pars["height"] / (pars["scale"]^pars["shape"] * gamma(pars["shape"])) *
      (x^(pars["shape"] - 1)) * exp(-x / pars["scale"]) + 1
  },
  oscillation = function(x, pars) {
    pars["amp"] * exp(-pars["decay"] * x) * sin(pars["rate"] * x) + 1
  }
)

fit_func <- function(x, y, type = "gamma", pars = NULL, show_plot = FALSE) {
  row_n <- grp <- n <- NULL # nonsense to avoid R CMD check notess
  df <- data.frame(x = x, y = y) |>
    dplyr::mutate(row_n = dplyr::row_number()) |>
    dplyr::filter(y > 1) |>
    dplyr::mutate(grp = row_n - dplyr::row_number()) |>
    dplyr::group_by(grp) |>
    dplyr::mutate(n = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(n == max(n)) |>
    dplyr::mutate(x = x - min(x)) |>
    dplyr::select(-grp, -row_n, -n)

  if (nrow(df) == 0) stop("No data after filtering.")

  # Model-specific parameters
  if (is.null(pars)) { # If pars is not NULL, use it as initial values
    if (type == "gamma") {
      y_base <- df$y - 1 # pcf is equal to dgamma + 1
      m <- mean(df$x * y_base) / mean(y_base) # weighted mean
      v <- sum(y_base * (df$x - m)^2) / sum(y_base) # weighted variance
      shape_init <- m^2 / v # From: shape = μ² / σ²
      scale_init <- v / m # From: scale = σ² / μ
      pars <- c(height = max(y_base), shape = shape_init, scale = scale_init)
    } else if (type == "oscillation") {
      amp <- (max(df$y) - min(df$y)) / 2
      rate <- 2 * pi / (df$x[which.max(df$y)] - df$x[which.max(df$y) - 1])
      pars <- c(amp = amp, decay = 1, rate = rate)
    }
  }
  print(paste0("Initial parameters: ", paste(pars, collapse = " ")))

  fit <- nls.lm(par = pars, fn = function(pars, x, y, model) y - model(x, pars),
    x = df$x, y = df$y, model = model_lib[[type]],
    control = nls.lm.control(maxiter = 1000), lower = c(0, 0, 0)
  )

  if (show_plot) {
    par(las = 1, bty = "l")
    plot(df$x, df$y, xlab = "r", ylab = "Pcf(r)", ylim = c(0, max(df$y) * 1.1))
    # draw the fitted curve directly using your model and fitted pars
    lines(df$x, model_lib[[type]](df$x, fit$par), col = "red", lwd = 2, lty = 2)
    # add the fitted parameters to the plot
    legend("topright", legend = c(
      paste0(sprintf("%s = ", names(fit$par)[1]), round(fit$par[1], 3)),
      paste0(sprintf("%s = ", names(fit$par)[2]), round(fit$par[2], 3)),
      paste0(sprintf("%s = ", names(fit$par)[3]), round(fit$par[3], 3))
    ), bty = "n", cex = 0.8)
    title(sprintf("Model Fit (%s) to Pair Correlation Function (Pcf)", type))
  }
  return(fit)
}

# Load data
pcf_plain <- read.csv("../data/pcf/sci_embryo/Cdh5+ Endothelium cells/pcf_plain.csv") # nolint
pcf_space <- read.csv("../data/pcf/sci_embryo/Cdh5+ Endothelium cells/pcf_space.csv") # nolint

# Fit models (specify type if not "oscillation")
pars <- c(height = 50, shape = 3, scale = 5) # if some understanding exists
fit <- fit_func(x = pcf_plain$r, y = pcf_plain$mean, type = "oscillation",
                pars = NULL, show_plot = TRUE)
print(fit)
