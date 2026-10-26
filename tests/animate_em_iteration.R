#!/usr/bin/env Rscript

# Animate the EM trajectory for the low-dispersion, early-dominant scenario.
# The iteration parameters come from the saved workbook. The original ICC sample
# is reconstructed from the simulation seed and checked against a fresh model fit.

args <- commandArgs(trailingOnly = TRUE)

extra_library <- Sys.getenv("CODEX_R_LIB", unset = "")
if (nzchar(extra_library)) {
  .libPaths(c(extra_library, "D:/R-4.4.3/R-4.4.3/library", .libPaths()))
}

source_xlsx <- if (length(args) >= 1) {
  args[[1]]
} else {
  "notes/Simulation_results/Simu_results_s100_r100_i300.xlsx"
}

output_gif <- if (length(args) >= 2) {
  args[[2]]
} else {
  "notes/Simulation_results/EM_low_dispersion_early_dominant_rep1.gif"
}

target_scenario_id <- if (length(args) >= 3) as.integer(args[[3]]) else 5L
target_rep <- if (length(args) >= 4) as.integer(args[[4]]) else 1L

if (is.na(target_scenario_id) || !target_scenario_id %in% 1:9) {
  stop("target_scenario_id must be an integer from 1 to 9.")
}
if (is.na(target_rep) || target_rep < 1 || target_rep > 100) {
  stop("target_rep must be an integer from 1 to 100.")
}

dispersion_index <- ((target_scenario_id - 1L) %/% 3L) + 1L
weight_index <- ((target_scenario_id - 1L) %% 3L) + 1L
mean_sd_values <- list(c(2, 2), c(15, 3), c(5, 10))
route_weight_values <- list(
  c(CP = 0.20, PS = 0.50, PT = 0.20, PQ = 0.10),
  c(CP = 0.40, PS = 0.35, PT = 0.15, PQ = 0.10),
  c(CP = 0.10, PS = 0.35, PT = 0.35, PQ = 0.20)
)
weight_names <- c("baseline", "early_dominant", "heavy_tailed")

target_true_mean <- mean_sd_values[[dispersion_index]][1]
target_true_sd <- mean_sd_values[[dispersion_index]][2]
target_true_weights <- route_weight_values[[weight_index]]
target_weight_name <- weight_names[[weight_index]]

required_packages <- c(
  "dplyr", "fdrtool", "furrr", "future", "gganimate",
  "ggplot2", "gifski", "readxl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    ". Install them before running this script."
  )
}

# readxl cannot open an xlsx whose normalized path contains non-ASCII characters
# in this Windows R setup. Copying to the session temp directory avoids that issue.
source_copy <- tempfile(fileext = ".xlsx")
if (!file.copy(source_xlsx, source_copy, overwrite = TRUE)) {
  stop("Could not copy the source workbook: ", source_xlsx)
}

iteration_history <- readxl::read_excel(source_copy, sheet = "Per_Iteration") |>
  dplyr::filter(
    .data$scenario_id == target_scenario_id,
    .data$true_mean == target_true_mean,
    .data$true_sd == target_true_sd,
    .data$weight_name == target_weight_name,
    .data$rep == target_rep
  ) |>
  dplyr::arrange(.data$iteration)

if (nrow(iteration_history) == 0) {
  stop("The requested scenario was not found in the Per_Iteration sheet.")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(fdrtool)
})

# Loading with devtools normalizes the package path, which fails for the Chinese
# characters in this workspace under the current Windows R locale. Sourcing the
# package's top-level R files uses the same implementation without that conversion.
for (r_file in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  sys.source(r_file, envir = .GlobalEnv)
}

rlnorm_natural <- function(n, mean, sd) {
  sdlog <- sqrt(log(1 + (sd / mean)^2))
  meanlog <- log(mean) - 0.5 * sdlog^2
  stats::rlnorm(n, meanlog, sdlog)
}

generate_sim_data <- function(N, true_mean, true_sd, route_weights) {
  n_routes <- round(N * route_weights)
  cp <- fdrtool::rhalfnorm(
    n_routes[1],
    theta = sqrt(pi / 2) / (sqrt(2) * true_sd)
  )
  ps <- rlnorm_natural(n_routes[2], true_mean, true_sd)
  pt <- rlnorm_natural(n_routes[3], true_mean, true_sd) +
    rlnorm_natural(n_routes[3], true_mean, true_sd)
  pq <- rlnorm_natural(n_routes[4], true_mean, true_sd) +
    rlnorm_natural(n_routes[4], true_mean, true_sd) +
    rlnorm_natural(n_routes[4], true_mean, true_sd)
  round(c(cp, ps, pt, pq))
}

# Reproduce the future_map seed stream used by run_full_simulation(). Each
# furrr call advances the main RNG once, so all preceding scenario calls must
# be made before reconstructing the requested scenario and repetition.
future::plan(future::sequential)
set.seed(2026)
sim_data <- NULL

for (scenario_index in seq_len(target_scenario_id)) {
  reconstructed <- furrr::future_map(
    seq_len(100),
    function(rep_index) {
      if (
        scenario_index == target_scenario_id &&
          rep_index == target_rep
      ) {
        generate_sim_data(
          N = 100,
          true_mean = target_true_mean,
          true_sd = target_true_sd,
          route_weights = unname(target_true_weights)
        )
      } else {
        NULL
      }
    },
    .options = furrr::furrr_options(seed = TRUE)
  )

  if (scenario_index == target_scenario_id) {
    sim_data <- reconstructed[[target_rep]]
  }
}

if (is.null(sim_data) || length(sim_data) != 100) {
  stop("The target ICC sample could not be reconstructed.")
}

# Confirm that the reconstructed ICC sample reproduces the stored trajectory.
fit_check <- si_estim(
  sim_data,
  n = 300,
  dist = "lognormal",
  init = c(mean(sim_data), stats::sd(sim_data))
)

check_columns <- c("iteration", "mean", "sd", "w_1", "w_2", "w_4", "w_6")
stored_check <- as.data.frame(iteration_history[, check_columns])
rerun_check <- as.data.frame(fit_check$history[, check_columns])

if (nrow(stored_check) != nrow(rerun_check)) {
  stop(
    "Reconstructed trajectory has ", nrow(rerun_check),
    " iterations, but the workbook contains ", nrow(stored_check), "."
  )
}

max_parameter_difference <- max(
  abs(as.matrix(stored_check[, -1]) - as.matrix(rerun_check[, -1])),
  na.rm = TRUE
)

if (!is.finite(max_parameter_difference) || max_parameter_difference > 1e-8) {
  stop(
    "The reconstructed ICC sample did not reproduce the saved trajectory. ",
    "Maximum parameter difference: ", signif(max_parameter_difference, 4)
  )
}

initial_state <- iteration_history[1, ]
initial_state$iteration <- 0
initial_state$mean <- mean(sim_data)
initial_state$sd <- stats::sd(sim_data)
initial_state$mu_change <- NA_real_
initial_state$sigma_change <- NA_real_
initial_state$w_change <- NA_real_
initial_state$w_1 <- 0.25
initial_state$w_2 <- 0.25
initial_state$w_4 <- 0.25
initial_state$w_6 <- 0.25

animation_history <- dplyr::bind_rows(initial_state, iteration_history)

component_density <- function(x, mu, sigma, component) {
  if (component == "CP") {
    return(ifelse(
      x >= 0,
      fdrtool::dhalfnorm(
        x,
        theta = sqrt(pi / 2) / (sqrt(2) * sigma)
      ),
      0
    ))
  }

  generation <- c(PS = 1, PT = 2, PQ = 3)[[component]]
  component_mean <- generation * mu
  component_variance <- generation * sigma^2
  sdlog <- sqrt(log(1 + component_variance / component_mean^2))
  meanlog <- log(component_mean) - 0.5 * sdlog^2
  ifelse(x > 0, stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog), 0)
}

x_max <- max(60, max(sim_data) + 3)
x_grid <- seq(0, x_max, length.out = 601)
component_names <- c("CP", "PS", "PT", "PQ")
weight_columns <- c(CP = "w_1", PS = "w_2", PT = "w_4", PQ = "w_6")

density_frames <- lapply(seq_len(nrow(animation_history)), function(row_index) {
  current <- animation_history[row_index, ]
  pieces <- lapply(component_names, function(component) {
    weight <- current[[weight_columns[[component]]]]
    data.frame(
      iteration = current$iteration,
      x = x_grid,
      component = component,
      density = weight * component_density(
        x_grid,
        current$mean,
        current$sd,
        component
      )
    )
  })
  dplyr::bind_rows(pieces)
}) |>
  dplyr::bind_rows()

mixture_frames <- density_frames |>
  dplyr::group_by(.data$iteration, .data$x) |>
  dplyr::summarise(density = sum(.data$density), .groups = "drop")

true_density <- data.frame(
  x = x_grid,
  density = Reduce(
    `+`,
    lapply(component_names, function(component) {
      target_true_weights[[component]] * component_density(
        x_grid,
        mu = target_true_mean,
        sigma = target_true_sd,
        component = component
      )
    })
  )
)

frame_summary <- animation_history |>
  dplyr::transmute(
    iteration = .data$iteration,
    mean = .data$mean,
    sd = .data$sd,
    label = sprintf(
      "mu = %.2f    sigma = %.2f\nCP = %.3f   PS = %.3f\nPT = %.3f   PQ = %.3f",
      .data$mean, .data$sd, .data$w_1, .data$w_2, .data$w_4, .data$w_6
    )
  )

histogram_peak <- max(table(sim_data)) / length(sim_data)
y_max <- max(
  histogram_peak,
  density_frames$density,
  mixture_frames$density,
  true_density$density,
  na.rm = TRUE
) * 1.16

frame_summary <- frame_summary |>
  dplyr::mutate(x_label = x_max * 0.985, y_label = y_max * 0.965)

component_colours <- c(
  CP = "#0072B2",
  PS = "#009E73",
  PT = "#E69F00",
  PQ = "#CC79A7"
)

animation_plot <- ggplot2::ggplot() +
  ggplot2::geom_histogram(
    data = data.frame(icc = sim_data),
    ggplot2::aes(x = .data$icc, y = ggplot2::after_stat(.data$density)),
    binwidth = 1,
    boundary = -0.5,
    fill = "#D9E8F5",
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_line(
    data = true_density,
    ggplot2::aes(x = .data$x, y = .data$density),
    colour = "grey45",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    data = density_frames,
    ggplot2::aes(
      x = .data$x,
      y = .data$density,
      colour = .data$component,
      group = .data$component
    ),
    linewidth = 1.0
  ) +
  ggplot2::geom_line(
    data = mixture_frames,
    ggplot2::aes(x = .data$x, y = .data$density),
    colour = "black",
    linewidth = 1.25
  ) +
  ggplot2::geom_vline(
    xintercept = target_true_mean,
    colour = "grey35",
    linetype = "dashed",
    linewidth = 0.7
  ) +
  ggplot2::geom_vline(
    data = frame_summary,
    ggplot2::aes(xintercept = .data$mean),
    colour = "#D55E00",
    linetype = "dotdash",
    linewidth = 0.8
  ) +
  ggplot2::geom_label(
    data = frame_summary,
    ggplot2::aes(
      x = .data$x_label,
      y = .data$y_label,
      label = .data$label
    ),
    hjust = 1,
    vjust = 1,
    family = "mono",
    size = 4.1,
    lineheight = 1.15,
    label.size = 0.25,
    fill = scales::alpha("white", 0.90)
  ) +
  ggplot2::scale_colour_manual(
    values = component_colours,
    breaks = component_names
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(0, y_max),
    expand = FALSE
  ) +
  ggplot2::labs(
    title = "EM mixture fit: iteration {current_frame}",
    subtitle = paste0(
      "Mean = ", target_true_mean,
      ", SD = ", target_true_sd,
      ", ", gsub("_", " ", target_weight_name),
      ", N = 100, repetition ", target_rep, ". ",
      "Black: fitted mixture; grey dashed: true mixture"
    ),
    x = "ICC interval (days)",
    y = "Density",
    colour = "Weighted component"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 18),
    plot.subtitle = ggplot2::element_text(size = 11, colour = "grey30"),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(15, 18, 12, 15)
  ) +
  gganimate::transition_manual(.data$iteration)

animation <- gganimate::animate(
  animation_plot,
  nframes = nrow(animation_history),
  fps = 10,
  width = 1200,
  height = 720,
  res = 110,
  renderer = gganimate::gifski_renderer(loop = TRUE)
)

gif_temp <- tempfile(fileext = ".gif")
gganimate::anim_save(gif_temp, animation = animation)

output_dir <- dirname(output_gif)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
if (!file.copy(gif_temp, output_gif, overwrite = TRUE)) {
  stop("Could not save the animation to: ", output_gif)
}

message("Animation saved to: ", output_gif)
message(
  "Frames animated: ", nrow(animation_history),
  " (initialization plus ", nrow(iteration_history), " EM iterations)"
)
message("Reconstruction max parameter difference: ", max_parameter_difference)
message(
  "Final weights: CP=", round(tail(iteration_history$w_1, 1), 4),
  ", PS=", round(tail(iteration_history$w_2, 1), 4),
  ", PT=", round(tail(iteration_history$w_4, 1), 4),
  ", PQ=", round(tail(iteration_history$w_6, 1), 4)
)
