# Run from the repository root or the folder containing the QMD report.
# Required packages: openxlsx, dplyr, tidyr, ggplot2
# Produces: ../Simulation_results/sample_size_effect_rmse_convergence.png

library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)

simulation_dir <- if (dir.exists(file.path("notes", "Simulation_results"))) {
  file.path("notes", "Simulation_results")
} else {
  file.path("..", "Simulation_results")
}
sample_sizes <- c(20, 50, 100, 200)

read_runs <- function(n) {
  path <- file.path(simulation_dir, sprintf("Simu_results_s%d_r100_i300.xlsx", n))
  if (!file.exists(path)) stop("Missing input workbook: ", path)
  dat <- openxlsx::read.xlsx(path, sheet = "Per_Repetition")
  required <- c("scenario_id", "true_mean", "true_sd", "weight_name",
                "est_mean", "est_sd", "converged")
  if (!all(required %in% names(dat))) {
    stop("Missing columns in ", path, ": ",
         paste(setdiff(required, names(dat)), collapse = ", "))
  }
  dat |>
    mutate(
      N = n,
      converged = as.character(converged) %in% c("TRUE", "T", "1")
    )
}

runs <- bind_rows(lapply(sample_sizes, read_runs))

rmse_finite <- function(est, truth) {
  err2 <- (est - truth)^2
  err2 <- err2[is.finite(err2)]
  if (!length(err2)) return(NA_real_)
  sqrt(mean(err2))
}

summary <- runs |>
  group_by(N, scenario_id, true_mean, true_sd, weight_name) |>
  summarise(
    n_reps = n(),
    n_converged = sum(converged),
    mean_rmse = rmse_finite(est_mean[converged], dplyr::first(true_mean)),
    sd_rmse = rmse_finite(est_sd[converged], dplyr::first(true_sd)),
    convergence_pct = 100 * mean(converged),
    valid_mean = sum(converged & is.finite(est_mean)),
    valid_sd = sum(converged & is.finite(est_sd)),
    .groups = "drop"
  )|>
  mutate(
    Dispersion = case_when(
      true_mean == 15 & true_sd == 3 ~ "Low (15, 3)",
      true_mean == 2 & true_sd == 2 ~ "Moderate (2, 2)",
      true_mean == 5 & true_sd == 10 ~ "High (5, 10)",
      TRUE ~ NA_character_
    ),
    Profile = recode(weight_name,
      baseline = "Baseline",
      early_dominant = "CP-heavy",
      heavy_tailed = "Higher-generation-heavy"
    ),
    Dispersion = factor(Dispersion,
      levels = c("Low (15, 3)", "Moderate (2, 2)", "High (5, 10)")),
    Profile = factor(Profile,
      levels = c("Baseline", "CP-heavy", "Higher-generation-heavy"))
  )

if (nrow(summary) != 36L) stop("Expected 36 combinations (4 sample sizes x 9 scenarios)")
if (anyNA(summary$Dispersion)) stop("Unexpected true mean/SD combination")
if (any(summary$n_reps != 100L)) warning("Some scenario groups do not contain 100 repetitions")
if (any(summary$valid_mean < summary$n_converged |
        summary$valid_sd < summary$n_converged)) {
  warning("Non-finite converged estimates found; only finite converged estimates are included in RMSE")
}

plot_data <- summary |>
  pivot_longer(
    cols = c(mean_rmse, sd_rmse, convergence_pct),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(
    Metric = factor(Metric,
      levels = c("mean_rmse", "sd_rmse", "convergence_pct"),
      labels = c("SI mean RMSE (days)", "SI SD RMSE (days)", "Convergence (%)"))
  )

# Free y-axis scales by metric row. All three dispersions in the same row
# share a scale, making their between-scenario values comparable.
p <- ggplot(plot_data, aes(x = N, y = Value, colour = Profile, group = Profile)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.4, na.rm = TRUE) +
  facet_grid(Metric ~ Dispersion, scales = "free_y", switch = "y") +
  scale_x_continuous(breaks = sample_sizes, limits = c(15, 205)) +
  scale_color_manual(values = c(
    "Baseline" = "#68BBA5",
    "CP-heavy" = "#F08F68",
    "Higher-generation-heavy" = "#839CC9"
  )) +
  labs(
    title = "Effect of sample size on SI parameter recovery and EM convergence",
    subtitle = "Nine scenarios across N = 20, 50, 100 and 200; 100 repetitions per setting",
    x = "ICC sample size (N)", y = NULL, colour = "Pathway profile",
    caption = paste0(
      "RMSE: finite converged estimates only. ",
      "Convergence: all repetitions. Maximum EM iterations = 300."
    )
  ) +
  theme_bw(base_size = 15) +
  theme(
    strip.background = element_rect(fill = "#F0F1F3", colour = NA),
    strip.text = element_text(face = "bold", size = 14),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 14),
    plot.caption = element_text(size = 11),
    panel.spacing = grid::unit(0.8, "lines")
  )

output <- file.path(simulation_dir, "sample_size_effect_rmse_convergence.png")
ggsave(output, plot = p, width = 16, height = 11.5, units = "in", dpi = 300,
       bg = "white", device = grDevices::png, type = "cairo")
message("Saved trend plot: ", normalizePath(output))

# Optional check table: avoids accidentally comparing RMSE from different subsets.
write.csv(summary,
          file.path(simulation_dir, "sample_size_effect_summary.csv"),
          row.names = FALSE)
