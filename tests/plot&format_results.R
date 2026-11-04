library(ggplot2)
library(dplyr)
library(readxl)
library(openxlsx)
library(tools)
library(cowplot)

profile_colours <- c(
  "Baseline" = "#66C2A5",
  "CP-heavy" = "#FC8D62",
  "Higher-gen-heavy" = "#8DA0CB"
)

# ==========================================
# 1. Helper function for plotting
# ==========================================
plot_scatter_bar <- function(data, summary_data, var_bias, var_mean_bias, var_lower, var_upper, title) {
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "#E41A1C", linewidth = 0.7) +
    geom_point(
      data = data,
      aes(x = as.numeric(weight_name) - 0.15, y = !!sym(var_bias), color = weight_name),
      position = position_jitter(width = 0.08, height = 0),
      alpha = 0.35,
      size = 1.0
    ) +
    geom_pointrange(
      data = summary_data,
      aes(x = as.numeric(weight_name) + 0.15,
          y = !!sym(var_mean_bias),
          ymin = !!sym(var_lower),
          ymax = !!sym(var_upper)),
      color = "black",
      fill = "black",
      shape = 21,
      size = 0.8,
      linewidth = 0.9
    ) +
    facet_wrap(~ dispersion, ncol = 1, scales = "free_y") + 
    scale_x_continuous(
      breaks = c(1, 2, 3),
      labels = c("Baseline", "CP-heavy", "Higher-gen-\nheavy"),
      limits = c(0.5, 3.5) 
    ) +
    scale_y_continuous(n.breaks = 6) +
    scale_color_manual(values = profile_colours) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 14),
      strip.background = element_rect(fill = "#F3F4F6", color = NA),
      axis.text = element_text(color = "#374151", size = 12.5),
      axis.text.x = element_text(size = 11.5),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 8)),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.45, "lines"),
      plot.title = element_text(face = "bold", size = 20, margin = margin(b = 8)),
      plot.margin = margin(7, 7, 7, 7)
    ) +
    labs(
      title = title,
      y = "Estimation Bias (est - true)"
    )
}

# ==========================================
# 2. Helper function for Excel table formatting
# ==========================================
format_for_excel <- function(data) {
  bind_rows(
    data %>%
      transmute(
        Scenario = scenario_id, Profile = profile, Parameter = "Mean",
        `Avg. estimate` = round(avg_est_mean, 2), Bias = round(bias_est_mean, 2),
        MSE = round(mse_est_mean, 2), RMSE = round(sqrt(mse_est_mean), 2),
        `Convergence (%)` = round(100 * convergence_rate, 0)
      ),
    data %>%
      transmute(
        Scenario = scenario_id, Profile = profile, Parameter = "SD",
        `Avg. estimate` = round(avg_est_sd, 2), Bias = round(bias_est_sd, 2),
        MSE = round(mse_est_sd, 2), RMSE = round(sqrt(mse_est_sd), 2),
        `Convergence (%)` = round(100 * convergence_rate, 0)
      )
  ) %>%
    arrange(Scenario, factor(Parameter, levels = c("Mean", "SD")))
}

plot_convergence_rate <- function(data) {
  convergence_data <- data %>%
    group_by(dispersion, weight_name) %>%
    summarise(
      n_converged = sum(converged %in% TRUE, na.rm = TRUE),
      n_repetitions = n(),
      convergence_rate = 100 * n_converged / n_repetitions,
      .groups = "drop"
    ) %>%
    mutate(
      profile_bar = factor(
        as.character(weight_name),
        levels = c("Higher-gen-heavy", "CP-heavy", "Baseline")
      ),
      rate_label = sprintf("%.0f%%", convergence_rate),
      label_x = if_else(convergence_rate >= 88, convergence_rate - 2, convergence_rate + 2),
      label_hjust = if_else(convergence_rate >= 88, 1, 0)
    )

  ggplot(convergence_data, aes(x = convergence_rate, y = profile_bar, fill = weight_name)) +
    geom_col(width = 0.62, alpha = 0.58) +
    geom_text(
      aes(x = label_x, label = rate_label, hjust = label_hjust),
      colour = "#26323A",
      fontface = "bold",
      size = 4.5
    ) +
    facet_wrap(~dispersion, ncol = 1) +
    scale_fill_manual(values = profile_colours, guide = "none") +
    scale_x_continuous(
      limits = c(0, 108),
      breaks = c(0, 50, 100),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = expansion(add = 0.45)) +
    labs(
      title = "Convergence rate",
      x = "Converged repetitions (%)",
      y = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 14),
      strip.background = element_rect(fill = "#F3F4F6", color = NA),
      panel.border = element_rect(fill = NA, colour = "#E5E7EB", linewidth = 0.35),
      panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.35),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.y = grid::unit(0.45, "lines"),
      axis.text = element_text(colour = "#374151", size = 12.5),
      axis.text.y = element_text(size = 11.5),
      axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 8)),
      plot.title = element_text(face = "bold", size = 20, margin = margin(b = 8)),
      plot.margin = margin(7, 7, 7, 7)
    )
}

# ==========================================
# 3. Main wrapper function for report generation
# ==========================================
generate_simulation_report <- function(
  input_path,
  output_dir,
  generate_tables = FALSE,
  combined_output_dir = file.path(output_dir, "Combined_summary_panels"),
  save_individual_plots = FALSE,
  save_combined = TRUE
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(combined_output_dir, recursive = TRUE, showWarnings = FALSE)

  base_name <- file_path_sans_ext(basename(input_path))
  metadata <- regmatches(
    base_name,
    regexec("_s([0-9]+)_r([0-9]+)_i([0-9]+)$", base_name)
  )[[1]]
  sample_size <- if (length(metadata) == 4) metadata[2] else "unknown"
  repetitions <- if (length(metadata) == 4) metadata[3] else "unknown"
  iteration_limit <- if (length(metadata) == 4) metadata[4] else "unknown"

  all_runs <- read_excel(input_path, sheet = "Per_Repetition") %>%
    mutate(
      bias_mean = est_mean - true_mean,
      bias_sd = est_sd - true_sd,
      dispersion = case_when(
        true_mean == 15 & true_sd == 3  ~ "Low dispersion (15, 3)",
        true_mean == 2  & true_sd == 2  ~ "Moderate dispersion (2, 2)",
        true_mean == 5  & true_sd == 10 ~ "High dispersion (5, 10)",
        TRUE ~ paste0("m=", true_mean, ", s=", true_sd)
      ),
      dispersion = factor(
        dispersion,
        levels = c(
          "Low dispersion (15, 3)",
          "Moderate dispersion (2, 2)",
          "High dispersion (5, 10)"
        )
      ),
      weight_name = factor(
        weight_name,
        levels = c("baseline", "early_dominant", "heavy_tailed"),
        labels = c("Baseline", "CP-heavy", "Higher-gen-heavy")
      )
    )

  converged_runs <- all_runs %>% filter(converged %in% TRUE)
  summary_stats <- converged_runs %>%
    group_by(dispersion, weight_name) %>%
    summarise(
      n_converged = n(),
      mean_bias_mean = mean(bias_mean, na.rm = TRUE),
      p05_bias_mean = quantile(bias_mean, probs = 0.05, na.rm = TRUE),
      p95_bias_mean = quantile(bias_mean, probs = 0.95, na.rm = TRUE),
      mean_bias_sd = mean(bias_sd, na.rm = TRUE),
      p05_bias_sd = quantile(bias_sd, probs = 0.05, na.rm = TRUE),
      p95_bias_sd = quantile(bias_sd, probs = 0.95, na.rm = TRUE),
      .groups = "drop"
    )

  p_mean <- plot_scatter_bar(
    converged_runs,
    summary_stats,
    var_bias = "bias_mean",
    var_mean_bias = "mean_bias_mean",
    var_lower = "p05_bias_mean",
    var_upper = "p95_bias_mean",
    title = "Mean bias"
  )
  p_sd <- plot_scatter_bar(
    converged_runs,
    summary_stats,
    var_bias = "bias_sd",
    var_mean_bias = "mean_bias_sd",
    var_lower = "p05_bias_sd",
    var_upper = "p95_bias_sd",
    title = "Standard-deviation bias"
  )
  p_convergence <- plot_convergence_rate(all_runs)

  if (save_individual_plots) {
    ggsave(
      file.path(output_dir, paste0(base_name, "_mean_summary.jpg")),
      p_mean,
      width = 7.4,
      height = 10.5,
      dpi = 300,
      bg = "white"
    )
    ggsave(
      file.path(output_dir, paste0(base_name, "_sd_summary.jpg")),
      p_sd,
      width = 7.4,
      height = 10.5,
      dpi = 300,
      bg = "white"
    )
  }

  three_panel_body <- cowplot::plot_grid(
    p_mean,
    p_sd,
    p_convergence,
    nrow = 1,
    rel_widths = c(1, 1, 0.9),
    align = "h",
    axis = "tb"
  )
  header <- cowplot::ggdraw() +
    cowplot::draw_label(
      sprintf("Simulation summary - sample size N = %s", sample_size),
      x = 0.01,
      y = 0.72,
      hjust = 0,
      fontface = "bold",
      size = 28,
      colour = "#1F2937"
    ) +
    cowplot::draw_label(
      sprintf(
        "Bias panels use converged runs only; convergence uses all %s repetitions | %s-iteration limit",
        repetitions,
        iteration_limit
      ),
      x = 0.01,
      y = 0.20,
      hjust = 0,
      size = 15,
      colour = "#5E6874"
    )
  combined_plot <- cowplot::plot_grid(
    header,
    three_panel_body,
    ncol = 1,
    rel_heights = c(0.10, 1)
  )

  combined_path <- file.path(
    combined_output_dir,
    paste0(base_name, "_mean_sd_convergence.jpg")
  )
  if (save_combined) {
    ggsave(
      combined_path,
      combined_plot,
      width = 22,
      height = 12.5,
      dpi = 300,
      bg = "white",
      device = grDevices::jpeg,
      type = "cairo",
      limitsize = FALSE
    )
  }

  # Give each dispersion level its own three-panel figure for legible report use.
  dispersion_files <- c(
    "Low dispersion (15, 3)" = "low_dispersion",
    "Moderate dispersion (2, 2)" = "moderate_dispersion",
    "High dispersion (5, 10)" = "high_dispersion"
  )
  dispersion_paths <- character(0)
  for (level in names(dispersion_files)) {
    level_runs <- filter(all_runs, as.character(dispersion) == level)
    level_converged <- filter(converged_runs, as.character(dispersion) == level)
    level_stats <- filter(summary_stats, as.character(dispersion) == level)

    mean_panel <- plot_scatter_bar(
      level_converged, level_stats, "bias_mean", "mean_bias_mean",
      "p05_bias_mean", "p95_bias_mean", "Mean bias"
    ) +
      labs(y = "Mean bias (est - true)") +
      theme(
        strip.text = element_blank(), strip.background = element_blank(),
        axis.text = element_text(size = 17), axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 18), plot.title = element_text(size = 24)
      )
    sd_panel <- plot_scatter_bar(
      level_converged, level_stats, "bias_sd", "mean_bias_sd",
      "p05_bias_sd", "p95_bias_sd", "Standard-deviation bias"
    ) +
      labs(y = "SD bias (est - true)") +
      theme(
        strip.text = element_blank(), strip.background = element_blank(),
        axis.text = element_text(size = 17), axis.text.x = element_text(size = 15),
        axis.title.y = element_text(size = 18), plot.title = element_text(size = 24)
      )
    convergence_panel <- plot_convergence_rate(level_runs) +
      theme(
        strip.text = element_blank(), strip.background = element_blank(),
        axis.text = element_text(size = 17), axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 18), plot.title = element_text(size = 24)
      )

    level_body <- cowplot::plot_grid(
      mean_panel, sd_panel, convergence_panel,
      nrow = 1, rel_widths = c(1, 1, 0.9), align = "h", axis = "tb"
    )
    level_header <- cowplot::ggdraw() +
      cowplot::draw_label(
        sprintf("%s - sample size N = %s", level, sample_size),
        x = 0.01, y = 0.72, hjust = 0, fontface = "bold",
        size = 28, colour = "#1F2937"
      ) +
      cowplot::draw_label(
        sprintf(
          "Bias: converged runs only | Convergence: all %s repetitions | %s-iteration limit",
          repetitions, iteration_limit
        ),
        x = 0.01, y = 0.20, hjust = 0, size = 17, colour = "#5E6874"
      )
    level_plot <- cowplot::plot_grid(
      level_header, level_body, ncol = 1, rel_heights = c(0.19, 1)
    )
    level_path <- file.path(
      combined_output_dir,
      paste0(base_name, "_", dispersion_files[[level]], "_mean_sd_convergence.jpg")
    )
    ggsave(
      level_path, level_plot, width = 18, height = 6.6, dpi = 300,
      bg = "white", device = grDevices::jpeg, type = "cairo",
      limitsize = FALSE
    )
    dispersion_paths <- c(dispersion_paths, level_path)
  }

  if (generate_tables) {
    scenario_results <- read.xlsx(input_path, sheet = "Per_Scenario") %>%
      mutate(
        profile = recode(
          weight_name,
          baseline = "Baseline",
          early_dominant = "CP-heavy",
          heavy_tailed = "Higher-generation-heavy"
        )
      )
    sheet_list <- list(
      "Mean2_SD2" = scenario_results %>%
        filter(true_mean == 2, true_sd == 2) %>%
        format_for_excel(),
      "Mean15_SD3" = scenario_results %>%
        filter(true_mean == 15, true_sd == 3) %>%
        format_for_excel(),
      "Mean5_SD10" = scenario_results %>%
        filter(true_mean == 5, true_sd == 10) %>%
        format_for_excel()
    )
    table_path <- file.path(output_dir, paste0(base_name, "_Formatted_Tables.xlsx"))
    write.xlsx(sheet_list, file = table_path, rowNames = FALSE)
  }

  if (save_combined) message(sprintf("Combined summary saved to: %s", combined_path))
  message(sprintf("Dispersion summaries saved to: %s", paste(dispersion_paths, collapse = ", ")))
  invisible(list(
    mean_plot = p_mean,
    sd_plot = p_sd,
    convergence_plot = p_convergence,
    combined_plot = combined_plot,
    combined_path = combined_path,
    dispersion_paths = dispersion_paths
  ))
}
