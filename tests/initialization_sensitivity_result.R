
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# ============================================================
# 1. Input
# ============================================================

input_file <- file.path(
  "notes",
  "Simulation_results",
  "Init_sensitivity",
  "Init_sensitivity_N100_R100_I300.xlsx"
)

out_dir <- dirname(input_file)

true_mean <- 15
true_sd <- 3
iteration_cap <- 300

# Display limits
sd_display <- c(-1.5, 14)
iter_display <- c(0, 200)

method_keys <- c(
  "true_value",
  "mean_misspecified",
  "sd_misspecified",
  "icc_based"
)

method_names <- c(
  "True-parameter start",
  "Mean misspecified",
  "Standard deviation misspecified",
  "ICC-based"
)

profile_keys <- c(
  "baseline",
  "early_dominant",
  "heavy_tailed"
)

profile_names <- c(
  "Baseline",
  "CP-heavy",
  "Higher-generation-heavy"
)

method_colors <- c(
  "True-parameter start" = "#15796E",
  "Mean misspecified" = "#5D73B5",
  "Standard deviation misspecified" = "#D17F3F",
  "ICC-based" = "#9A5B92"
)

# ============================================================
# 2. Read data
# ============================================================

# A temporary ASCII-path copy avoids Excel zip-reading failures caused by
# non-ASCII project paths on some Windows R installations.
input_copy <- tempfile(fileext = ".xlsx")
stopifnot(file.copy(input_file, input_copy, overwrite = TRUE))

dat <- read_excel(
  input_copy,
  sheet = "Per_Repetition"
) %>%
  filter(
    profile %in% profile_keys,
    init_method %in% method_keys
  ) %>%
  mutate(
    profile = factor(
      profile,
      levels = profile_keys,
      labels = profile_names
    ),
    init_method = factor(
      init_method,
      levels = method_keys,
      labels = method_names
    ),
    est_mean = as.numeric(est_mean),
    est_sd = as.numeric(est_sd),
    iterations = as.numeric(iterations),
    converged = as.logical(converged),
    bias_mean = est_mean - true_mean,
    bias_sd = est_sd - true_sd
  )

# ============================================================
# 3. Summary statistics
# ============================================================

percentile_boxes <- function(df, variable) {
  df %>%
    group_by(profile, init_method) %>%
    summarise(
      ymin = quantile(
        .data[[variable]], 0.05, na.rm = TRUE
      ),
      lower = quantile(
        .data[[variable]], 0.25, na.rm = TRUE
      ),
      middle = median(
        .data[[variable]], na.rm = TRUE
      ),
      upper = quantile(
        .data[[variable]], 0.75, na.rm = TRUE
      ),
      ymax = quantile(
        .data[[variable]], 0.95, na.rm = TRUE
      ),
      .groups = "drop"
    )
}

mean_boxes <- percentile_boxes(dat, "bias_mean")
sd_boxes <- percentile_boxes(dat, "bias_sd")

iteration_summary <- dat %>%
  group_by(profile, init_method) %>%
  summarise(
    median_iter = median(
      iterations, na.rm = TRUE
    ),
    q25_iter = quantile(
      iterations, 0.25, na.rm = TRUE
    ),
    q75_iter = quantile(
      iterations, 0.75, na.rm = TRUE
    ),
    nonconverged = sum(
      !converged, na.rm = TRUE
    ),
    .groups = "drop"
  )

# Count observations hidden by SD-axis zoom.
# They are NOT removed from the statistics.

sd_outside <- dat %>%
  filter(
    !is.na(bias_sd),
    bias_sd < sd_display[1] |
      bias_sd > sd_display[2]
  )

if (nrow(sd_outside) > 0) {
  message(
    "SD observations outside plot range: ",
    nrow(sd_outside)
  )

  print(
    sd_outside %>%
      select(profile, rep, init_method, bias_sd)
  )
}

make_limits <- function(x) {
  r <- range(c(0, x), na.rm = TRUE)
  pad <- max(0.08 * diff(r), 0.5)
  c(r[1] - pad, r[2] + pad)
}

mean_display <- make_limits(dat$bias_mean)

# ============================================================
# 4. Theme
# ============================================================

base_theme <- theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 17,
      hjust = 0.5
    ),

    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    panel.grid.major.y = element_line(
      color = "#E5E9F0",
      linewidth = 0.4
    ),

    panel.background = element_rect(
      fill = "#F9FAFC",
      color = NA
    ),

    # Remove ALL repeated bottom x-axis labels.
    # A single four-color legend will be used.
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),

    axis.text.y = element_text(
      size = 12.5,
      color = "#364152"
    ),

    axis.title.y = element_text(size = 14),

    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 5. Mean / SD bias panels
# ============================================================

make_bias_panel <- function(
    pr,
    metric,
    boxes,
    title = NULL,
    y_label = NULL,
    y_window = NULL,
    show_legend = FALSE
) {

  sub <- filter(dat, profile == pr)
  sums <- filter(boxes, profile == pr)

  ggplot(
    sub,
    aes(
      x = init_method,
      y = .data[[metric]]
    )
  ) +

    geom_hline(
      yintercept = 0,
      color = "#D62728",
      linewidth = 0.75,
      linetype = "dashed"
    ) +

    geom_jitter(
      aes(color = init_method),
      width = 0.17,
      height = 0,
      alpha = 0.28,
      size = 0.9,
      show.legend = FALSE
    ) +

    geom_boxplot(
      data = sums,

      aes(
        x = init_method,
        ymin = ymin,
        lower = lower,
        middle = middle,
        upper = upper,
        ymax = ymax,
        fill = init_method,
        color = init_method
      ),

      stat = "identity",
      inherit.aes = FALSE,
      width = 0.55,
      alpha = 0.20,
      linewidth = 0.65,
      show.legend = show_legend
    ) +

    scale_color_manual(
      values = method_colors,
      guide = "none"
    ) +

    scale_fill_manual(
      values = method_colors,
      guide = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(alpha = 0.85)
      )
    ) +

    # Zoom without deleting observations or
    # changing percentile calculations.
    coord_cartesian(
      ylim = y_window
    ) +

    labs(
      title = title,
      x = NULL,
      y = y_label,
      fill = NULL
    ) +

    base_theme
}

# ============================================================
# 6. EM iterations: median bars + IQR
# ============================================================

make_iteration_panel <- function(
    pr,
    title = NULL
) {

  sub <- filter(
    iteration_summary,
    profile == pr
  )

  ggplot(
    sub,
    aes(
      x = init_method,
      y = median_iter,
      fill = init_method
    )
  ) +

    geom_col(
      width = 0.63,
      alpha = 0.78,
      show.legend = FALSE
    ) +

    geom_errorbar(
      aes(
        ymin = q25_iter,
        ymax = q75_iter
      ),
      width = 0.20,
      linewidth = 0.8,
      color = "#253041",
      show.legend = FALSE
    ) +

    # Count runs that reached the 300-iteration cap.
    geom_text(
      data = filter(
        sub,
        nonconverged > 0
      ),

      aes(
        y = pmax(
          q75_iter,
          median_iter
        ) + 10,

        label = paste0(
          "cap: ",
          nonconverged
        )
      ),

      color = "#B93737",
      size = 4.2,
      show.legend = FALSE
    ) +

    scale_fill_manual(
      values = method_colors,
      guide = "none"
    ) +

    scale_y_continuous(
      breaks = seq(0, 200, 50),
      expand = expansion(
        mult = c(0, 0.02)
      )
    ) +

    # Display upper limit is 200, NOT 300.
    coord_cartesian(
      ylim = iter_display
    ) +

    labs(
      title = title,
      x = NULL,
      y = "EM iterations"
    ) +

    base_theme
}

# ============================================================
# 7. Assemble nine panels
# ============================================================

rows <- lapply(
  seq_along(profile_names),
  function(i) {

    pr <- profile_names[i]
    first_row <- i == 1

    # Column 1: Mean Bias
    p_mean <- make_bias_panel(
      pr,
      metric = "bias_mean",
      boxes = mean_boxes,

      title = if (first_row) {
        "SI Mean Bias"
      } else {
        NULL
      },

      y_label = paste0(
        pr,
        "\n\nEstimation bias"
      ),

      y_window = mean_display,

      # Generate the legend only ONCE.
      show_legend = first_row
    )

    # Column 2: SD Bias
    p_sd <- make_bias_panel(
      pr,
      metric = "bias_sd",
      boxes = sd_boxes,

      title = if (first_row) {
        "SI SD Bias"
      } else {
        NULL
      },

      y_label = "Estimation bias",
      y_window = sd_display,
      show_legend = FALSE
    )

    # Column 3: Iterations
    p_iter <- make_iteration_panel(
      pr,

      title = if (first_row) {
        "EM Iterations"
      } else {
        NULL
      }
    )

    p_mean + p_sd + p_iter +
      plot_layout(
        widths = c(1, 1, 1)
      )
  }
)

# ============================================================
# 8. Caption
# ============================================================

clip_caption <- if (nrow(sd_outside) > 0) {

  paste0(
    "SD-bias axis zoomed to [",
    sd_display[1], ", ", sd_display[2],
    "]; ",
    nrow(sd_outside),
    " extreme observation(s) outside the plotted ",
    "range (retained in summaries)."
  )

} else {

  paste0(
    "SD-bias axis shown within [",
    sd_display[1], ", ", sd_display[2],
    "]."
  )

}

# ============================================================
# 9. Final plot
# ============================================================

final_plot <- (
  rows[[1]] /
    rows[[2]] /
    rows[[3]]
) +

  plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) +

  plot_annotation(
    title = paste(
      "Initialization Sensitivity",
      "Under Low SI Dispersion"
    ),

    subtitle = paste0(
      "True SI mean = 15, SD = 3",
      "   |   N = 100",
      "   |   100 repetitions per pathway profile"
    ),

    caption = paste0(
      "SI bias = estimate - truth; ",
      "boxes = IQR; whiskers = 5th-95th percentiles; ",
      "dots = repetitions; red dashed line = zero bias.\n",

      "EM iterations: columns = median, ",
      "error bars = 25th-75th percentiles; ",
      "maximum EM iterations = ",
      iteration_cap,
      "; cap = non-converged runs.\n",

      clip_caption,

      " All fits, including iteration-capped runs, ",
      "are retained."
    ),

    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 24
      ),
      plot.subtitle = element_text(
        size = 14,
        color = "#667085"
      ),
      plot.caption = element_text(
        size = 10.5,
        color = "#667085",
        hjust = 0
      ),
      plot.margin = margin(
        14, 18, 14, 18
      )
    )
  ) &

  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_blank(),

    legend.text = element_text(
      size = 13
    ),

    legend.spacing.x = grid::unit(
      13, "pt"
    )
  )

print(final_plot)

# ============================================================
# 10. Export
# ============================================================

ggsave(
  filename = file.path(
    out_dir,
    "Init_sensitivity_3x3_Bias_optimized.png"
  ),
  plot = final_plot,
  width = 18,
  height = 13.5,
  dpi = 300,
  bg = "white",
  device = grDevices::png,
  type = "cairo"
)

ggsave(
  filename = file.path(
    out_dir,
    "Init_sensitivity_3x3_Bias_optimized.pdf"
  ),
  plot = final_plot,
  width = 18,
  height = 13.5,
  bg = "white"
)

message("Figure files saved to: ", out_dir)
