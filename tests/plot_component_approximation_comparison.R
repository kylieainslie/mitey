
# ============================================================
# Component-Density Approximation Assessment
#
# CP: Exact |X1-X2| vs Half-Normal approximation
# PS: Exact Lognormal (no approximation)
# PT: Exact sum of 2 Lognormals vs Fenton-Wilkinson
# PQ: Exact sum of 3 Lognormals vs Fenton-Wilkinson
#
# Three scenarios:
# Low      (mean = 15, SD = 3)
# Moderate (mean = 2,  SD = 2)
# High     (mean = 5,  SD = 10)
#
# Output:
#   component_approximation_comparison_corrected.png
#   component_approximation_comparison_corrected.pdf
#   component_approximation_overlap_corrected.csv
# ============================================================

library(ggplot2)


# ------------------------------------------------------------
# 1. Find mitey package root
# ------------------------------------------------------------

find_pkg_root <- function() {

  file_arg <- grep(
    "^--file=",
    commandArgs(FALSE),
    value = TRUE
  )

  script_dir <- if (length(file_arg)) {
    dirname(sub("^--file=", "", file_arg[1]))
  } else {
    NULL
  }

  starts <- unique(c(getwd(), script_dir))

  for (start in starts) {

    current <- normalizePath(
      start,
      winslash = "/",
      mustWork = FALSE
    )

    for (i in 0:6) {

      if (file.exists(file.path(current, "DESCRIPTION"))) {
        return(current)
      }

      parent <- dirname(current)

      if (identical(parent, current)) break

      current <- parent
    }
  }

  stop(
    "Could not find mitey project root. ",
    "Please set the working directory to the mitey project."
  )
}


# ------------------------------------------------------------
# 2. Lognormal parameter conversion
# ------------------------------------------------------------

# Convert natural-scale mean and SD into meanlog/sdlog.

lnorm_parameters <- function(mean, sd) {

  sdlog <- sqrt(
    log1p((sd / mean)^2)
  )

  meanlog <- log(mean) - 0.5 * sdlog^2

  c(
    meanlog = meanlog,
    sdlog = sdlog
  )
}


# ------------------------------------------------------------
# 3. FFT convolution for PT and PQ
# ------------------------------------------------------------

# Retained from your original script.

fft_convolution_power <- function(probability_mass, power) {

  output_length <- power *
    (length(probability_mass) - 1) + 1

  fft_length <- 2^ceiling(
    log2(output_length)
  )

  transformed <- fft(
    c(
      probability_mass,
      rep(
        0,
        fft_length - length(probability_mass)
      )
    )
  )

  result <- Re(
    fft(
      transformed^power,
      inverse = TRUE
    )
  ) / fft_length

  pmax(
    result[seq_len(output_length)],
    0
  )
}


# ------------------------------------------------------------
# 4. NEW: Exact CP density using FFT autocorrelation
# ------------------------------------------------------------

# For CP:
#
#   C_CP = |X1 - X2|
#
# where X1, X2 are iid Lognormal.
#
# Exact continuous density:
#
#   f_CP(c) = 2 * integral f_X(x) f_X(x+c) dx
#
# This is calculated numerically by autocorrelation
# of discretised Lognormal probability masses.
#
# Unlike the previous script, CP no longer compares
# the Half-Normal approximation with itself.

fft_absolute_difference_density <- function(
    probability_mass,
    dx
) {

  L <- length(probability_mass)

  fft_length <- 2^ceiling(
    log2(2 * L - 1)
  )

  padded <- c(
    probability_mass,
    rep(0, fft_length - L)
  )

  transformed <- fft(padded)

  # Autocorrelation
  correlation <- Re(
    fft(
      Conj(transformed) * transformed,
      inverse = TRUE
    )
  ) / fft_length

  # Fold positive and negative differences
  # to obtain the density of |X1-X2|.
  #
  # The factor 2 applies to the density at c >= 0.

  density <- 2 *
    correlation[seq_len(L)] / dx

  pmax(density, 0)
}


# ------------------------------------------------------------
# 5. Numerical integration helper
# ------------------------------------------------------------

trapz <- function(x, y) {

  sum(
    diff(x) *
      (
        head(y, -1) +
          tail(y, -1)
      ) / 2
  )
}


# ------------------------------------------------------------
# 6. Compare exact and approximating densities
# ------------------------------------------------------------

compare_component <- function(
    mean,
    sd,
    path,
    plot_probability = 0.995,
    integration_probability = 0.999999,
    bins = 32768,
    plot_points = 1800
) {

  pars <- lnorm_parameters(mean, sd)

  meanlog <- pars["meanlog"]
  sdlog <- pars["sdlog"]

  # ----------------------------------------------------------
  # CP: Exact absolute difference vs Half-Normal
  # ----------------------------------------------------------

  if (path == "CP") {

    # Truncate the underlying Lognormal distribution
    # sufficiently far into its right tail.

    base_probability <- 1 -
      (1 - integration_probability) / 3

    base_max <- qlnorm(
      base_probability,
      meanlog = meanlog,
      sdlog = sdlog
    )

    edges <- seq(
      0,
      base_max,
      length.out = bins + 1
    )

    dx <- edges[2] - edges[1]

    # Probability in each discretisation bin
    base_mass <- diff(
      plnorm(
        edges,
        meanlog = meanlog,
        sdlog = sdlog
      )
    )

    x_full <- (
      seq_along(base_mass) - 1
    ) * dx

    # Exact CP density
    exact_full <- fft_absolute_difference_density(
      base_mass,
      dx
    )

    # Normalize small truncation/numerical errors
    exact_full <- exact_full /
      trapz(x_full, exact_full)

    # Half-Normal working approximation
    #
    # f(c) = exp[-c^2 / (4*s^2)] / (s*sqrt(pi))

    approximation_full <- exp(
      -x_full^2 / (4 * sd^2)
    ) / (sd * sqrt(pi))

    approximation_full <- approximation_full /
      trapz(x_full, approximation_full)

    # Full-domain density overlap
    overlap <- trapz(
      x_full,
      pmin(
        exact_full,
        approximation_full
      )
    )

    # Choose a readable plotting range.
    # Overlap itself is calculated BEFORE this truncation.

    exact_cdf <- cumsum(exact_full) * dx

    exact_limit <- x_full[
      which(exact_cdf >= plot_probability)[1]
    ]

    halfnormal_limit <- qnorm(
      (1 + plot_probability) / 2,
      mean = 0,
      sd = sqrt(2) * sd
    )

    plot_limit <- max(
      exact_limit,
      halfnormal_limit
    )

    available <- which(
      x_full <= plot_limit
    )

    selected <- unique(
      round(
        seq(
          min(available),
          max(available),
          length.out = min(
            plot_points,
            length(available)
          )
        )
      )
    )

    x <- x_full[selected]

    exact <- exact_full[selected]

    approximation <- approximation_full[selected]

  # ----------------------------------------------------------
  # PS: Exact Lognormal
  # ----------------------------------------------------------

  } else if (path == "PS") {

    x_max <- qlnorm(
      plot_probability,
      meanlog = meanlog,
      sdlog = sdlog
    )

    x <- seq(
      0,
      x_max,
      length.out = plot_points
    )

    exact <- dlnorm(
      x,
      meanlog = meanlog,
      sdlog = sdlog
    )

    approximation <- exact

    overlap <- 1

  # ----------------------------------------------------------
  # PT/PQ: Exact convolution vs FW approximation
  # ----------------------------------------------------------

  } else {

    generations <- unname(
      c(PT = 2, PQ = 3)[path]
    )

    # Truncate each base SI sufficiently far
    # into its right tail.

    base_probability <- 1 -
      (1 - integration_probability) / generations

    base_max <- qlnorm(
      base_probability,
      meanlog = meanlog,
      sdlog = sdlog
    )

    edges <- seq(
      0,
      base_max,
      length.out = bins + 1
    )

    dx <- edges[2] - edges[1]

    base_mass <- diff(
      plnorm(
        edges,
        meanlog = meanlog,
        sdlog = sdlog
      )
    )

    # Exact sum using FFT convolution
    exact_mass <- fft_convolution_power(
      base_mass,
      generations
    )

    x_full <- (
      (seq_along(exact_mass) - 1) +
        generations / 2
    ) * dx

    exact_full <- exact_mass / dx

    # Fenton-Wilkinson approximation:
    # matches exact natural-scale mean and SD.

    approx_pars <- lnorm_parameters(
      mean = generations * mean,
      sd = sqrt(generations) * sd
    )

    approximation_full <- dlnorm(
      x_full,
      meanlog = approx_pars["meanlog"],
      sdlog = approx_pars["sdlog"]
    )

    # Normalize numerical and tail errors
    exact_full <- exact_full /
      (sum(exact_full) * dx)

    approximation_full <- approximation_full /
      (sum(approximation_full) * dx)

    # Full-domain overlap
    overlap <- sum(
      pmin(
        exact_full,
        approximation_full
      )
    ) * dx

    # Select plotting domain from the distributions

    exact_cdf <- cumsum(exact_full) * dx

    exact_limit <- x_full[
      which(exact_cdf >= plot_probability)[1]
    ]

    approx_limit <- qlnorm(
      plot_probability,
      meanlog = approx_pars["meanlog"],
      sdlog = approx_pars["sdlog"]
    )

    plot_limit <- max(
      exact_limit,
      approx_limit
    )

    available <- which(
      x_full <= plot_limit
    )

    selected <- unique(
      round(
        seq(
          min(available),
          max(available),
          length.out = min(
            plot_points,
            length(available)
          )
        )
      )
    )

    x <- c(0, x_full[selected])

    exact <- c(0, exact_full[selected])

    approximation <- c(
      0,
      approximation_full[selected]
    )
  }

  # ----------------------------------------------------------
  # Return comparison
  # ----------------------------------------------------------

  list(

    curves = data.frame(
      x = rep(x, 2),

      density = c(
        exact,
        approximation
      ),

      distribution = rep(
        c(
          "Ideal construction",
          "Working approximation"
        ),
        each = length(x)
      )
    ),

    difference = data.frame(
      x = x,
      lower = pmin(
        exact,
        approximation
      ),
      upper = pmax(
        exact,
        approximation
      )
    ),

    overlap = min(
      max(overlap, 0),
      1
    )
  )
}


# ------------------------------------------------------------
# 7. Simulation scenarios
# ------------------------------------------------------------

scenarios <- data.frame(
  scenario = factor(
    c(
      "Low dispersion",
      "Moderate dispersion",
      "High dispersion"
    ),
    levels = c(
      "Low dispersion",
      "Moderate dispersion",
      "High dispersion"
    )
  ),

  mean = c(15, 2, 5),
  sd = c(3, 2, 10)
)

paths <- data.frame(
  path = c(
    "CP",
    "PS",
    "PT",
    "PQ"
  ),

  path_name = c(
    "Co-primary",
    "Primary-secondary",
    "Primary-tertiary",
    "Primary-quaternary"
  )
)


# ------------------------------------------------------------
# 8. Build data for all 12 panels
# ------------------------------------------------------------

curve_list <- list()
difference_list <- list()
annotation_list <- list()
overlap_list <- list()

panel_levels <- character()

index <- 1

for (i in seq_len(nrow(scenarios))) {

  for (j in seq_len(nrow(paths))) {

    scenario <- scenarios[i, ]

    path <- paths$path[j]

    comparison <- compare_component(
      mean = scenario$mean,
      sd = scenario$sd,
      path = path
    )

    panel <- sprintf(
      "%s (m=%g, s=%g)\n%s - %s",
      as.character(scenario$scenario),
      scenario$mean,
      scenario$sd,
      path,
      paths$path_name[j]
    )

    panel_levels <- c(
      panel_levels,
      panel
    )

    curve_list[[index]] <- transform(
      comparison$curves,
      panel = panel
    )

    difference_list[[index]] <- transform(
      comparison$difference,
      panel = panel
    )

    annotation_list[[index]] <- data.frame(
      panel = panel,

      overlap_label = sprintf(
        "Overlap: %.2f%%",
        100 * comparison$overlap
      )
    )

    overlap_list[[index]] <- data.frame(
      scenario = as.character(scenario$scenario),
      mean = scenario$mean,
      sd = scenario$sd,
      path = path,
      overlap = comparison$overlap,
      overlap_percent = 100 * comparison$overlap
    )

    index <- index + 1
  }
}


# ------------------------------------------------------------
# 9. Combine panel data
# ------------------------------------------------------------

curve_data <- do.call(
  rbind,
  curve_list
)

difference_data <- do.call(
  rbind,
  difference_list
)

annotation_data <- do.call(
  rbind,
  annotation_list
)

overlap_data <- do.call(
  rbind,
  overlap_list
)

curve_data$panel <- factor(
  curve_data$panel,
  levels = panel_levels
)

difference_data$panel <- factor(
  difference_data$panel,
  levels = panel_levels
)

annotation_data$panel <- factor(
  annotation_data$panel,
  levels = panel_levels
)

curve_data$distribution <- factor(
  curve_data$distribution,
  levels = c(
    "Ideal construction",
    "Working approximation"
  )
)


# ------------------------------------------------------------
# 10. Plotting
# ------------------------------------------------------------

comparison_colors <- c(
  "Ideal construction" = "#1F2937",
  "Working approximation" = "#D55E00"
)

comparison_linetypes <- c(
  "Ideal construction" = "solid",
  "Working approximation" = "dashed"
)

p <- ggplot() +

  # Shaded density disagreement
  geom_ribbon(
    data = difference_data,

    aes(
      x = x,
      ymin = lower,
      ymax = upper
    ),

    fill = "#CC6677",
    alpha = 0.25
  ) +

  # Exact and approximating densities
  geom_line(
    data = curve_data,

    aes(
      x = x,
      y = density,
      colour = distribution,
      linetype = distribution
    ),

    linewidth = 1.25
  ) +

  # Overlap percentage
  geom_text(
    data = annotation_data,

    aes(
      x = Inf,
      y = Inf,
      label = overlap_label
    ),

    hjust = 1.05,
    vjust = 1.4,

    colour = "#374151",
    size = 4,
    fontface = "bold"
  ) +

  # 3 rows x 4 columns
  facet_wrap(
    ~panel,
    ncol = 4,
    scales = "free"
  ) +

  scale_colour_manual(
    values = comparison_colors,
    breaks = names(comparison_colors)
  ) +

  scale_linetype_manual(
    values = comparison_linetypes,
    breaks = names(comparison_linetypes)
  ) +

  scale_x_continuous(
    expand = expansion(
      mult = c(0, 0.03)
    )
  ) +

  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +

  labs(
    title = "Accuracy of Component-Density Approximations",

    subtitle = paste(
      "Shading indicates density disagreement;",
      "overlap is the shared area under both curves"
    ),

    x = "ICC interval (days)",
    y = "Density",

    colour = NULL,
    linetype = NULL
  ) +

  # ----------------------------------------------------------
  # Larger fonts for Quarto PDF
  # ----------------------------------------------------------

  theme_minimal(base_size = 18) +

  theme(

    # Main title
    plot.title = element_text(
      face = "bold",
      size = 27,
      colour = "#111827",
      margin = margin(b = 10)
    ),

    # Subtitle
    plot.subtitle = element_text(
      size = 17,
      colour = "#4B5563",
      margin = margin(b = 16)
    ),

    # Panel headings
    strip.text = element_text(
      face = "bold",
      size = 14,
      lineheight = 1.08,
      colour = "#20242A",
      margin = margin(7, 5, 7, 5)
    ),

    strip.background = element_rect(
      fill = "#F3F4F6",
      colour = NA
    ),

    # Axis titles
    axis.title.x = element_text(
      size = 21,
      margin = margin(t = 13)
    ),

    axis.title.y = element_text(
      size = 21,
      margin = margin(r = 13)
    ),

    # Axis tick labels
    axis.text.x = element_text(
      size = 17,
      colour = "#374151"
    ),

    axis.text.y = element_text(
      size = 17,
      colour = "#374151"
    ),

    # Grid
    panel.grid.minor = element_blank(),

    panel.grid.major = element_line(
      colour = "#E5E7EB",
      linewidth = 0.40
    ),

    # Panel separation
    panel.spacing.x = grid::unit(0.8, "lines"),
    panel.spacing.y = grid::unit(1.1, "lines"),

    # Legend
    legend.position = "bottom",

    legend.direction = "horizontal",

    legend.text = element_text(
      size = 19
    ),

    legend.key.width = grid::unit(
      2.2,
      "lines"
    ),

    legend.spacing.x = grid::unit(
      15,
      "pt"
    ),

    # Outer margins
    plot.margin = margin(
      18,
      22,
      18,
      18
    )
  ) +

  guides(
    colour = guide_legend(
      nrow = 1,

      override.aes = list(
        linewidth = 1.8
      )
    )
  )


# ------------------------------------------------------------
# 11. Export to PNG, PDF, CSV
# ------------------------------------------------------------

pkg_root <- find_pkg_root()

output_dir <- file.path(
  pkg_root,
  "notes",
  "Simulation_results"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_png <- file.path(
  output_dir,
  "component_approximation_comparison_corrected.png"
)

output_pdf <- file.path(
  output_dir,
  "component_approximation_comparison_corrected.pdf"
)

output_csv <- file.path(
  output_dir,
  "component_approximation_overlap_corrected.csv"
)

# High-resolution raster image
ggsave(
  filename = output_png,
  plot = p,
  width = 16,
  height = 12,
  dpi = 320,
  bg = "white"
)

# Vector image for PDF reports
ggsave(
  filename = output_pdf,
  plot = p,
  width = 16,
  height = 12,
  device = grDevices::cairo_pdf,
  bg = "white"
)

# Save numerical overlap results
write.csv(
  overlap_data,
  file = output_csv,
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Print results
# ------------------------------------------------------------

print(p)

print(
  overlap_data[
    ,
    c(
      "scenario",
      "path",
      "overlap_percent"
    )
  ],
  row.names = FALSE
)

message("PNG saved: ", output_png)
message("PDF saved: ", output_pdf)
message("CSV saved: ", output_csv)
