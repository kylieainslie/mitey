# Sensitivity Analyses

# we will calculate the index case-to-case (ICC) interval for each person by class
# the person with the greatest value for number of days since symptom onset will
# be considered the index case. The rest of the class members will have an ICC
# interval calculated as the number of days between their symptom onset and the
# index case.

# read in data
si_data <- readRDS("inst/extdata/data/si_data.rds")

# assume a gamma distribution
result_gam <- si_data %>%
  select(icc_interval, study) %>%
  group_by(study) %>%
  summarise(result = list(si_estim(icc_interval, dist = "gamma"))) %>%
  mutate(
    mean = map_dbl(result, "mean"),
    sd = map_dbl(result, "sd"),
    wts = map(result, "wts")  # Store wts as a list-column
  ) %>%
  select(-result) %>%
  unnest(wts) %>% # Unnest the wts column if needed %>%
  pivot_longer(
    cols = c(mean, sd, wts),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  group_by(study, statistic) %>%
  mutate(
    occurrence = row_number(),
    statistic = if_else(statistic == "wts", paste0("weight_", occurrence), statistic)
  ) %>%
  filter(statistic != "mean" | occurrence == 1) %>%
  filter(statistic != "sd" | occurrence == 1) %>%
  select(-occurrence) %>%
  ungroup()

# Plot serial interval curves
# Reshape results from long to wide format
result_gam_wide <- result_gam %>%
  pivot_wider(
    names_from = statistic,
    values_from = value
  )

# merge si_data and result_norm_wide for plotting
df_merged_gam <- si_data %>%
  select(study, icc_interval) %>%
  left_join(result_gam_wide, by = "study", relationship = "many-to-many")

# Apply the plot_si_fit function by study
plots_gam <- df_merged_gam %>%
  group_by(study) %>%
  group_map(~ plot_si_fit(
    dat = .x$icc_interval,
    mean = .x$mean[1],
    sd = .x$sd[1],
    weights = c(.x$weight_1[1], .x$weight_2[1], .x$weight_4[1]),
    dist = "gamma"
  ))

# Annotate plots with study names and labels
group_order_gam <- df_merged_gam %>%
  group_by(study) %>%
  group_keys()

labeled_plots_gam <- lapply(seq_along(plots_gam), function(i) {
  plots[[i]] +
    ggtitle(group_order_gam[i,1]) +            # Add study names as titles
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
})

# Combine plots into a multi-pane figure
final_plot_gam <- plot_grid(
  plotlist = labeled_plots_gam,
  labels = "AUTO",      # Automatically adds labels (A, B, C, etc.)
  label_size = 14,      # Size of the labels
  ncol = 2              # Number of columns; adjust as needed
)

# Display the final combined plot
print(final_plot_gam)
