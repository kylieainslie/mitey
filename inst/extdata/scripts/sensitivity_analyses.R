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
  plotlist = labeled_plots_gam[-c(3,5)],
  labels = "AUTO",      # Automatically adds labels (A, B, C, etc.)
  label_size = 14,      # Size of the labels
  ncol = 2              # Number of columns; adjust as needed
)

# Display the final combined plot
print(final_plot_gam)

# save plot
ggsave("vignettes/figures/gam_epidemic_curve_density_plot.png", plot = final_plot_gam,
       width = 7, height = 5, dpi = 300, units = "in")

# ------------------------------------------------------------------------------
# Perform a Bayesian meta-analysis
df_ma_gam <- df_merged_gam %>%
  group_by(study) %>%
  mutate(n = n()) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(se = sd/sqrt(n)) %>%
  select(study, n, mean, sd, se) %>%
  # don't include Larrosa and DPBH because their study periods are too short
  filter(study %in% c("Akunzirwe et al.",
                      "Ariza et al.",
                      "Kaburi et al.",
                      "Tjon-Kon-Fat et al"))

# we will perform a Bayesian meta-analysis using the {brms} package
# specify priors
priors <- c(prior(normal(100,50), class = Intercept),
            prior(cauchy(0,1), class = sd))

# fit a random effects model
# Fit the random effects model with adjusted control parameters
m.brm.g <- brm(
  mean | se(se) ~ 1 + (1 | study),
  data = df_ma_gam,
  prior = priors,
  iter = 8000,  # Increased number of iterations
  warmup = 4000,  # Increased warmup
  control = list(adapt_delta = 0.999, max_treedepth = 20)  # Increased adapt_delta and max_treedepth
)

# check convergence
# we want to verify that Rhat = 1 (signifying convergence)
summary(m.brm.g)
# Multilevel Hyperparameters:
#   ~study (Number of levels: 4)
# Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    34.40     17.55    14.82    79.67 1.00     4335     5394
#
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept   117.94     18.80    76.45   151.84 1.00     4487     4554
# posterior predictive check
pp_check(m.brm.g)

# Look at posterior distributions
post.samples <- posterior_samples(m.brm.g, c("^b", "^sd"))
names(post.samples)
names(post.samples) <- c("smd", "tau")

ggplot(aes(x = smd), data = post.samples) +
  geom_density(fill = "lightblue",                # set the color
               color = "lightblue", alpha = 0.7) +
  geom_point(y = 0,                               # add point at mean
             x = mean(post.samples$smd)) +
  labs(x = expression(italic(SMD)),
       y = element_blank()) +
  theme_minimal()

ggplot(aes(x = tau), data = post.samples) +
  geom_density(fill = "lightgreen",               # set the color
               color = "lightgreen", alpha = 0.7) +
  geom_point(y = 0,
             x = mean(post.samples$tau)) +        # add point at mean
  labs(x = expression(tau),
       y = element_blank()) +
  theme_minimal()

# Create forest plot with posteriors

# get posterior draws from each study
study.draws <- spread_draws(m.brm.g, r_study[study,], b_Intercept) %>%
  mutate(b_Intercept = r_study + b_Intercept)

# get pooled posterior draws
pooled.effect.draws <- spread_draws(m.brm.g, b_Intercept) %>%
  mutate(study = "Pooled Effect")

# combine posterior draws from each study and pooled
forest.data <- bind_rows(study.draws,
                         pooled.effect.draws) %>%
  ungroup() %>%
  mutate(study = str_replace_all(study, "[.]", " ")) %>%
  mutate(study = reorder(study, b_Intercept))

# calculate mean and credible intervals
forest.data.summary <- group_by(forest.data, study) %>%
  mean_qi(b_Intercept)

# plot
forest_plot <- ggplot(aes(b_Intercept,
                          relevel(study, "Pooled Effect",
                                  after = Inf)),
                      data = forest.data) +

  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(m.brm.g)[1, 1],
             color = "grey", size = 1) +
  geom_vline(xintercept = fixef(m.brm.g)[1, 3:4],
             color = "grey", linetype = 2) +
  #geom_vline(xintercept = 0, color = "black",
  #           size = 1) +

  # Add densities
  geom_density_ridges(fill = "blue",
                      rel_min_height = 0.01,
                      col = NA, scale = 1,
                      alpha = 0.8) +

  geom_pointinterval(aes(y = study,
                         x = b_Intercept,
                         xmin = .lower,
                         xmax = .upper),
                     data = forest.data.summary,
                     size = 1,
                     orientation = "horizontal") +

  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary,
                             is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),
                x = Inf), hjust = "inward") +
  labs(x = "Mean Serial Interval (days)", # summary measure
       y = element_blank()) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    axis.line = element_line(color = "black", linewidth = 0.5), # Add axis lines
    axis.ticks = element_line(color = "black"), # Add axis ticks
    axis.title = element_text(size = 12, face = "bold"), # Customize axis titles
    axis.text = element_text(size = 10) # Customize axis text
  )

print(forest_plot)

# save plot
ggsave("vignettes/figures/gam_meta-analysis_forest_plot.png", plot = forest_plot,
       width = 7, height = 5, dpi = 300, units = "in")

