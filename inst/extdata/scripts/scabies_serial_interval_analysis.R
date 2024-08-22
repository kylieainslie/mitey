# Determining serial interval of scabies using epidemic curves from different
# sources

# load required packages
library(tidyverse)
library(brms)
library(devtools)
load_all()

# we will calculate the index case-to-case (ICC) interval for each person by class
# the person with the greatest value for number of days since symptom onset will
# be considered the index case. The rest of the class members will have an ICC
# interval calculated as the number of days between their symptom onset and the
# index case.

# read in data
si_data <- readRDS("inst/extdata/data/si_data.rds")

# use method from Vink et al. to estimate SI for each study
# assume a normal distribution, then do some wrangling
result_norm <- si_data %>%
  select(icc_interval, study) %>%
  group_by(study) %>%
  summarise(result = list(si_estim(icc_interval))) %>%
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


# perform sensitivity analysis
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
# Data preparation
x <- seq(0, 250, by = 1)
y1 <- dnorm(x, mean = 167.34442, sd = 9.71763)
y2 <- dnorm(x, mean = 98.4, sd = 8.542332)
y3 <- dnorm(x, mean = 122.92385, sd = 26.92035)
y4 <- dnorm(x, mean = 21.91776, sd = 15.23666)
y5 <- dnorm(x, mean = 110.71571, sd = 16.13879)
y6 <- dnorm(x, mean = 16.106246, sd = 2.421762)
# Create a data frame
my_data <- data.frame(
  x = rep(x, 6),
  y = c(y1, y2, y3, y4, y5, y6),
  group = factor(rep(c("Ghana outbreak", "Germany outbreak", "Uganda outbreak",
                       "Nevada (US) outbreak", "Netherlands outbreak", "Spain outbreak"), each = length(x)),
                 levels = c("Ghana outbreak", "Germany outbreak", "Uganda outbreak",
                            "Nevada (US) outbreak", "Netherlands outbreak", "Spain outbreak"))
)

# Create the plot
ggplot(my_data, aes(x = x, y = y, color = group)) +
  geom_line(linewidth = 1) +
  labs(title = "Estimated Serial Interval Distributions",
       x = "Time (days)",
       y = "Density") +
  scale_color_viridis_d() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )

# ------------------------------------------------------------------------------

# SI Meta-analysis

# create data frame of effect sizes
df_effect_sizes <- data.frame(
  mean_si = c(167.34442,
              98.4,
              122.92385,
              #21.91776,
              110.71571#,
              #16.106246
              ),
  sd_si = c(9.71763,
            8.542332,
            26.92035,
            #15.23666,
            16.13879#,
            #2.421762
            ),
  n = c(nrow(ghana_df),
        nrow(germany_df),
        nrow(uganda_long_df),
        #nrow(nevada_long_df),
        nrow(dutch_long_df)#,
        #nrow(spain_long_df)
        ),
  country = c("Ghana",
              "Germany",
              "Uganda",
              #"Nevada (USA)",
              "Netherlands"#,
              #"Spain"
              ),
  article = c("Kaburi et al.",
              "Ariza et al.",
              "Akunzirwe et al.",
              #"Division of Public and Behavioral Health",
              "Tjon-Kon-Fat et al."#,
              #"Larrosa et al."
              )
) %>%
  mutate(se_si = sd_si/sqrt(n))

# we will perform a Bayesian meta-analysis using the {brms} package
# specify priors
priors <- c(prior(normal(100,50), class = Intercept),
            prior(cauchy(0,1), class = sd))

# fit a random effects model
# Fit the random effects model with adjusted control parameters
m.brm <- brm(
  mean_si | se(se_si) ~ 1 + (1 | article),
  data = df_effect_sizes,
  prior = priors,
  iter = 8000,  # Increased number of iterations
  warmup = 4000,  # Increased warmup
  control = list(adapt_delta = 0.999, max_treedepth = 20)  # Increased adapt_delta and max_treedepth
)

# check convergence
# we want to verify that Rhat = 1 (signifying convergence)
summary(m.brm)
# pooled mean SI = 95.57; 95% CI = (65.27, 124.83); sd = 15.17
# sd(Intercept) = 58.31, represents the standard deviation of the random
# intercepts (the between-study variability). This large estimate indicates that
# there is substantial heterogeneity among the studies.

# posterior predictive check
pp_check(m.brm)

# Look at posterior distributions
post.samples <- posterior_samples(m.brm, c("^b", "^sd"))
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
library(tidybayes)
library(dplyr)
library(ggplot2)
library(ggridges)
library(glue)
library(stringr)
library(forcats)

# get posterior draws from each study
study.draws <- spread_draws(m.brm, r_country[country,], b_Intercept) %>%
  mutate(b_Intercept = r_country + b_Intercept)

# get pooled posterior draws
pooled.effect.draws <- spread_draws(m.brm, b_Intercept) %>%
  mutate(country = "Pooled Effect")

# combine posterior draws from each study and pooled
forest.data <- bind_rows(study.draws,
                         pooled.effect.draws) %>%
  ungroup() %>%
  mutate(country = str_replace_all(country, "[.]", " ")) %>%
  mutate(country = reorder(country, b_Intercept))

# calculate mean and credible intervals
forest.data.summary <- group_by(forest.data, country) %>%
  mean_qi(b_Intercept)

# plot
ggplot(aes(b_Intercept,
           relevel(country, "Pooled Effect",
                   after = Inf)),
       data = forest.data) +

  # Add vertical lines for pooled effect and CI
  geom_vline(xintercept = fixef(m.brm)[1, 1],
             color = "grey", size = 1) +
  geom_vline(xintercept = fixef(m.brm)[1, 3:4],
             color = "grey", linetype = 2) +
  geom_vline(xintercept = 0, color = "black",
             size = 1) +

  # Add densities
  geom_density_ridges(fill = "blue",
                      rel_min_height = 0.01,
                      col = NA, scale = 1,
                      alpha = 0.8) +
  geom_pointintervalh(data = forest.data.summary,
                      size = 1) +

  # Add text and labels
  geom_text(data = mutate_if(forest.data.summary,
                             is.numeric, round, 2),
            aes(label = glue("{b_Intercept} [{.lower}, {.upper}]"),
                x = Inf), hjust = "inward") +
  labs(x = "Standardized Mean Difference", # summary measure
       y = element_blank()) +
  theme_minimal()

# ---------
# Create a forest plot with each study's estimate of mean SI and the pooled estimate
# Create the data frame with the study estimates and the pooled estimate
forest_data <- df_effect_sizes %>%
  mutate(
    lower = mean_si - 1.96 * se_si,
    upper = mean_si + 1.96 * se_si,
    type = "Study"
  ) %>%
  select(country, mean_si, lower, upper, type) %>%
  rename(estimate = mean_si) %>%
  bind_rows(
    data.frame(
      country = "Pooled Estimate",
      estimate = posterior_summary(m.brm, variable = "Intercept")[, "Estimate"],
      lower = posterior_summary(m.brm, variable = "Intercept")[, "Q2.5"],
      upper = posterior_summary(m.brm, variable = "Intercept")[, "Q97.5"],
      type = "Pooled"
    )
  ) %>%
  # Ensure the pooled estimate is positioned last
  mutate(country = factor(country, levels = c(setdiff(country, "Pooled Estimate"), "Pooled Estimate")))

saveRDS(forest_data, "vignettes/data_for_forest_plot.RDS")
# Create the forest plot using ggplot2
ggplot(forest_data, aes(x = estimate, y = country, xmin = lower, xmax = upper, color = type)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = posterior_summary(m.brm, variable = "Intercept")[, "Estimate"], linetype = "dashed", color = "red") +
  labs(x = "Effect Size", y = "") +
  theme_minimal() +
  scale_color_manual(values = c("Study" = "blue", "Pooled" = "red")) +
  theme(legend.position = "none")
