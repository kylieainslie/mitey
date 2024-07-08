# Determining serial interval of scabies using epidemic curves from different
# sources

# load required packages
library(tidyverse)
library(ggplot2)
library(readxl)
library(devtools)
library(lubridate)
library(brms)
load_all()

# ------------------------------------------------------------------------------
# Kaburi et al. ----------------------------------------------------------------
# ------------------------------------------------------------------------------
# This paper describes an outbreak of scabies in a preschool in Ghana
# Kaburi et al. 2019 BMC Public Health https://doi.org/10.1186/s12889-019-7085-6

# read in epidemic curve data file
ghana_df <- read_xlsx("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies/data/Kaburi_et_al_data_scabies.xlsx")

# we will calculate the index case-to-case (ICC) interval for each person by class
# the person with the greatest value for number of days since symptom onset will
# be considered the index case. The rest of the class members will have an ICC
# interval calculated as the number of days between their symptom onset and the
# index case.

# icc_df <- kaburi_df %>%
#   # for now, we only need the class ID and the number of days since symptom onset
#   select(`Class (0=Creche; 1=Nursery 1; 2=Nursery 2; 3=KG1; 4=KG2)`,
#          `number of days since onset`) %>%
#   # rename variables for ease
#   rename(class = `Class (0=Creche; 1=Nursery 1; 2=Nursery 2; 3=KG1; 4=KG2)`,
#          no_days_since_onset = `number of days since onset`) %>%
#   # make an identification variable for whether or not an individual is an index case
#   group_by(class) %>%
#   mutate(index_case = if_else(no_days_since_onset == max(no_days_since_onset), 1, 0),
#          # calculate ICC intervals
#          icc_interval = abs(no_days_since_onset - max(no_days_since_onset))
#          )

# analysis without splitting observations into classes
icc_df2 <- ghana_df %>%
  # select the number of days since symptom onset
  select(`number of days since onset`) %>%
  # rename variables for ease
  rename(no_days_since_onset = `number of days since onset`) %>%
  # make an identification variable for whether or not an individual is an index case
  mutate(index_case = if_else(no_days_since_onset == max(no_days_since_onset), 1, 0),
  # calculate ICC intervals
         icc_interval = abs(no_days_since_onset - max(no_days_since_onset))
  )

# use method from Vink et al. to estimate SI
# si_estim(icc_df$icc_interval)
# [1] mean = 78.01803 sd = 67.18990   different results than script estimates
# [1] mean = 81.36419, sd = 62.89626

# without classes
si_estim(icc_df2$icc_interval)
# [1] mean = 165.95206  sd = 19.65646
# [1] mean = 167.34442  sd = 9.71763 script estimate
data <- icc_df2$icc_interval
source("./inst/extdata/scripts/SI_estimation_method_from_Vink_et_al.R")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Ariza et al. -----------------------------------------------------------------
# ------------------------------------------------------------------------------
# this data is from a pre-school in Germany
germany_df <- data.frame(
  id = c(seq(1, 16, by = 1)),
  date_onset = as.Date(c("01/12/2011", "25/02/2012", "01/03/2012", "02/03/2012",
                 "05/03/2012", "08/03/2012", "10/03/2012", "10/03/2012",
                 "13/03/2012", "13/03/2012", "14/03/2012", rep("16/03/2012",3),
                 "17/02/2012", "18/03/2012"), format = "%d/%m/%Y"),
  index_case = c(1, rep(0,15))
)

# calculate icc intervals from date of symptom onset
germany_df2 <- germany_df %>%
  mutate(icc_interval = as.integer(date_onset - min(date_onset)))

# use method from Vink et al. to estimate SI
si_estim(germany_df2$icc_interval)
# [1] mean = 98.4, sd = 8.542332

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Uganda outbreak --------------------------------------------------------------
# ------------------------------------------------------------------------------
# Akunzirwe et al 2023.
# https://uniph.go.ug/an-outbreak-of-scabies-in-a-fishing-community-in-hoima-district-uganda-february%E2%88%92june-2022/
# outbreak in a fishing community in 2022. Data are provided weekly, so all
# infections for each week are attributed to the first day of the week

uganda_df <- data.frame(
  date_onset = as.Date(c("01/01/2022", "08/01/2022", "15/01/2022", "29/01/2022",
                         "05/02/2022", "12/02/2022", "19/02/2022", "26/02/2022",
                         "05/03/2022", "12/03/2022", "19/03/2022", "26/03/2022",
                         "02/04/2022", "09/04/2022", "16/04/2022", "23/04/2022",
                         "30/04/2022", "07/05/2022", "14/05/2022", "21/05/2022",
                         "28/05/2022", "11/06/2022", "18/06/2022", "25/06/2022",
                         "02/07/2022"), format = "%d/%m/%Y"),
  num_cases = c(3, 4, 2, 5, 3, 2, 3, 6, 6, 4, 3, 2, 6, 20, 7, 4, 14, 8, 18, 2,
                10, 5, 10, 5, 3)
)

# Transform the data frame to long format
uganda_long_df <- uganda_df %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         icc_interval = as.integer(date_onset - min(date_onset)) # calculate ICC interval
  )

# use method from Vink et al. to estimate SI
si_estim(uganda_long_df$icc_interval)
# [1] mean = 122.92385  sd = 26.92035
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Nevada outbreak
# ------------------------------------------------------------------------------
# SCABIES OUTBREAK AMONG RESIDENTS OF A LONG TERM CARE FACILITY, CLARK COUNTY, NEVADA, 2015
# file:///Users/kylieainslie/Downloads/Scabies%20%E2%80%93%20Clark%20[v%202015%20i%2021%20e%201.0]_BP.pdf
nevada_df <- data.frame(
  date_onset = as.Date(c("07/04/2015", "08/04/2015", "16/04/2015", "17/04/2015",
                         "22/04/2015", "01/05/2015", "14/05/2015", "16/05/2015",
                         "18/05/2015"), format = "%d/%m/%Y"),
  num_cases = c(1, 1, 2, 1, 1, 1, 1, 1, 2)
)
# Transform the data frame to long format
nevada_long_df <- nevada_df %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         icc_interval = as.integer(date_onset - min(date_onset)) # calculate ICC interval
  )

# use method from Vink et al. to estimate SI
si_estim(nevada_long_df$icc_interval)
# [1] 21.91776 15.23666
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Dutch care home outbreak
# ------------------------------------------------------------------------------
# Tjon-Kon-Fat et al. (2021) Short report: The potential of PCR on skin flakes
# from bed linens for diagnosis of scabies in an outbreak.
# PLoS Negl Trop Dis 15(6): e0009485. https://doi.org/10.1371/journal.pntd.0009485
# cases identified by week, date of onset is attributed to first day of the week
#
dutch_df <- data.frame(
  week_onset = c(24, 28, 31, 33, 36, 37, 38, 40, 41, 42),
  num_cases = c(1, 1, 1, 2, 1, 2, 1, 1, 9, 2)
)

# Define the year
year <- 2018

# Convert week number to the first day of each week using piping
dutch_df <- dutch_df %>%
  mutate(date_onset = ymd(paste0(year, "-01-01")) + weeks(week_onset - 1)) %>%
  mutate(date_onset = floor_date(date_onset, unit = "week", week_start = 1))
# Convert to long format
dutch_long_df <- dutch_df %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         icc_interval = as.integer(date_onset - min(date_onset)) # calculate ICC interval
  )

# use method from Vink et al. to estimate SI
si_estim(dutch_long_df$icc_interval)
# [1] mean = 110.71571  sd = 16.13879
# ------------------------------------------------------------------------------
# img <- image_read("./inst/extdata/ethiopia_epidemic_curve.png")
# img_processed <- image_convert(img, "gray")  # Convert to grayscale
# tess <- tesseract(options = list(psm = 6))
# text <- ocr(img_processed, engine = tess)
# cat(text)

# ------------------------------------------------------------------------------
# Spain outbreak in a hospital
# Larrosa A., et al. Nosocomial outbreak of scabies in a hospital in Spain.
# Euro Surveill. 2003;8(10):pii=429. https://doi.org/10.2807/esm.08.10.00429-en
# ------------------------------------------------------------------------------
spain_df <- data.frame(
  day_onset = c(1, 2, 5, 9, 15, 16, 19, 21, 29, 44, 51, 61, 62),
  num_cases = c(1, 1, 2, 1, 1, 3, 1, 1, 1, 2, 1, 1, 1)
)

# Define start date
start_date <- as.Date("2002-11-05")

# Create sequence of dates based on day_onset
spain_df$date_onset <- start_date + (spain_df$day_onset - 1)

# Convert to long format and calculate icc_interval
spain_long_df <- spain_df %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         icc_interval = as.integer(date_onset - min(date_onset)) # calculate ICC interval
  )
# use method from Vink et al. to estimate SI
si_estim(spain_long_df$icc_interval)
# [1] 16.106246  2.421762
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
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
  mean_si = c(167.34442, 98.4, 122.92385, 21.91776, 110.71571, 16.106246),
  sd_si = c(9.71763, 8.542332, 26.92035, 15.23666, 16.13879, 2.421762),
  n = c(nrow(ghana_df), nrow(germany_df), nrow(uganda_long_df), nrow(nevada_long_df),
        nrow(dutch_long_df), nrow(spain_long_df)),
  country = c("Ghana", "Germany", "Uganda", "Nevada (USA)", "Netherlands", "Spain")
) %>%
  mutate(se_si = sd_si/sqrt(n))

# we will perform a Bayesian meta-analysis using the {brms} package
# specify priors
priors <- c(prior(normal(100,50), class = Intercept),
            prior(cauchy(0,1), class = sd))

# fit a random effects model
# convert sd to se
m.brm <- brm(mean_si|se(se_si) ~ 1 + (1|country),
             data = df_effect_sizes,
             prior = priors,
             iter = 5000,
             warmup = 2000,
             control = list(adapt_delta = 0.99, max_treedepth = 15))

# check convergence
# we want to verify that Rhat = 1 (signifying convergence)
summary(m.brm)
# pooled mean SI = 95.57; 95% CI = (65.27, 124.83); sd = 15.17
# sd(Intercept) = 58.31, represents the standard deviation of the random
# intercepts (the between-study variability). This large estimate indicates that
# there is substantial heterogeneity among the studies.

# posterior predictive check
pp_check(m.brm)

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

# Create the forest plot using ggplot2
ggplot(forest_data, aes(x = estimate, y = country, xmin = lower, xmax = upper, color = type)) +
  geom_point() +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = posterior_summary(m.brm, variable = "Intercept")[, "Estimate"], linetype = "dashed", color = "red") +
  labs(x = "Effect Size", y = "", title = "Forest Plot of Study Estimates and Pooled Estimate") +
  theme_minimal() +
  scale_color_manual(values = c("Study" = "blue", "Pooled" = "red")) +
  theme(legend.position = "none")
