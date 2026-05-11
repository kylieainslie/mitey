library(dplyr)
library(e1071)  # pour skewness

stats <- validation_data %>%
  group_by(Pathogen, Country, Author) %>%
  summarise(
    mean = mean(ICC_interval, na.rm = TRUE),
    sd = sd(ICC_interval, na.rm = TRUE),
    cv = sd / mean,
    skewness = skewness(ICC_interval, na.rm = TRUE),
    median = median(ICC_interval, na.rm = TRUE),
    max = max(ICC_interval, na.rm = TRUE),
    max_median_ratio = max / median,
    .groups = "drop"
  )

print(stats)
