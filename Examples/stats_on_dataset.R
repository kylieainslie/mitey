library(dplyr)
library(e1071)  # pour skewness

stats <- validation_data %>%
  group_by(Pathogen, Country, Author) %>%
  summarise(
    mean = mean(ICC_interval, na.rm = TRUE),
    sd = sd(ICC_interval, na.rm = TRUE),
    cv = sd / mean,
    median = median(ICC_interval, na.rm = TRUE),
    max = max(ICC_interval, na.rm = TRUE),
    .groups = "drop"
  )

print(stats)
