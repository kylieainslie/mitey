library(dplyr)
library(ggplot2)

# définir les colonnes qui identifient une étude
list_of_studies <- validation_data %>%
  select(Pathogen, Country, Author) %>%
  distinct()

# boucle
for(i in 1:nrow(list_of_studies)) {

  study <- list_of_studies[i, ]

  subset_data <- validation_data %>%
    filter(Pathogen == study$Pathogen,
           Country == study$Country,
           Author == study$Author)

  icc <- subset_data[[5]]  # colonne ICC

  df <- data.frame(value = icc)

  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(
      title = paste("ICC -", study$Pathogen, "|", study$Country, "|", study$Author),
      x = "ICC",
      y = "Fréquence"
    ) +
    theme_minimal()

  print(p)
}
