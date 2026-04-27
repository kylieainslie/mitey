library(dplyr)
library(ggplot2)

# définir les colonnes qui identifient une étude
list_of_studies <- si_data2 %>%
  select(study) %>%
  distinct()

# boucle
for(i in 1:nrow(list_of_studies)) {

  study_list <- list_of_studies[i, ]

  subset_data <- si_data2 %>%
    filter(study == study_list$study
           )

  icc <- subset_data[[4]]  # colonne ICC

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

count(list_of_studies)
