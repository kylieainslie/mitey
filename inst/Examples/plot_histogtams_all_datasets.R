library(dplyr)
library(ggplot2)

validation_data <- readRDS("vignettes/articles/validation_data.rds")

# Créer un identifiant unique par étude
studies <- validation_data %>%
  distinct(Pathogen, Country, Author)

# Boucle sur toutes les études
plot_list <- vector("list", nrow(studies))

for (i in seq_len(nrow(studies))) {

  # Extraire les données de l'étude i
  subset_study <- validation_data %>%
    filter(
      Pathogen == studies$Pathogen[i],
      Country  == studies$Country[i],
      Author   == studies$Author[i]
    )

  # Extraire les ICC intervals (colonne 5)
  icc <- subset_study[5] %>% unlist(use.names = FALSE)
  icc <- icc[!is.na(icc)]

  if (length(icc) == 0) next

  # Titre du plot
  title <- paste0(
    studies$Pathogen[i], " | ",
    studies$Country[i],  " | ",
    studies$Author[i]
  )

  # Histogramme
  plot_list[[i]] <- ggplot(data.frame(icc = icc), aes(x = icc)) +
    geom_histogram(
      breaks   = seq(min(icc) - 0.5, max(icc) + 0.5, by = 1),
      fill     = "lightblue",
      color    = "black"
    ) +
    labs(
      title = title,
      x     = "Index-case to case interval (days)",
      y     = "Count"
    ) +
    theme_minimal()

  print(plot_list[[i]])
}
