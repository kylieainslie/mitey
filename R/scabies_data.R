# Load scabies ICC interval data
library(here)
library(purrr)
file_path <-  here("vignettes", "data", "si_data.rds")
scabies_si_data <- readRDS(file_path)

scabies_si_data


# Estimate serial interval for each study
result_by_study <- scabies_si_data %>%
  group_by(study) %>%
  summarise(result = list(si_estim(icc_interval))) %>%
  mutate(
    mean = map_dbl(result, "mean"),
    sd = map_dbl(result, "sd"),
    wts = map(result, "wts")
  ) %>%
  select(-result)

# Display results
result_by_study %>%
  select(study, mean, sd) %>%
  mutate(across(c(mean, sd), round, 2)) %>%
  arrange(mean) %>%
  knitr::kable(caption = "Estimated mean and standard deviation of serial interval (days) by study")




#plot scabies hist
# définir les colonnes qui identifient une étude
list_of_studies_scabies <- scabies_si_data %>%
  select(study) %>%
  distinct()

list_of_studies_scabies
# boucle
for(i in 1:nrow(list_of_studies_scabies)) {

  study1 <- list_of_studies_scabies[i, ]

  subset_data <- scabies_si_data %>%
    filter(study == study1)

  icc <- subset_data[[4]]  # colonne ICC

  df <- data.frame(value = icc)

  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(
      title = paste("ICC -"),
      x = "ICC",
      y = "Fréquence"
    ) +
    theme_minimal()

  print(p)
}

scabies_subset <- scabies_si_data%>% filter (study == "Tjon-Kon-Fat et al")
icc_scabies<-scabies_subset[4] %>% unlist(,use.names = FALSE)
result_scabies<-si_estim(icc_scabies)

plot_si_fit_result(result_scabies,icc_scabies)
