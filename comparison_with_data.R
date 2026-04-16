library(dplyr)
library(purrr)
library(tidyr)

# 1. Chargement des données
# On utilise le chemin relatif tel que défini dans votre fichier essai_avec_data.R
validation_data <- readRDS("vignettes/articles/validation_data.rds")

# 2. Définition manuelle des cibles
targets <- tribble(
  ~Pathogen,                ~Country,  ~Author,
  "Measles",                "Kenya",   "Aaby",
  "Influenza A(H3N2)",      "France",  "Viboud",
  "Influenza A(H1N1)pdm09", "Canada",  "Papenbourg",
  "Influenza A(H1N1)pdm09", "USA",     "France",
  "Influenza A(H1N1)pdm09", "USA",     "Morgan",
  "Influenza A(H1N1)pdm09", "USA",     "Cauchemez"
)

# 3. Fonction d'analyse par route
analyze_routes <- function(p, c, a, routes_to_test = 2:6) {
  # Filtrage des données (colonne 5 pour les ICC)
  subset_data <- validation_data %>%
    filter(Pathogen == p, Country == c, Author == a)

  icc_values <- subset_data[[5]] %>% unlist(use.names = FALSE)

  if (length(icc_values) < 2) return(NULL)

  message(paste("Analyse de", p, "au", c, "..."))

  map_df(routes_to_test, function(nr) {
    res <- tryCatch({
      # Utilisation de si_estim du package mitey
      si_estim(icc_values, n_routes = nr, dist = "normal")
    }, error = function(e) return(NULL))

    if (is.null(res)) return(data.frame())

    data.frame(
      n_routes = nr,
      Mean     = round(res$mean, 3),
      SD       = round(res$sd, 2),
      LogLik   = round(res$loglik, 2),
      AIC      = round((2 * (2 * nr)) - (2 * res$loglik), 2),
      Converged = res$converged
    )
  })
}

# 4. Exécution et séparation en plusieurs tableaux
all_results <- targets %>%
  pmap(function(Pathogen, Country, Author) {
    res_table <- analyze_routes(Pathogen, Country, Author)
    if (!is.null(res_table)) {
      # On ajoute le nom pour identifier le tableau dans la liste
      attr(res_table, "title") <- paste(Pathogen, "-", Country)
      return(res_table)
    }
    return(NULL)
  }) %>%
  compact() # Supprime les éléments vides

# 5. Affichage des tableaux individuels
for (tbl in all_results) {
  cat("\n" , paste(rep("=", 30), collapse = ""), "\n")
  cat("TABLEAU :", attr(tbl, "title"), "\n")
  cat(paste(rep("-", 30), collapse = ""), "\n")
  print(tbl)
}
