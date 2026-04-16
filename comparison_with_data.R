library(dplyr)
library(purrr)
library(tidyr)

# 1. Data loading
# Use the relative path as defined in your essai_avec_data.R file
validation_data <- readRDS("vignettes/articles/validation_data.rds")

# 2. Manual definition of targets
targets <- tribble(
  ~Pathogen,                ~Country,  ~Author,
  "Measles",                "Kenya",   "Aaby",
  "Influenza A(H3N2)",      "France",  "Viboud",
  "Influenza A(H1N1)pdm09", "Canada",  "Papenbourg",
  "Influenza A(H1N1)pdm09", "USA",     "France",
  "Influenza A(H1N1)pdm09", "USA",     "Morgan",
  "Influenza A(H1N1)pdm09", "USA",     "Cauchemez"
)

# 3. Route analysis function
analyze_routes <- function(p, c, a, routes_to_test = 2:6) {
  # Data filtering (column 5 for ICC)
  subset_data <- validation_data %>%
    filter(Pathogen == p, Country == c, Author == a)

  icc_values <- subset_data[[5]] %>% unlist(use.names = FALSE)

  if (length(icc_values) < 2) return(NULL)

  message(paste("Analyse de", p, "au", c, "..."))

  map_df(routes_to_test, function(nr) {
    res <- tryCatch({
      # Use si_estim from the mitey package
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

# 4. Execution and separation into multiple tables
all_results <- targets %>%
  pmap(function(Pathogen, Country, Author) {
    res_table <- analyze_routes(Pathogen, Country, Author)
    if (!is.null(res_table)) {
      # Add the name to identify the table in the list
      attr(res_table, "title") <- paste(Pathogen, "-", Country)
      return(res_table)
    }
    return(NULL)
  }) %>%
  compact() # Remove empty elements

# 5. Display individual tables
for (tbl in all_results) {
  cat("\n" , paste(rep("=", 30), collapse = ""), "\n")
  cat("TABLEAU :", attr(tbl, "title"), "\n")
  cat(paste(rep("-", 30), collapse = ""), "\n")
  print(tbl)
}
