library(extraDistr)

# ============================================================================
# PARAMÈTRES GÉNÉRAUX
# ============================================================================
N <- 300
mu <- 33      # moyenne d'intervalle de référence
sigma <- 10    # écart-type de référence

# Intervalle d'observation (en jours)
intervalle_observation <- 7

# ============================================================================
# PARAMÈTRES DES ROUTES
# ============================================================================
# Chaque route a une distribution différente d'intervalle
routes <- list(
  route1 = list(mu = mu, sigma = sigma, label = "Route 1 "),
  route2 = list(mu = 2*mu, sigma = sqrt(2)*sigma, label = "Route 2 "),
  route3 = list(mu = 3*mu, sigma = sqrt(3)*sigma, label = "Route 3 "),
  route4 = list(mu = 4*mu, sigma = sqrt(4)*sigma, label = "Route 4 ")
)

# POIDS DES ROUTES (doivent sommer à 1)
# Par exemple : 50% route1, 30% route2, 15% route3, 5% route4
poids <- c(
  route1 = 0.70,
  route2 = 0.15,
  route3 = 0.1,
  route4 = 0.05
)

# Vérifier que les poids somment à 1
if (abs(sum(poids) - 1) > 1e-10) {
  stop("ERREUR : Les poids ne somment pas à 1. Somme actuelle:", sum(poids))
}


# ============================================================================
# GÉNÉRATION DU DATASET
# ============================================================================

# P1 : temps d'apparition du 1er événement
P1 <- rdunif(n = N, min = 1, max = 14)

# Créer L : mélange des intervalles selon les routes et poids
L <- numeric(N)
route_assignment <- numeric(N)  # Pour tracker quelle route a été prise

# Assigner chaque observation à une route selon les poids
routes_indices <- sample(1:4, size = N, replace = TRUE, prob = poids)

for (i in 1:N) {
  route_idx <- routes_indices[i]
  route_name <- names(routes)[route_idx]
  route_params <- routes[[route_idx]]



  # Générer l'intervalle selon les paramètres de la route
  L[i] <- trunc(rnorm(n = 1, mean = route_params$mu, sd = route_params$sigma))

  route_assignment[i] <- route_idx
}

# P2 : temps d'apparition du 2e événement
P2 <- P1 + L

# I : intervalle réel
I <- P2 - P1

# ============================================================================
# AFFICHAGE DES STATISTIQUES
# ============================================================================
cat("=== DONNÉES ORIGINALES (Timing exact) ===\n\n")
cat("P1 (1er événement)   - Moy:", round(mean(P1), 2),
    "Min:", round(min(P1), 2), "Max:", round(max(P1), 2), "\n")
cat("P2 (2e événement)    - Moy:", round(mean(P2), 2),
    "Min:", round(min(P2), 2), "Max:", round(max(P2), 2), "\n")
cat("I (intervalle réel)  - Moy:", round(mean(I), 2),
    "Écart-type:", round(sd(I), 2),
    "Min:", round(min(I), 2), "Max:", round(max(I), 2), "\n\n")


# ============================================================================
# FONCTION DE TRANSFORMATION (INCERTITUDE D'OBSERVATION)
# ============================================================================
transformer_selon_intervalle <- function(x, intervalle = 7) {
  reste <- x %% intervalle
  x_transforme <- ifelse(reste > 0,
                         x + (intervalle - reste),
                         x)
  return(x_transforme)
}

# ============================================================================
# DONNÉES OBSERVÉES
# ============================================================================
cat("\n=== DONNÉES OBSERVÉES (Intervalle observation: ", intervalle_observation, " jours) ===\n\n")

P1_obs <- transformer_selon_intervalle(P1, intervalle = intervalle_observation)
P2_obs <- transformer_selon_intervalle(P2, intervalle = intervalle_observation)
I_obs <- P2_obs - P1_obs

cat("P1_obs (1er événement)   - Moy:", round(mean(P1_obs), 2), "\n")
cat("P2_obs (2e événement)    - Moy:", round(mean(P2_obs), 2), "\n")
cat("I_obs (intervalle obs)   - Moy:", round(mean(I_obs), 2),
    "Écart-type:", round(sd(I_obs), 2), "\n\n")

# ============================================================================
# VISUALISATIONS
# ============================================================================
par(mfrow = c(2, 3))

# Données exactes
breaks1 <- seq(0.5, max(P1) + 0.5, 1)
hist(P1, main = "P1 : 1er événement (Exact)",
     xlab = "Jours", col = "skyblue", breaks = breaks1)

hist(P2, main = "P2 : 2e événement (Exact)",
     xlab = "Jours", col = "lightgreen", breaks = 40)

hist(I, main = "I : Intervalle réel (Exact)",
     xlab = "Jours", col = "cyan", breaks = 40)

# Données observées
hist(P1_obs, main = paste("P1_obs : Obs. tous les", intervalle_observation, "j"),
     xlab = "Jours", col = "lightblue", breaks = breaks1)

hist(P2_obs, main = paste("P2_obs : Obs. tous les", intervalle_observation, "j"),
     xlab = "Jours", col = "lightyellow", breaks = 40)

hist(I_obs, main = "I_obs : Intervalle observé",
     xlab = "Jours", col = "lightcoral", breaks = 40)

par(mfrow = c(1, 1))

# Comparaison des distributions d'intervalle


plot(density(I), main = "Distribution de l'intervalle (Exact)",
     xlab = "Jours", col = "blue", lwd = 2)
lines(density(I_obs), col = "red", lwd = 2, lty = 2)
legend("topright", c("Exact", "Observé"), col = c("blue", "red"), lwd = 2, lty = c(1, 2))

# Distribution par route
colors <- c("red", "orange", "green", "purple")
plot(density(L[route_assignment == 1]), main = "Distribution des intervalles par route",
     xlab = "Jours", col = colors[1], lwd = 2, xlim = c(0, max(L)))
for (r in 2:4) {
  if (sum(route_assignment == r) > 1) {
    lines(density(L[route_assignment == r]), col = colors[r], lwd = 2)
  }
}
legend("topright", paste("Route", 1:4, sprintf("(%.0f%%)", poids*100)),
       col = colors, lwd = 2)


# ============================================================================
# RÉSUMÉ COMPARATIF
# ============================================================================
cat("\n=== RÉSUMÉ COMPARATIF ===\n\n")
cat("Métrique                 | Exact    | Observé  | Erreur\n")
cat(sprintf("-%.0s", 1:55), "\n")
cat(sprintf("Moyenne intervalle       | %8.2f | %8.2f | %+6.2f\n",
            mean(I), mean(I_obs), mean(I_obs) - mean(I)))
cat(sprintf("Écart-type intervalle    | %8.2f | %8.2f | %+6.2f\n",
            sd(I), sd(I_obs), sd(I_obs) - sd(I)))
cat(sprintf("Médiane intervalle       | %8.2f | %8.2f | %+6.2f\n",
            median(I), median(I_obs), median(I_obs) - median(I)))

# ============================================================================
# EXPORT DU DATASET
# ============================================================================
dataset <- data.frame(
  P1 = P1,
  P1_obs = P1_obs,
  P2 = P2,
  P2_obs = P2_obs,
  L = L,
  I = I,
  I_obs = I_obs
)

# Aperçu du dataset
cat("\n=== APERÇU DU DATASET ===\n\n")
print(head(dataset, 10))
cat("\nDimensions:", nrow(dataset), "lignes,", ncol(dataset), "colonnes\n")


ICC_dataset<-dataset$I_obs

ICC_real <-dataset$I
