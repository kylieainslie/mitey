library(extraDistr)

# ============================================================================
# PARAMÈTRES GÉNÉRAUX
# ============================================================================
N <- 10000
mu <- 22      # moyenne de l'intervalle réel entre événements
sigma <- 9    # écart-type de l'intervalle réel

# Intervalle d'observation (en jours)
# Par exemple : 1 = quotidien, 7 = hebdomadaire, 30 = mensuel
intervalle_observation <- 7  # données récupérées 1 fois par semaine

# ============================================================================
# FONCTION GÉNÉRALE DE TRANSFORMATION
# ============================================================================
# Cette fonction simule l'incertitude due à l'intervalle d'observation
transformer_selon_intervalle <- function(x, intervalle = 7) {
  # Calcule le reste de la division par l'intervalle d'observation
  reste <- x %% intervalle

  # Arrondit vers le haut pour les valeurs proches de la limite
  # (simule le fait qu'on rate les événements entre deux observations)
  x_transforme <- ifelse(reste > 0 & x > intervalle,
                         x + (intervalle - reste),  # arrondit au prochain multiple
                         x)
  return(x_transforme)
}

# ============================================================================
# DONNÉES ORIGINALES (données "vraies" avec timing exact)
# ============================================================================
cat("=== DONNÉES ORIGINALES (Timing exact) ===\n\n")

# P1 : temps d'apparition du 1er événement (jours, entre 1 et 14)
P1 <- rdunif(n = N, min = 1, max = 14)

# norm : intervalle réel entre les deux événements (loi normale tronquée)
norm <- trunc(rnorm(n = N, mean = mu, sd = sigma))

# P2 : temps d'apparition du 2e événement
P2 <- P1 + norm

# I : intervalle réel entre les deux événements
I <- P2 - P1


# Visualisation données originales

breaks1 <- seq(0.5, max(P1) + 0.5, 1)
hist(P1, main = "P1 : 1er événement (Exact)",
     xlab = "Jours", col = "skyblue", breaks = breaks1)
hist(P2, main = "P2 : 2e événement (Exact)",
     xlab = "Jours", col = "lightgreen", breaks = 40)
hist(I, main = "I : Intervalle réel (Exact)",
     xlab = "Jours", col = "cyan", breaks = 40)


# ============================================================================
# DONNÉES OBSERVÉES (avec incertitude d'observation)
# ============================================================================
cat("=== DONNÉES OBSERVÉES (Intervalle observation: ", intervalle_observation, " jours) ===\n\n")

# Appliquer la transformation
P1_obs <- transformer_selon_intervalle(P1, intervalle = intervalle_observation)
P2_obs <- transformer_selon_intervalle(P2, intervalle = intervalle_observation)
I_obs <- P2_obs - P1_obs


# Visualisation données observées

breaks1 <- seq(0.5, max(P1_obs) + 0.5, intervalle_observation)
hist(P1_obs, main = paste("P1_obs : 1er événement (Obs. tous les", intervalle_observation, "j)"),
     xlab = "Jours", col = "lightblue", breaks = breaks1)
hist(P2_obs, main = paste("P2_obs : 2e événement (Obs. tous les", intervalle_observation, "j)"),
     xlab = "Jours", col = "lightyellow", breaks = 40)
hist(I_obs, main = "I_obs : Intervalle observé",
     xlab = "Jours", col = "lightcoral", breaks = 40)

mean(I)
sd(I)
mean(I_obs)
sd(I_obs)
