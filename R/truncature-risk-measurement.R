# DÉTECTION DE TRONCATURE À DROITE - AVEC SEULEMENT LES INTERVALLES
# (sans dates ni contexte temporel)

library(ggplot2)
library(dplyr)

# ============================================================================
# TEST 1: ANALYSE DE LA DISTRIBUTION ELLE-MÊME
# ============================================================================
# Idée : si troncature, la queue droite est coupée
# → distribution asymétrique, skewness positif anormal, max trop bas

test_distribution_shape <- function(delays) {

  # Statistiques descriptives
  mean_d <- mean(delays, na.rm=TRUE)
  median_d <- median(delays, na.rm=TRUE)
  sd_d <- sd(delays, na.rm=TRUE)
  max_d <- max(delays, na.rm=TRUE)
  min_d <- min(delays, na.rm=TRUE)

  # Skewness (asymétrie) : mesure la queue
  # skewness > 1 = queue droite prononcée (attendu pour délais)
  # skewness > 2 = très asymétrique (possible troncature!)

  skewness <- function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    m <- mean(x)
    s <- sd(x)
    m3 <- mean((x - m)^3)
    return(m3 / s^3)
  }

  skew <- skewness(delays)

  # Ratio : écart-type / moyenne
  # Délais naturels : ratio ~0.8-1.5
  # Si ratio << 1 (ex: 0.3) = peu de variabilité = troncature probable!

  cv <- sd_d / mean_d

  # Ratio : max / médiane
  # Attendu : 2-4× (délais ont longue queue)
  # Si < 2 = queue coupée?

  ratio_max_median <- max_d / median_d

  cat("\n════════════════════════════════════════════════════════\n")
  cat("TEST 1: ANALYSE DE LA DISTRIBUTION\n")
  cat("════════════════════════════════════════════════════════\n\n")

  cat("Statistiques de base:\n")
  cat("  Minimum        :", round(min_d, 1), "j\n")
  cat("  Q1 (25%)       :", round(quantile(delays, 0.25), 1), "j\n")
  cat("  Médiane        :", round(median_d, 1), "j\n")
  cat("  Moyenne        :", round(mean_d, 1), "j\n")
  cat("  Q3 (75%)       :", round(quantile(delays, 0.75), 1), "j\n")
  cat("  Maximum        :", round(max_d, 1), "j\n")
  cat("  Écart-type     :", round(sd_d, 1), "j\n\n")

  cat("Indicateurs de TRONCATURE à droite:\n")
  cat("  1. Coefficient de variation (CV) :", round(cv, 2), "\n")
  cat("     Interprétation:\n")
  cat("       - 0.8-1.5  : Normal\n")
  cat("       - 0.4-0.8  : ⚠️ Troncature POSSIBLE\n")
  cat("       - < 0.4    : 🚨 Troncature PROBABLE\n\n")

  cat("  2. Skewness (asymétrie) :", round(skew, 2), "\n")
  cat("     Interprétation:\n")
  cat("       - 1.0-2.0  : Normal (longue queue attendue)\n")
  cat("       - 0.3-1.0  : ⚠️ Queue courte, possible troncature\n")
  cat("       - < 0.3    : 🚨 Très peu de variabilité\n\n")

  cat("  3. Ratio max/médiane :", round(ratio_max_median, 1), "x\n")
  cat("     Interprétation:\n")
  cat("       - 2.5-4.0  : Normal\n")
  cat("       - 1.5-2.5  : ⚠️ Queue courte, possible troncature\n")
  cat("       - < 1.5    : 🚨 Troncature PROBABLE\n\n")

  # Score composite
  score <- 0
  if (cv < 0.4) score <- score + 3
  else if (cv < 0.8) score <- score + 1.5

  if (skew < 0.3) score <- score + 3
  else if (skew < 1.0) score <- score + 1.5

  if (ratio_max_median < 1.5) score <- score + 3
  else if (ratio_max_median < 2.5) score <- score + 1.5

  cat("SCORE DE RISQUE (0-9):", round(score, 1), "\n")
  cat("  0-2    : ✓ Peu probable\n")
  cat("  2-4    : ⚠️  Modéré\n")
  cat("  4-6    : 🚨 Important\n")
  cat("  > 6    : 🚨🚨 Très important\n\n")

  return(list(
    cv = cv,
    skewness = skew,
    ratio_max_median = ratio_max_median,
    risk_score = score
  ))
}

# ============================================================================
# TEST 2: COMPARAISON AVEC DISTRIBUTIONS THÉORIQUES
# ============================================================================
# Fit plusieurs distributions et voir lesquelles s'ajustent bien
# Si aucune n'ajuste la queue droite = troncature

test_goodness_of_fit <- function(delays) {

  library(fitdistrplus)

  cat("\n════════════════════════════════════════════════════════\n")
  cat("TEST 2: AJUSTEMENT À DES DISTRIBUTIONS THÉORIQUES\n")
  cat("════════════════════════════════════════════════════════\n\n")

  # Fit plusieurs distributions
  fit_gamma <- fitdist(delays, "gamma", method="mle")
  fit_lnorm <- fitdist(delays, "lnorm", method="mle")
  fit_weibull <- fitdist(delays, "weibull", method="mle")

  # Goodness-of-fit tests
  # Kolmogorov-Smirnov : teste si données = distribution
  # p-value élevée = bon ajustement

  ks_gamma <- ks.test(delays, "pgamma",
                      shape=fit_gamma$estimate["shape"],
                      rate=fit_gamma$estimate["rate"])

  ks_lnorm <- ks.test(delays, "plnorm",
                      meanlog=fit_lnorm$estimate["meanlog"],
                      sdlog=fit_lnorm$estimate["sdlog"])

  ks_weibull <- ks.test(delays, "pweibull",
                        shape=fit_weibull$estimate["shape"],
                        scale=fit_weibull$estimate["scale"])

  cat("Test de Kolmogorov-Smirnov (p > 0.05 = bon ajustement):\n\n")

  cat("Gamma:\n")
  cat("  p-value :", format.pval(ks_gamma$p.value, digits=3), "\n")
  cat("  shape =", round(fit_gamma$estimate["shape"], 2),
      ", rate =", round(fit_gamma$estimate["rate"], 2), "\n\n")

  cat("Lognormal:\n")
  cat("  p-value :", format.pval(ks_lnorm$p.value, digits=3), "\n")
  cat("  meanlog =", round(fit_lnorm$estimate["meanlog"], 2),
      ", sdlog =", round(fit_lnorm$estimate["sdlog"], 2), "\n\n")

  cat("Weibull:\n")
  cat("  p-value :", format.pval(ks_weibull$p.value, digits=3), "\n")
  cat("  shape =", round(fit_weibull$estimate["shape"], 2),
      ", scale =", round(fit_weibull$estimate["scale"], 2), "\n\n")

  cat("Interprétation:\n")
  cat("  Si p-value < 0.05 pour TOUTES => mauvais ajustement partout\n")
  cat("  => possible que la queue droite soit coupée (troncature)\n\n")

  return(list(
    fit_gamma = fit_gamma,
    fit_lnorm = fit_lnorm,
    fit_weibull = fit_weibull,
    ks_results = list(gamma = ks_gamma, lnorm = ks_lnorm, weibull = ks_weibull)
  ))
}

# ============================================================================
# TEST 3: Q-Q PLOT - VISUALISATION GRAPHIQUE
# ============================================================================
# Si troncature = la queue droite s'écarte de la ligne diagonale

plot_qq_diagnosis <- function(delays) {

  p1 <- ggplot(data.frame(x = delays), aes(sample = x)) +
    stat_qq(distribution = qnorm, color = "#185FA5", size = 2, alpha = 0.6) +
    stat_qq_line(distribution = qnorm, color = "#E24B4A", linewidth = 1) +
    labs(title = "Q-Q plot vs Normal", x = "Théorique", y = "Observé") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11))

  p2 <- ggplot(data.frame(x = log(delays)), aes(sample = x)) +
    stat_qq(distribution = qnorm, color = "#185FA5", size = 2, alpha = 0.6) +
    stat_qq_line(distribution = qnorm, color = "#E24B4A", linewidth = 1) +
    labs(title = "Q-Q plot vs Normal (log-échelle)", x = "Théorique", y = "Observé") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11))

  p3 <- ggplot(data.frame(delay = delays), aes(x = delay)) +
    geom_histogram(bins = 30, fill = "#378ADD", alpha = 0.7,
                   color = "#185FA5", linewidth = 0.5) +
    labs(title = "Histogramme des intervalles", x = "Intervalle (jours)", y = "Fréquence") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 11))

  gridExtra::grid.arrange(p1, p2, p3, nrow = 1)
}

# ============================================================================
# TEST 4: ANALYSE DE LA QUEUE (TAIL INDEX)
# ============================================================================
# Les intervalles sériels ont naturellement une longue queue (queue-heavy)
# Si queue-index anormalement bas = troncature

test_tail_heaviness <- function(delays) {

  cat("\n════════════════════════════════════════════════════════\n")
  cat("TEST 4: ANALYSE DE LA QUEUE (TAIL HEAVINESS)\n")
  cat("════════════════════════════════════════════════════════\n\n")

  # Hill estimator pour index de queue
  # Ordonne les valeurs et examine les plus grandes

  sorted <- sort(delays, decreasing = TRUE)
  n <- length(sorted)

  # Nombre de points dans la queue (heuristique : sqrt(n))
  k <- max(10, floor(sqrt(n)))

  # Hill estimator
  log_ratio <- log(sorted[1:k] / sorted[k+1])
  tail_index <- mean(log_ratio)

  cat("Hill tail index (pour les", k, "plus grandes valeurs):",
      round(tail_index, 2), "\n\n")

  cat("Interprétation:\n")
  cat("  > 2.0  : Queue très légère (lognormal)\n")
  cat("  1.0-2.0: Queue modérée (gamma)\n")
  cat("  0.5-1.0: Queue lourde (Pareto)\n")
  cat("  < 0.5  : ⚠️ Queue très légère ou TRONCATURE!\n\n")

  # Proportion dans la queue
  p90 <- quantile(delays, 0.90)
  p95 <- quantile(delays, 0.95)
  p99 <- quantile(delays, 0.99)

  cat("Percentiles extrêmes:\n")
  cat("  90e percentile :", round(p90, 1), "j\n")
  cat("  95e percentile :", round(p95, 1), "j\n")
  cat("  99e percentile :", round(p99, 1), "j\n\n")

  prop_above_p90 <- mean(delays > p90)
  prop_above_p95 <- mean(delays > p95)

  cat("Proportion au-delà des percentiles:\n")
  cat("  Attendu > P90 : ~10%, observé :", round(100*prop_above_p90, 1), "%\n")
  cat("  Attendu > P95 : ~5%, observé :", round(100*prop_above_p95, 1), "%\n\n")

  if (prop_above_p95 < 0.02) {
    cat("🚨 Seulement", round(100*prop_above_p95, 1),
        "% au-delà P95 (attendu 5%)\n")
    cat("   => Possible troncature à droite!\n\n")
  }

  return(list(
    tail_index = tail_index,
    p90 = p90,
    p95 = p95,
    p99 = p99,
    prop_above_p95 = prop_above_p95
  ))
}

# ============================================================================
# FONCTION PRINCIPALE : RAPPORT COMPLET
# ============================================================================

diagnose_truncation <- function(delays) {

  cat("\n\n")
  cat("╔════════════════════════════════════════════════════════╗\n")
  cat("║  DIAGNOSTIC DE TRONCATURE À DROITE - INTERVALLES SEULS║\n")
  cat("╚════════════════════════════════════════════════════════╝\n")

  cat("\nVotre dataset :\n")
  cat("  n =", length(delays), "intervalles\n")
  cat("  Plage : [", min(delays, na.rm=T), ",", max(delays, na.rm=T), "] jours\n\n")

  # TEST 1
  result1 <- test_distribution_shape(delays)

  # TEST 4 (rapide, pas de fitting)
  result4 <- test_tail_heaviness(delays)

  # Décision finale
  cat("\n════════════════════════════════════════════════════════\n")
  cat("CONCLUSION\n")
  cat("════════════════════════════════════════════════════════\n\n")

  score <- result1$risk_score

  if (score < 2) {
    cat("✓ TRONCATURE IMPROBABLE\n")
    cat("  → Vous pouvez utiliser les données telles quelles\n")
    cat("  → Fit lognormal standard recommandé\n")
  } else if (score < 4) {
    cat("⚠️  TRONCATURE MODÉRÉE (risque faible-modéré)\n")
    cat("  → Vérifiez visuellement avec les Q-Q plots\n")
    cat("  → Considérez une correction légère (filtrage)\n")
  } else if (score < 6) {
    cat("🚨 TRONCATURE IMPORTANTE (risque modéré-élevé)\n")
    cat("  → RECOMMANDÉ: correction par régression Weibull\n")
    cat("  → ou fit de mélange (mixture model)\n")
  } else {
    cat("🚨🚨 TRONCATURE SÉVÈRE (risque très élevé)\n")
    cat("  → CRITIQUE: vos données sont fortement biaisées\n")
    cat("  → Nécessite approche paramétrique robuste\n")
    cat("  → Considérez de recueillir plus de données\n")
  }

  cat("\n\nProchain pas recommandé:\n")
  if (score > 3) {
    cat("  1. Exécuter test_goodness_of_fit()\n")
    cat("  2. Exécuter plot_qq_diagnosis()\n")
    cat("  3. Considérer la correction (voir script correction_truncation.R)\n")
  } else {
    cat("  Vos données semblent OK. Vous pouvez utiliser\n")
    cat("  fitdist() ou fit_disc_gamma_mle() directement.\n")
  }

  return(invisible(list(
    dist_shape = result1,
    tail = result4,
    risk_score = score
  )))
}

# ============================================================================
# EXEMPLE D'UTILISATION
# ============================================================================

# Générer des données de test
set.seed(42)

# Scénario 1: Pas de troncature (données complètes)
delays_complete <- rlnorm(200, meanlog=0.8, sdlog=0.6)

# Scénario 2: Avec troncature (queue coupée)
delays_true <- rlnorm(200, meanlog=0.8, sdlog=0.6)
delays_truncated <- delays_true[delays_true <= 10]  # couper à 10 jours

cat("\n\n")
cat("════════════════════════════════════════════════════════\n")
cat("EXEMPLE 1: DONNÉES SANS TRONCATURE\n")
cat("════════════════════════════════════════════════════════\n")

diagnose_truncation(delays_complete)

cat("\n\n")
cat("════════════════════════════════════════════════════════\n")
cat("EXEMPLE 2: DONNÉES AVEC TRONCATURE\n")
cat("════════════════════════════════════════════════════════\n")

diagnose_truncation(delays_truncated)



