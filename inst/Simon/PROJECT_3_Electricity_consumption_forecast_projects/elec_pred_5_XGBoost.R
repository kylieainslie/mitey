#-------------------------------------------------------------------------------
# 1. BIBLIOTHÈQUES ET CONFIGURATION
#-------------------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(lubridate)
library(forecast)
library(xgboost)
library(MLmetrics)
library(readr)
library(pdp)
library(zoo)

setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/inst/Simon/PROJECT_3_Electricity_consumption_forecast_projects")

#-------------------------------------------------------------------------------
# 2. CHARGEMENT ET NETTOYAGE DES DONNÉES
#-------------------------------------------------------------------------------
# Chargement
df_elec  <- read.csv("consommation-quotidienne-brute.csv", sep = ";")
df_meteo <- read.csv("meteo_paris_2025.csv", sep = ",")

# Nettoyage Consommation
df_hourly <- df_elec %>%
  select(Date = 1, Value = 9) %>%
  mutate(
    Date = as.POSIXct(gsub("T", " ", substr(Date, 1, 16)), format = "%Y-%m-%d %H:%M", tz = "Europe/Paris"),
    Hour = as.POSIXct(format(Date, "%Y-%m-%d %H:00:00"), tz = "Europe/Paris")
  ) %>%
  filter(!is.na(Date)) %>%
  group_by(Hour) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()

# Nettoyage Météo
df_meteo <- df_meteo %>%
  mutate(Hour = parse_date_time(time, orders = c("ymd HMS", "ymd"), tz = "Europe/Paris")) %>%
  distinct(Hour, .keep_all = TRUE) %>%
  select(Hour, temp = temperature_2m)

# Fusion et Interpolation
df_final <- data.frame(Hour = seq(min(df_hourly$Hour), max(df_hourly$Hour), by = "hour")) %>%
  left_join(df_hourly, by = "Hour") %>%
  left_join(df_meteo, by = "Hour") %>%
  filter(Hour >= as.POSIXct("2025-01-01 00:00:00", tz = "Europe/Paris")) %>%
  mutate(
    Value = as.numeric(na.interp(Value / 1000)), # Conversion MW -> GW
    temp  = as.numeric(na.interp(temp))
  )

#-------------------------------------------------------------------------------
# 3. FEATURE ENGINEERING (LE CŒUR DU MODÈLE)
#-------------------------------------------------------------------------------
# A. Multi-Saisonnalité de Fourier
msts_obj <- msts(df_final$Value, seasonal.periods = c(24, 168))
fourier_df <- as.data.frame(fourier(msts_obj, K = c(10, 5)))

# B. Création de toutes les variables explicatives
df_xgb_data <- df_final %>%
  cbind(fourier_df) %>%
  arrange(Hour) %>%
  mutate(
    # Variables Temporelles Cycliques
    sin_hour = sin(2 * pi * hour(Hour) / 24),
    cos_hour = cos(2 * pi * hour(Hour) / 24),
    Hour_v   = hour(Hour),
    Wday     = wday(Hour),
    Month    = month(Hour),
    IsWeekend = ifelse(Wday %in% c(1, 7), 1, 0),

    # Variables Thermiques (Physique du bâtiment)
    Froid       = pmax(0, 15 - temp),
    Froid_2     = Froid^2,
    Grand_Froid = pmax(0, 5 - temp),
    Chaud       = pmax(0, temp - 25),
    temp_roll3  = as.numeric(rollapplyr(temp, 3, mean, fill = "extend")),

    # Lags (Mémoire du réseau)
    Lag24  = lag(Value, 24),
    Lag168 = lag(Value, 168)
  ) %>%
  filter(!is.na(Lag168)) # On retire le début sans historique

#-------------------------------------------------------------------------------
# 4. PRÉPARATION ET ENTRAÎNEMENT
#-------------------------------------------------------------------------------
# Liste propre des features (On inclut tout ce qui a survécu à tes tests)
features <- c("temp", "temp_roll3", "Froid", "Froid_2",
              "sin_hour", "cos_hour", "Lag24", "Lag168",
              "Hour_v", "Wday", "IsWeekend", colnames(fourier_df))

# Split Train/Test (Dernières 24h)
h <- 24
n <- nrow(df_xgb_data)
train_df <- df_xgb_data[1:(n - h), ]
test_df  <- df_xgb_data[(n - h + 1):n, ]

train_x <- as.matrix(train_df[, features]); train_y <- train_df$Value
test_x  <- as.matrix(test_df[, features]);  test_y  <- test_df$Value

# Modèle XGBoost
fit_xgb <- xgboost(
  data = train_x, label = train_y,
  nrounds = 1000,
  max_depth = 6,
  eta = 0.05,
  subsample = 0.7,
  colsample_bytree = 0.7,
  min_child_weight = 15,
  objective = "reg:squarederror",
  verbosity = 0
)

#-------------------------------------------------------------------------------
# 5. RÉSULTATS ET DIAGNOSTICS
#-------------------------------------------------------------------------------
preds <- predict(fit_xgb, test_x)

# Métriques
cat("\n--- RÉSULTATS ---\n")
cat("MAPE Test :", round(MAPE(preds, test_y) * 100, 2), "%\n")
cat("MAPE Train:", round(MAPE(predict(fit_xgb, train_x), train_y) * 100, 2), "%\n")

# Importance des variables
xgb.importance(feature_names = features, model = fit_xgb) %>%
  xgb.plot.importance(top_n = 10)

# Visualisation Prévision
results_plot <- data.frame(Hour = test_df$Hour, Reel = test_y, Pred = preds)
ggplot(results_plot, aes(x = Hour)) +
  geom_line(aes(y = Reel, color = "Réalité"), size = 1) +
  geom_line(aes(y = Pred, color = "XGBoost"), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Réalité" = "black", "XGBoost" = "red")) +
  theme_minimal() + labs(title = "Prévision 24h", y = "GW")

# Analyse des Résidus
residus <- test_y - preds
plot(test_df$temp, residus, main="Résidus vs Température", xlab="Temp", ylab="Erreur (GW)")
abline(h=0, col="red")

# PDP : Influence Température
partial(fit_xgb, pred.var = "temp", train = train_x, grid.resolution = 30) %>%
  autoplot(color = "darkorange", size = 1) + theme_minimal() + labs(title = "PDP Température")
