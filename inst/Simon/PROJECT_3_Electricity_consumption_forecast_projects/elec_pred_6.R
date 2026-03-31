#-------------------------------------------------------------------------------
# 1. BIBLIOTHÈQUES ET ENVIRONNEMENT
#-------------------------------------------------------------------------------
library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)
library(prophet)
library(lubridate)
library(readr)
library(tictoc)
library(xgboost)
library(MLmetrics)
library(zoo)

model_times <- list()
setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/inst/Simon/PROJECT_3_Electricity_consumption_forecast_projects")

#-------------------------------------------------------------------------------
# 2. CHARGEMENT ET NETTOYAGE DES DONNÉES
#-------------------------------------------------------------------------------
df_elec  <- read.csv("consommation-quotidienne-brute.csv", sep = ";")
df_meteo <- read.csv("meteo_paris_2025.csv", sep = ",")

# Nettoyage Elec
df_elec <- df_elec[, c(1, 9)]
colnames(df_elec) <- c("Date", "Value")
df_elec$Date <- as.POSIXct(gsub("T", " ", substr(df_elec$Date, 1, 16)),
                           format = "%Y-%m-%d %H:%M", tz = "Europe/Paris")
df_elec <- na.omit(df_elec) %>% arrange(Date)
df_elec$Hour <- as.POSIXct(format(df_elec$Date, "%Y-%m-%d %H:00:00"), tz = "Europe/Paris")

df_hourly <- df_elec %>%
  group_by(Hour) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()

# Nettoyage Météo
df_meteo <- df_meteo %>%
  mutate(Hour = parse_date_time(time, orders = c("ymd HMS", "ymd"), tz = "Europe/Paris")) %>%
  distinct(Hour, .keep_all = TRUE) %>%
  select(Hour, temperature_2m)

# Fusion et Interpolation
df_final <- data.frame(Hour = seq(min(df_hourly$Hour), max(df_hourly$Hour), by = "hour")) %>%
  left_join(df_hourly, by = "Hour") %>%
  left_join(df_meteo, by = "Hour") %>%
  filter(Hour >= as.POSIXct("2025-01-01 00:00:00", tz = "Europe/Paris")) %>%
  mutate(
    Value = as.numeric(na.interp(Value / 1000)), # MW -> GW
    temp  = as.numeric(na.interp(temperature_2m))
  )

#-------------------------------------------------------------------------------
# 3. FEATURE ENGINEERING (POUR TOUS LES MODÈLES)
#-------------------------------------------------------------------------------
# Variables de base
df_final$Froid <- pmax(0, 15 - df_final$temp)
df_final$Chaud <- pmax(0, df_final$temp - 25)

# Multi-Saisonnalité de Fourier
msts_total <- msts(df_final$Value, seasonal.periods = c(24, 168))
fourier_df <- as.data.frame(fourier(msts_total, K = c(4, 2)))

# Variables spécifiques XGBoost
df_xgb_data <- df_final %>%
  cbind(fourier_df) %>%
  arrange(Hour) %>%
  mutate(
    sin_hour = sin(2 * pi * hour(Hour) / 24),
    cos_hour = cos(2 * pi * hour(Hour) / 24),
    temp_roll3 = as.numeric(rollapplyr(temp, 3, mean, fill = "extend")),
    Froid_2 = Froid^2,
    Grand_Froid = pmax(0, 5 - temp),
    Lag24  = lag(Value, 24),
    Lag168 = lag(Value, 168),
    Hour_v = hour(Hour),
    Wday   = wday(Hour),
    IsWeekend = ifelse(Wday %in% c(1, 7), 1, 0)
  ) %>%
  filter(!is.na(Lag168)) # Nécessaire pour XGBoost

#-------------------------------------------------------------------------------
# 4. PRÉPARATION DES MATRICES ET OBJETS TS
#-------------------------------------------------------------------------------
h <- 24
n_total <- nrow(df_xgb_data)

# Pour les modèles Classiques (sur df_xgb_data pour que le test soit identique)
ts_train   <- ts(df_xgb_data$Value[1:(n_total - h)], frequency = 24)
ts_test    <- ts(df_xgb_data$Value[(n_total - h + 1):n_total], frequency = 24)
msts_train <- msts(df_xgb_data$Value[1:(n_total - h)], seasonal.periods = c(24, 168))
msts_test  <- msts(df_xgb_data$Value[(n_total - h + 1):n_total], seasonal.periods = c(24, 168))

# Regresseurs ARIMAX
xreg_train_arimax <- as.matrix(df_xgb_data[1:(n_total - h), c("Froid", "Chaud")])
xreg_test_arimax  <- as.matrix(df_xgb_data[(n_total - h + 1):n_total, c("Froid", "Chaud")])

# Fourier pour ARIMAX
f_train <- as.matrix(df_xgb_data[1:(n_total - h), colnames(fourier_df)])
f_test  <- as.matrix(df_xgb_data[(n_total - h + 1):n_total, colnames(fourier_df)])

# Matrice XGBoost
features_xgb <- c("temp", "temp_roll3", "Froid", "Froid_2", "Grand_Froid",
                  "sin_hour", "cos_hour", "Lag24", "Lag168",
                  "Hour_v", "Wday", "IsWeekend", colnames(fourier_df))
train_x <- as.matrix(df_xgb_data[1:(n_total - h), features_xgb])
test_x  <- as.matrix(df_xgb_data[(n_total - h + 1):n_total, features_xgb])
train_y <- df_xgb_data$Value[1:(n_total - h)]

#--------------------- 5. MODÉLISATION (1 à 11) --------------------------------
# Modèles 1 à 10 (Classiques)
tic(); fc1 <- meanf(ts_train, h = h); t <- toc(quiet = TRUE); model_times[["Historic Mean"]] <- t$toc - t$tic
tic(); fc2 <- naive(ts_train, h = h); t <- toc(quiet = TRUE); model_times[["Naive"]] <- t$toc - t$tic
tic(); fc3 <- snaive(ts_train, h = h); t <- toc(quiet = TRUE); model_times[["Seasonal Naive"]] <- t$toc - t$tic

tic(); fit4 <- auto.arima(ts_train, seasonal = FALSE, lambda = "auto")
fc4 <- forecast(fit4, h = h); t <- toc(quiet = TRUE); model_times[["Arima Auto"]] <- t$toc - t$tic

tic(); fit5 <- Arima(ts_train, order = c(0,0,3), seasonal = list(order = c(0,1,1), period = 24), lambda = "auto")
fc5 <- forecast(fit5, h = h); t <- toc(quiet = TRUE); model_times[["SARIMA Man"]] <- t$toc - t$tic

tic(); fit6 <- bats(ts_train, use.parallel = TRUE, num.cores = 6)
fc6 <- forecast(fit6, h = h); t <- toc(quiet = TRUE); model_times[["BATS"]] <- t$toc - t$tic

tic(); fit7 <- tbats(msts_train, use.parallel = TRUE, num.cores = 6, use.trend=FALSE)
fc7 <- forecast(fit7, h = h); t <- toc(quiet = TRUE); model_times[["TBATS"]] <- t$toc - t$tic

tic(); fit8 <- stlm(msts_train, method = "arima")
fc8 <- forecast(fit8, h = h); t <- toc(quiet = TRUE); model_times[["STLM"]] <- t$toc - t$tic

tic(); fit9 <- Arima(ts_train, xreg = xreg_train_arimax, order= c(5,0,0),
                     seasonal = list(order = c(1, 1, 2), period = 24), lambda = "auto")
fc9 <- forecast(fit9, xreg = xreg_test_arimax, h = h); t <- toc(quiet = TRUE); model_times[["ARIMAX (Temp)"]] <- t$toc - t$tic

tic(); fit10 <- Arima(ts_train, xreg = cbind(xreg_train_arimax, f_train), order = c(4, 1, 4), lambda = "auto")
fc10 <- forecast(fit10, xreg = cbind(xreg_test_arimax, f_test), h = h); t <- toc(quiet = TRUE); model_times[["ARIMAX (Full)"]] <- t$toc - t$tic

# Modèle 11 : XGBoost
tic()
fit_xgb <- xgboost(data = train_x, label = train_y, nrounds = 1000, max_depth = 6, eta = 0.05,
                   subsample = 0.7, colsample_bytree = 0.7, min_child_weight = 15,
                   objective = "reg:squarederror", verbosity = 0)
preds_xgb <- predict(fit_xgb, test_x)
fc11 <- list(mean = preds_xgb) # Emballage pour compatibilité
t <- toc(quiet = TRUE)
model_times[["XGBoost"]] <- t$toc - t$tic

#-------------------------------------------------------------------------------
# 6. COMPARAISON ET RÉSULTATS
#-------------------------------------------------------------------------------
model_list <- list(fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8, fc9, fc10, fc11)
names(model_list) <- c("Historic Mean", "Naive", "Seasonal Naive", "Arima Auto", "SARIMA Man",
                       "BATS", "TBATS", "STLM", "ARIMAX (Temp)", "ARIMAX (Full)", "XGBoost")

results <- lapply(names(model_list), function(nm) {
  actuals <- as.numeric(ts_test)
  preds   <- as.numeric(model_list[[nm]]$mean)
  acc     <- accuracy(preds, actuals)
  exec_time <- if(!is.null(model_times[[nm]])) model_times[[nm]] else NA

  data.frame(Model = nm, MAPE = acc[1, "MAPE"], RMSE = acc[1, "RMSE"], Time_Sec = round(exec_time, 3))
}) %>% bind_rows()

print(results %>% arrange(MAPE))

#-------------------------------------------------------------------------------
# 7. VISUALISATION FINALE (Comparaison Top Modèles)
#-------------------------------------------------------------------------------
ggplot() +
  geom_line(aes(x = df_xgb_data$Hour[(n_total-h+1):n_total], y = as.numeric(ts_test), color = "Réalité"), size = 1) +
  geom_line(aes(x = df_xgb_data$Hour[(n_total-h+1):n_total], y = as.numeric(fc11$mean), color = "XGBoost"), linetype = "dashed") +
  geom_line(aes(x = df_xgb_data$Hour[(n_total-h+1):n_total], y = as.numeric(fc10$mean), color = "ARIMAX"), linetype = "dotted") +
  labs(title = "Top 2 Modèles vs Réalité", y = "GW", x = "Heure") +
  theme_minimal()
