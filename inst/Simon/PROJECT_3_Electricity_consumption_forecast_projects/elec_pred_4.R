library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)
library(prophet)
library(lubridate)
library(readr)
library(tictoc)
model_times <- list()

#---------------------- 1. Récupération des données ----------------------------
setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/inst/Simon/PROJECT_3_Electricity_consumption_forecast_projects")

df_elec <- read.csv("consommation-quotidienne-brute.csv", sep = ";")
df_meteo <- read.csv("meteo_paris_2025.csv", sep = ",")

#--------------------- 2. Nettoyage Consommation --------------------------------
df_elec <- df_elec[, c(1, 9)]
colnames(df_elec) <- c("Date", "Value")

df_elec$Date <- substr(df_elec$Date, 1, 16)
df_elec$Date <- gsub("T", " ", df_elec$Date)
df_elec$Date <- as.POSIXct(df_elec$Date, format = "%Y-%m-%d %H:%M", tz = "Europe/Paris")

df_elec <- na.omit(df_elec) %>% arrange(Date)
df_elec$Hour <- as.POSIXct(format(df_elec$Date, "%Y-%m-%d %H:00:00"), tz = "Europe/Paris")

df_hourly <- df_elec %>%
  group_by(Hour) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()

# Création du calendrier complet pour éviter les trous
full_hours <- data.frame(
  Hour = seq(from = min(df_hourly$Hour), to = max(df_hourly$Hour), by = "hour")
)

df_hourly_full <- full_hours %>% left_join(df_hourly, by = "Hour")

#--------------------- 3. Nettoyage Météo & Fusion -----------------------------
df_meteo <- df_meteo %>%
  mutate(Hour = parse_date_time(time, orders = c("ymd HMS", "ymd"), tz = "Europe/Paris")) %>%
  distinct(Hour, .keep_all = TRUE) %>%
  select(Hour, temperature_2m)

df_final <- df_hourly_full %>%
  left_join(df_meteo, by = "Hour") %>%
  # On filtre à partir de 2025 comme dans ton second code
  filter(Hour >= as.POSIXct("2025-01-01 00:00:00", tz = "Europe/Paris")) %>%
  # Scaling immédiat (MW -> GW)
  mutate(Value = Value / 1000)

# Création de la variable "Froid" (Heating Degree Days)
df_final$Froid <- pmax(0, 15 - df_final$temperature_2m)
df_final$Chaud <- pmax(0, df_final$temperature_2m-25)

# Remplissage des éventuels NA restants (interpolation)
df_final$Value <- na.interp(df_final$Value)
df_final$Froid <- na.interp(df_final$Froid)
df_final$Chaud <- na.interp(df_final$Chaud)

#--------------------- 4. Création des objets Time Series ----------------------
h <- 24
n_ts <- nrow(df_final)

# Objets TS classiques (Fréquence 24)
ts_total <- ts(df_final$Value, frequency = 24)
ts_train <- window(ts_total, end = time(ts_total)[n_ts - h])
ts_test  <- window(ts_total, start = time(ts_total)[n_ts - h + 1])

# Objets MSTS (Multi-saisonnalité 24h et 7j)
msts_total <- msts(df_final$Value, seasonal.periods = c(24, 24*7))
msts_train <- window(msts_total, end = time(msts_total)[n_ts - h])
msts_test  <- window(msts_total, start = time(msts_total)[n_ts - h + 1])


# Objets Regresseurs (Température : Froid ET Chaud)
xreg_total <- as.matrix(df_final[, c("Froid", "Chaud")]) # On prend les deux colonnes
xreg_train <- xreg_total[1:(n_ts - h), , drop = FALSE]
xreg_test  <- xreg_total[(n_ts - h + 1):n_ts, , drop = FALSE]


# Utilisation de msts_train qui contient les fréquences 24 et 168 (24*7)
# K = c(10, 5) signifie 10 paires pour la journée et 5 pour la semaine
fourier_train <- fourier(msts_train, K = c(4, 2))
fourier_test  <- fourier(msts_train, K = c(4, 2), h = h)

# Le reste de ton code ARIMAX reste identique
xreg_train_full <- cbind(xreg_train, fourier_train)
xreg_test_full  <- cbind(xreg_test, fourier_test)



#--------------------- 5. Modélisation (1 à 9) ---------------------------------
tic()
# 1-3 : Modèles basiques
fc1 <- meanf(ts_train, h = h)
t <- toc(quiet = TRUE)
model_times[["Historic Mean"]] <- t$toc - t$tic

tic()
fc2 <- naive(ts_train, h = h)
t <- toc(quiet = TRUE)
model_times[["Naive"]] <- t$toc - t$tic

tic()
fc3 <- snaive(ts_train, h = h)
t <- toc(quiet = TRUE)
model_times[["Seasonal Naive"]] <- t$toc - t$tic

# 4 : ARIMA Auto (Sans saisonnalité)
tic()
fit4 <- auto.arima(ts_train, seasonal = FALSE, lambda = "auto")
fc4  <- forecast(fit4, h = h)
t <- toc(quiet = TRUE)
model_times[["Arima Auto"]] <- t$toc - t$tic

# 5 : SARIMA Manuel
tic()
fit5 <- Arima(ts_train, order = c(0,0,3), seasonal = list(order = c(0,1,1), period = 24), lambda = "auto")
fc5  <- forecast(fit5, h = h)
t <- toc(quiet = TRUE)
model_times[["SARIMA Man"]] <- t$toc - t$tic

# 6 : BATS
tic("fit modèle 6")
fit6 <- bats(ts_train, use.parallel = TRUE, num.cores =6)
fc6  <- forecast(fit6, h = h)
t <- toc(quiet = TRUE)
model_times[["BATS"]] <- t$toc - t$tic

# 7 : TBATS (MSTS)
tic()
fit7 <- tbats(msts_train,use.parallel = TRUE, num.cores = 6, use.trend=FALSE)
fc7  <- forecast(fit7, h = h)
t <- toc(quiet = FALSE)
model_times[["TBATS"]] <- t$toc - t$tic

# 8 : STLM (MSTS)
tic()
fit8 <- stlm(msts_train, method = "arima")
fc8  <- forecast(fit8, h = h)
t <- toc(quiet = TRUE)
model_times[["STLM"]] <- t$toc - t$tic

# 9 : ARIMAX 1 (ARIMA + Température)
tic()
fit9 <- Arima(ts_train, xreg = xreg_train, order= c(5,0,0),seasonal = list(order = c(1, 1, 2), period = 24),lambda = "auto")
fc9  <- forecast(fit9, xreg = xreg_test, h = h)
t <- toc(quiet = TRUE)
model_times[["ARIMAX (Temp)"]] <- t$toc - t$tic

# 10 : ARIMAX 2 (ARIMA + Température + décomposition fourier)
tic()
fit10 <- Arima(ts_train, xreg = xreg_train_full, order = c(4, 1, 4),lambda = "auto")
fc10  <- forecast(fit10, xreg = xreg_test_full, h = h)
acc10<-accuracy(fc10$mean,ts_test)
acc10
t <- toc(quiet = TRUE)
model_times[["ARIMAX (Temp + Fourier)"]] <- t$toc - t$tic

#--------------------- 6. Comparaison des Résultats ----------------------------

# Liste des modèles pour automatiser
model_list <- list(fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8, fc9, fc10)
names(model_list) <- c("Historic Mean", "Naive", "Seasonal Naive", "Arima Auto",
                       "SARIMA Man", "BATS", "TBATS", "STLM", "ARIMAX (Temp)", "ARIMAX (Temp + Fourier)")

# Calcul robuste des métriques (on compare les vecteurs numériques)
results <- lapply(names(model_list), function(nm) {
  # Pour TBATS/STLM on compare à msts_test, pour les autres à ts_test
  actuals <- if(nm %in% c("TBATS", "STLM")) as.numeric(msts_test) else as.numeric(ts_test)
  preds <- as.numeric(model_list[[nm]]$mean)

  acc <- accuracy(preds, actuals)

  exec_time <- if(!is.null(model_times[[nm]])) model_times[[nm]] else NA

  data.frame(Model = nm, MAPE = acc[1, "MAPE"], RMSE = acc[1, "RMSE"],Time_Sec = round(exec_time, 3))
}) %>% bind_rows()


print(results %>% arrange(MAPE))

#--------------------- 7. Visualisation (Exemple ARIMAX) -----------------------

# On utilise le zoom de 3 jours comme précédemment
autoplot(fc10, include = 24*3) +
  autolayer(ts_test, series = "Réalité", color = "red") +
  ggtitle("") +
  theme_minimal()

checkresiduals(fit10)
