library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)
library(prophet)
library(lubridate)
#----------------------Data recuperation----------------------------------------
setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/PROJECT_3_Electricity_consumption_forecast_projects")
df_elec <- read.csv("consommation-quotidienne-brute.csv", sep = ";")

#---------------------Data cleaning---------------------------------------------
#removing useless columns and renaming them
df_elec <- df_elec[, c(1, 9)]
colnames(df_elec) <- c("Date", "Value")

# Timestamp cleaning
df_elec$Date <- substr(df_elec$Date, 1, 16)
df_elec$Date <- gsub("T", " ", df_elec$Date)
df_elec$Date <- as.POSIXct(df_elec$Date, format = "%Y-%m-%d %H:%M")


# NA cleaning
df_elec <- na.omit(df_elec)

# Chronological ordering
df_elec <- df_elec %>% arrange(Date)

# Rounding hour
df_elec$Hour <- as.POSIXct(format(df_elec$Date, "%Y-%m-%d %H:00:00"))

# Agregation by hour
df_hourly <- df_elec %>%
  group_by(Hour) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()

#creating robust hour column
full_hours <- data.frame(
  Hour = seq(from = min(df_hourly$Hour),
             to   = max(df_hourly$Hour),
             by   = "hour")
)

df_hourly_full <- full_hours %>%
  left_join(df_hourly, by = "Hour")


#SCALING
scale_factor <-1000
df_hourly_full$Value <-df_hourly_full$Value/scale_factor



#-------------- Use clean dataframe to create time serie object-----------------
freq<-24
ts_hourly<-ts(df_hourly_full$Value,frequency=freq) %>% na.interp()
plot(ts_hourly, main = "Hourly electricity series")


decomp_hourly <- stl(ts_hourly, s.window = "periodic")
plot(decomp_hourly)


msts_hourly <- msts(as.numeric(ts_hourly), seasonal.periods = c(24, 24*7))

h<-24
n_ts <- length(ts_hourly)

#---------------------Subsetting TS---------------------------------------------
ts_train <- window(ts_hourly, end = time(ts_hourly)[n_ts - h])
ts_test  <- window(ts_hourly, start = time(ts_hourly)[n_ts - h + 1])

# Dernier index temporel du train
end_time_train <- time(msts_hourly)[length(msts_hourly) - h]

# Premier index temporel du test
start_time_test <- time(msts_hourly)[length(msts_hourly) - h + 1]

# Subset avec window()
msts_train <- window(msts_hourly, end = end_time_train)
msts_test  <- window(msts_hourly, start = start_time_test)

#----------MODEL 1 : HISTORIC MEAN----------------------------------------------
fc1<-meanf(ts_train,h)
acc1<- accuracy(fc1,ts_test)


#----------MODEL 2 : NAIVE REPEATING LAST OBS-----------------------------------
fc2<-naive(ts_train,h)
acc2<- accuracy(fc2,ts_test)


#----------MODEL 3 : NAIVE REPEATING THE LAST VALUES OVER THE LAST DAY----------
fc3<-snaive(ts_train,24)
acc3<- accuracy(fc3,ts_test)


#----------MODEL 4 : ARIMA MODEL WITHOUT SEASONING - AUTO ----------------------
fit_arima <- auto.arima(ts_train, seasonal=FALSE, lambda="auto")
fc4<-forecast(fit_arima, h)
acc4<- accuracy(fc4,ts_test)


#----------MODEL 5 : ARIMA MODEL WITH SEASONING  - MANUAL ----------------

#after tests, the best model is : ARIMA(0,0,3)(0,1,1)
order=c(0,0,3)
seasonal=c(0,1,1)

fit_arima2 <- Arima(
  ts_train,
  order = c(0,0,3),
  seasonal = list(order = c(0,1,1), period = 24),
  lambda="auto"
)
fc5<-forecast(fit_arima2, h)
acc5<- accuracy(fc5,ts_test)


#----------MODEL 6 : BATS MODEL----------------------------------------
fit_bats <- bats(ts_train)
fc6  <- forecast(fit_bats, h)
acc6 <- accuracy(fc6, ts_test)


#-----------------------Ploting models 1 to 6 results ------------------------

fc_list <- list(fc1, fc2, fc3, fc4, fc5, fc6)
noms    <- c("Mean", "Naive", "S.Naive", "ARIMA_Auto", "SARIMA_Manual", "BATS")

number_of_days<-3
for(i in 1:length(fc_list)) {
  p <- autoplot(fc_list[[i]], include = 24*number_of_days) + # Affiche les 3 derniers jours
    autolayer(ts_test, series = "Réalité", color = "red") +

    ggtitle(paste("Modèle", i, ":", noms[i])) +
    theme_minimal()

  print(p) # Affiche le graph dans l'onglet Plots
}







#-----------------MODELE 7 TBATS----------------------------------------
fit_tbats <- tbats(msts_train)
fc7 <- forecast(fit_tbats, h = h)
acc7 <- accuracy(fc7, msts_test)

autoplot(fc7, include = 24*number_of_days) +
  autolayer(msts_test, series = "Réalité", color = "red", lwd = 1) +
  labs(
    title = "Modèle 7 : TBATS vs Réalité",
    subtitle = "Contexte : 1 semaine d'historique affichée",
    x = "Temps (Unités msts)",
    y = "Consommation (GW)",
    colour = "Légende"
  ) +
  theme_minimal()


#-------------------MODELE 8 STLM--------------------------------------------------

fit_stlm <- stlm(msts_train, method = "arima")
fc8 <- forecast(fit_stlm, h = 24)
acc8 <- accuracy(fc8, msts_test)

autoplot(fc8, include = 24*number_of_days) +
  autolayer(msts_test, series = "Réalité", color = "red", lwd = 1) +
  labs(
    title = "Modèle 8 : STLM vs Réalité",
    subtitle = "Méthode de décomposition STL + ARIMA",
    x = "Temps (Unités msts)",
    y = "Consommation (GW)"
  ) +
  theme_minimal()


#-------------------------------show results and compare models--------------------
model_list <- list(fc1, fc2, fc3, fc4, fc5, fc6, fc7, fc8)
names(model_list) <- c("Historic Mean", "Naive", "Seasonal Naive", "Arima without seasoning", "Arima with seasoning", "BATS", "TBATS", "STLM")


# Création du dataframe final proprement


results <- data.frame(
  Model = names(model_list),
  MAPE  = c(acc1[2, "MAPE"], acc2[2, "MAPE"], acc3[2, "MAPE"], acc4[2, "MAPE"],
            acc5[2, "MAPE"], acc6[2, "MAPE"], acc7[2, "MAPE"], acc8[2, "MAPE"]),
  RMSE  = c(acc1[2, "RMSE"], acc2[2, "RMSE"], acc3[2, "RMSE"], acc4[2, "RMSE"],
            acc5[2, "RMSE"], acc6[2, "RMSE"], acc7[2, "RMSE"], acc8[2, "RMSE"])
)

print(results)


