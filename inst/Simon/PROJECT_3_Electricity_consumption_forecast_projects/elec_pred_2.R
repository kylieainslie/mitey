library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)

setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/PROJECT_3_Electricity_consumption_forecast_projects")

df_elec <- read.csv("consommation-quotidienne-brute.csv", sep = ";")

# Garder uniquement date + valeur
df_elec <- df_elec[, c(1, 9)]
colnames(df_elec) <- c("Date", "Value")

# Nettoyage du timestamp
df_elec$Date <- substr(df_elec$Date, 1, 16)
df_elec$Date <- gsub("T", " ", df_elec$Date)
df_elec$Date <- as.POSIXct(df_elec$Date, format = "%Y-%m-%d %H:%M")


# Supprimer NA
df_elec <- na.omit(df_elec)

# Trier chronologiquement
df_elec <- df_elec %>% arrange(Date)

# Arrondir à l'heure inférieure (ex: 14:30 -> 14:00)
df_elec$Hour <- as.POSIXct(format(df_elec$Date, "%Y-%m-%d %H:00:00"))

# Agrégation horaire
df_hourly <- df_elec %>%
  group_by(Hour) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()

delta_hour <- diff(df_hourly$Hour)
table(delta_hour)


full_hours <- data.frame(
  Hour = seq(from = min(df_hourly$Hour),
             to   = max(df_hourly$Hour),
             by   = "hour")
)

df_hourly_full <- full_hours %>%
  left_join(df_hourly, by = "Hour")

sum(is.na(df_hourly_full$Value))


freq<-24
ts_hourly<-ts(df_hourly_full$Value,frequency=freq)
ts_hourly <- na.interp(ts_hourly)
plot(ts_hourly, main = "Hourly electricity series")


decomp_hourly <- stl(ts_hourly, s.window = "periodic")
plot(decomp_hourly)


msts_hourly <- msts(as.numeric(ts_hourly), seasonal.periods = c(24, 168))
decomp_multi <- mstl(msts_hourly)
plot(decomp_multi)

h<-24
n_ts <- length(ts_hourly)

ts_train <- window(ts_hourly, end = time(ts_hourly)[n_ts - h])
ts_test  <- window(ts_hourly, start = time(ts_hourly)[n_ts - h + 1])



#----------MODEL 1 : HISTORIC MEAN---------------
fc1<-meanf(ts_train,h)
acc1<- accuracy(fc1,ts_test)
plot(fc1$mean,main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)

#----------MODEL 2 : NAIVE REPEATING LAST OBS----
fc2<-naive(ts_train,h)
acc2<- accuracy(fc2,ts_test)
plot(fc2$mean, main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)

#----------MODEL 3 : NAIVE REPEATING THE LAST VALUES OVER THE LAST DAY----
fc3<-snaive(ts_train,h)
acc3<- accuracy(fc3,ts_test)
plot(fc3$mean,main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)

#----------MODEL 4 : ARIMA MODEL WITHOUT SEASONING - AUTO ---
fit_arima <- auto.arima(ts_train, seasonal=FALSE)
fc4<-forecast(fit_arima, h)
acc4<- accuracy(fc4,ts_test)
plot(fc4$mean,main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)




#----------MODEL 5 : ARIMA MODEL WITH SEASONING  - MANUAL --

#after tests, the best model is : ARIMA(0,0,3)(0,1,1)
order=c(0,0,3)
seasonal=c(0,1,1)

fit_arima2 <- Arima(
  ts_train,
  order = c(0,0,3),
  seasonal = list(order = c(0,1,1), period = 24)
)
fc5<-forecast(fit_arima2, h)
acc5<- accuracy(fc5,ts_test)
plot(fc5$mean,main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)


#----------MODEL 6 : BATS MODEL----

fit_bats <- bats(ts_train)
fc6  <- forecast(fit_bats, h)
acc6 <- accuracy(fc6, ts_test)
plot(fc6$mean,main = "forecast vs reality over the last day")
lines(ts_test, col = "red")
legend("bottom", legend = c("Forecast", "Test"),
       col = c("black", "red"), lty = 1)


#-----------------MODELE 7 TBATS----------------------------------------
# Créer msts sur toutes les données
df_hourly_full$Value<-na.interp(df_hourly_full$Value)
msts_elec <- msts(df_hourly_full$Value, seasonal.periods = c(24,7*24))

# Définir horizon
h <- 24  # horizon = 24h

# Dernier index temporel du train
end_time_train <- time(msts_elec)[length(msts_elec) - h]

# Premier index temporel du test
start_time_test <- time(msts_elec)[length(msts_elec) - h + 1]

# Subset avec window()
msts_train <- window(msts_elec, end = end_time_train)
msts_test  <- window(msts_elec, start = start_time_test)

# Fit TBATS
fit_tbats <- tbats(msts_train/1000)

# Forecast
fc7 <- forecast(fit_tbats, h = h)
fc7$mean <- fc7$mean * 1000
fc7$lower <- fc7$lower * 1000
fc7$upper <- fc7$upper * 1000

# Comparer sur les mêmes types d'objets
acc7 <- accuracy(as.numeric(fc7$mean), as.numeric(msts_test))

# Plot
plot(fc7$mean, main = "TBATS forecast vs reality")
lines(as.numeric(msts_test), col = "red")
legend("bottom", legend = c("Forecast", "Test"), col = c("black","red"), lty = 1)






#MODELE 8 STLM

fit_stlm <- stlm(msts_train / 1000, method = "arima")
fc8 <- forecast(fit_stlm, h = 24)
fc8$mean <- fc8$mean * 1000
fc8$lower <- fc8$lower * 1000
fc8$upper <- fc8$upper * 1000


acc8 <- accuracy(as.numeric(fc8$mean), as.numeric(msts_test))

# Plot
plot(fc8$mean, main = "TBATS forecast vs reality")
lines(as.numeric(ts_test), col = "red")
legend("bottom", legend = c("Forecast", "Test"), col = c("black","red"), lty = 1)



results <- data.frame(
  Model = c("Historic Mean", "Naive", "Seasonal Naive","Arima without seasoning","Arima with seasoning","BATS","TBATS","STLM"),
  ME    = c(acc1["Test set", "ME"], acc2["Test set", "ME"], acc3["Test set", "ME"],acc4["Test set", "ME"],acc5["Test set", "ME"],acc6["Test set", "ME"],acc7["Test set", "ME"],acc8["Test set", "ME"]),
  RMSE  = c(acc1["Test set", "RMSE"], acc2["Test set", "RMSE"], acc3["Test set", "RMSE"],acc4["Test set", "RMSE"],acc5["Test set", "RMSE"],acc6["Test set", "RMSE"],acc7["Test set", "RMSE"],acc8["Test set", "RMSE"])
)
print(results)


