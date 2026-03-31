library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)
library(prophet)
library(lubridate)
library(readr)
#----------------------Data recuperation----------------------------------------
setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/inst/Simon/PROJECT_3_Electricity_consumption_forecast_projects")
df_elec <- read.csv("consommation-quotidienne-brute.csv", sep = ";")
df_meteo <- read.csv("meteo_paris_2025.csv",sep=",")
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



#--------------------Adding Weather to the df--------------
df_meteo <- df_meteo %>%
  mutate(
    Hour = parse_date_time(
      time,
      orders = c("ymd HMS", "ymd"),
      tz = "Europe/Paris"
    )
  ) %>%
  select(-time)

# Vérification
sum(is.na(df_meteo$Hour))
head(df_meteo$Hour, 30)



df_hourly_full$Hour <- as.POSIXct(df_hourly_full$Hour, tz = "Europe/Paris")
df_meteo$Hour       <- as.POSIXct(df_meteo$Hour, tz = "Europe/Paris")



df_meteo <- df_meteo %>% distinct(Hour, .keep_all = TRUE)


df_final <- df_hourly_full %>%
  left_join(df_meteo, by = "Hour")

df_final <- df_final %>%
  filter(Hour >= as.POSIXct("2025-01-01 00:00:00", tz = "Europe/Paris"))



#créer variables température pertinentes
# On considère qu'on commence à chauffer en dessous de 15°C
df_final$Froid <- pmax(0, 15 - df_final$temperature_2m)

ts_hourly <- ts(df_final$Value, frequency = 24)
xreg_temp <- as.matrix(df_final$Froid)
sum(is.na(xreg_temp))
xreg_temp <- na.interp(xreg_temp)
n_ts <- length(ts_hourly)

ts_train <- window(ts_hourly, end = time(ts_hourly)[n_ts - h])
ts_test  <- window(ts_hourly, start = time(ts_hourly)[n_ts - h + 1])

xreg_train <- xreg_temp[1:(n_ts - h), drop = FALSE]
xreg_test  <- xreg_temp[(n_ts - h + 1):n, drop = FALSE]

fit_arimax <- auto.arima(
  ts_train,
  xreg = xreg_train,
  seasonal = TRUE,trace=TRUE,lambda="auto",stepwise = TRUE
)

fc_arimax <- forecast(
  fit_arimax,
  xreg = xreg_test,
  h = h
)


plot(fc_arimax)
accuracy(fc_arimax$mean, ts_test)