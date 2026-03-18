library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)

setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/gold vs inflation")
df_elec<-read.csv("elec_consumption.csv",sep=";")

col<-colnames(df_elec)
col[3:6]
df_elec[,col[1:2]]<-list(NULL)
df_elec[,col[7:37]]<-list(NULL)


colnames(df_elec)
sapply(df_elec,class)
df_elec$Date.et.Heure[1]
df_elec$Date_only <- substr(df_elec$Date.et.Heure, 1, 10)
df_elec[,col[1:5]]<-list(NULL)


daily_sum <- aggregate(Consommation..MW. ~ Date_only,
                       data = df_elec,
                       FUN = sum,
                       na.rm = TRUE)
nrow(daily_sum)
colnames(daily_sum)<-c("Date","Value")

freq<-7
conso_ts<-ts(daily_sum$Value,frequency=freq)
plot(conso_ts)



#take a training part and a test part :
ending_date_training<-"2022-01-01"
daily_sum$Date[3654]
ending_row_training<-3654

conso_training_ts<-ts(daily_sum$Value[1:ending_row_training],frequency=freq)
plot(conso_training_ts)

conso_test_ts<-ts(daily_sum$Value[ending_row_training+1:nrow(daily_sum)],frequency=freq)
plot(conso_test_ts)


adf_value<-adf.test(conso_training_ts)
adf_value
#p-value <0.05 ==> data is stationary (not for seasonality), no need for differenciation ==> d=0


Acf(conso_training_ts)
Pacf(conso_training_ts)


fit_auto <- auto.arima(diff_ts, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
summary(fit_auto)
checkresiduals(fit_auto)

#essai avec qq valeurs


frequency(conso_training_ts)
length(conso_training_ts)
sum(is.na(conso_training_ts))
