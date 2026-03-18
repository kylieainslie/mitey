#importing library
library(ggplot2)
library(tseries)
library(forecast)
library(dplyr)

#importing data
data("AirPassengers")
ts_data <- AirPassengers
plot(ts_data)

#splitting data
ts_data_train<-subset(ts_data,end = length(ts_data) - 12)
ts_data_test<-subset(ts_data,start = length(ts_data) - 11)

#decomposition
fit <- stl(ts_data_train,s.window=10)
plot(fit)

#First model : Holt Winters
fit1 <- HoltWinters(ts_data_train)
fc1<-forecast(fit1,12)
plot(fc1)


acc_model_1<-accuracy(fc1,ts_data_test)

#Second Model : ARIMA
#checking for dependency to lagged values
Acf(ts_data)
Pacf(ts_data)


fit2<-auto.arima(ts_data_train)
summary(fit2)

fc2<-forecast(fit2,h=length(ts_data_test))
plot(fc2)
acc_model_2<-accuracy(fc2,ts_data_test)


#Third Model : ETS
fit3 <- ets(ts_data_train)
fc3<-forecast(fit3,12)
plot(fc3)

acc_model_3<-accuracy(fc3,ts_data_test)


