library(ggplot2)
library(tseries)
library(forecast)


setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/gold vs inflation")
df_inflation<-read.csv("inflation_rate_US.csv")

df_inflation<-na.omit(df_inflation)
colnames(df_inflation)<-c("Date","Value")
df_inflation$Date<-as.Date(df_inflation$Date)

nrow(df_inflation)
df_inflation$Date[1]
df_inflation$Date[nrow(df_inflation)]
freq <- nrow(df_inflation)/5
freq

inflation_ts<-ts(df_inflation$Value,start=c(2021,1),frequency=freq)
plot(inflation_ts)


#take a training part and a test part :
ending_date_training<-"2025-01-01"
df_inflation$Date[994]
ending_row_training<-994

inflation_training_ts<-ts(df_inflation$Value[1:ending_row_training],start=c(2021,1),frequency=freq)
plot(inflation_training_ts)

inflation_test_ts<-ts(df_inflation$Value[ending_row_training+1:nrow(df_inflation)],start=c(2025,1),frequency=freq)
plot(inflation_test_ts)

inflation_training_ts

# decomposed components
ts_decomposed <- decompose(na.omit(inflation_training_ts))
plot(ts_decomposed)

inflation_training_ts

#Auto correlation
plot(inflation_training_ts)
acf(inflation_training_ts,lag.max=500)

# Forecasting
forecast <- forecast(inflation_training_ts, h = 31)

plot(forecast, main = "Forecasting for Next 12 Months")
lines(inflation_test_ts)




#Comparing forecast to real values
accuracy(forecast, inflation_test_ts)
