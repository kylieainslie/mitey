library(ggplot2)

setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/gold vs inflation")
df_gold<-read.csv("gold.csv")
df_inf_rate<-read.csv("inflation_rate_US.csv")

df_gold$Date<-as.Date(df_gold$Date)
df_inf_rate$observation_date<-as.Date(df_inf_rate$observation_date)

df_gold[,c("Open","High","Low","Volume","Currency")]<-list(NULL)

colnames(df_inf_rate)<-c("Date","Value")
colnames(df_gold)<-c("Date","Value")


df_gold$date_list<-1:nrow(df_gold)

f<- function (x){
  return (0.284 *x+230.359)
}



df_gold$y1 <-sapply(df_gold$date_list,f)
df_gold$y2 <-df_gold$Value-df_gold$y1
df_gold$y2
lm(df_gold$Value ~ df_gold$date_list)

graph_gold<-ggplot(df_gold,aes(x=date_list,y=Value))+
  geom_line()+
  ggtitle("Evolution du cours de l'or")
graph_gold








df_merge<-merge(df_gold,df_inf_rate,by="Date")
sapply(df_merge,class)

library(modelr)
library(tidyverse)
library(ggplot2)



df_inf_rate[1305,1]
df_gold[5703,1]


starting_date<-"2021-03-12"
ending_date<-"2022-09-02"

starting_row_gold<-which(df_gold==starting_date, arr.ind=TRUE)[1]
starting_row_inf<-which(df_inf_rate==starting_date, arr.ind=TRUE)[1]

ending_row_gold<-which(df_gold==ending_date, arr.ind=TRUE)[1]
ending_row_inf<-which(df_inf_rate==ending_date, arr.ind=TRUE)[1]

df_gold<-df_gold[starting_row_gold:ending_row_gold,]
df_inf_rate<-df_inf_rate[starting_row_inf:ending_row_inf,]


number_rows<-max(nrow(df_gold),nrow(df_inf_rate))

df_gold$type<-rep(c("cours de l'or"),nrow(df_gold))

df_inf_rate$type<-rep(c("inflation rate"),nrow(df_inf_rate))

df_rbind<-rbind(df_gold,df_inf_rate)

graph_gold<-ggplot(df_gold,aes(x=Date,y=Value))+
  geom_line()+
  ggtitle("Evolution du cours de l'or")
graph_gold
graph_inf_rate<-ggplot(df_inf_rate,aes(x=Date,y=Value))+
  geom_line()+
  ggtitle("Evolution du taux d'inflation")
graph_inf_rate




