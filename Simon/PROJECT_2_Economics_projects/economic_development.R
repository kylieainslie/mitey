library(tidyr)
library(zoo)
library(ggplot2)
library(dplyr)

setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon/PROJECT_2_Economics_projects")
data<-read.csv("economic_dataset.csv")
#récupérer le csv

df<-as.data.frame(data)
#le transformer en dataframe

df[df==".."]<-"NA"
#identifier les valeurs manquantes

df$df_emptyness=rowSums(df=="NA")
len<-ncol(data)
df$df_completeness<-100 - df$df_emptyness/len*100

df<-subset(df, df_completeness>10)
#créer un df qui supprime les lignes avec trop de valeurs manquantes
df$df_completeness<-NULL
df$df_emptyness<-NULL
df<-df[c(-2,-3,-4)]

colnames(df)[2:ncol(df)] <- substr(colnames(df)[2:ncol(df)], 2, 7)

df_long <- pivot_longer(
  df,
  cols = -Country,
  names_to = "Quarter",
  values_to = "Value"
)

df_long$Quarter<-as.Date(as.yearqtr(df_long$Quarter))
df_long$Value<-as.numeric(df_long$Value)


graph<-ggplot(df_long, aes(x=Quarter, y =Value, group=Country,color=Country))+
  geom_line(show.legend = FALSE)+
  ggtitle("Country's GDP in millions USD")
print(graph)

last_date<-max(df_long$Quarter)
last_date

quarters<-unique(df_long$Quarter)
quarters


last_record<-df_long %>% filter(Quarter=="2025-04-01")
last_record %>%
  arrange(desc(Value)) %>%

slice(1:10)
