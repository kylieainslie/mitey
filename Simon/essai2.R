setwd("C:/Users/carre/OneDrive/Bureau/Mitey_project/mitey/Simon")

data<-read.csv("adult_train.csv")
data[data==""]<-"NA"
colors<-c("Male"="red", "Female"="blue")
workclass<-unique(data$Workclass)
symbol <-c(" State-gov"=0, " Self-emp-not-inc"=1," Private"=2," Federal-gov"=3,
" Local-gov"=4,"NA"=5," Self-emp-inc"=6," Without-pay"=7," Never-worked"=8   )


data$Sex <- trimws(data$Sex)       # supprime les espaces
data$Sex <- factor(data$Sex)       # optionnel mais sécurise les niveaux
#model<-lm(Hours_per_week  ~ Education_Num, data=data)
#summary(model)

#plot(data$Age,
#     data$Capital_Gain,
#     col=colors[data$Sex],
#     pch = symbol[data$Workclass],
#     main = "Capital gain as a function of Age",
#     xlab = "Education",
#     ylab = "Capital Gain")

colnames(data)

#model<-lm(Capital_Gain  ~ Education_Num+Hours_per_week, data=data)
#summary(model)

hist(data$Education_Num, data=data, breaks=20)
