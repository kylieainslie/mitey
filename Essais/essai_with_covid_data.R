file_name <-"ICC_covid.csv"
ICC<-read.csv(file_name)
df<-as.data.frame(ICC)
df<-abs(df)

ggplot(df, aes(x= ICC.interval)) +
  geom_histogram(
    stat = "count",
    fill = "lightblue",
    color = "black"
  ) +
  labs(
    title = "Histogramme des valeurs",
    x = "Valeurs",
    y = "Fréquence"
  ) +
  theme_minimal()

ICC2 <- as.vector(df)

result_8<-si_estim(ICC2$ICC.interval,n_routes=8)
plot_si_fit_result(result_8,ICC2$ICC.interval)

result_8
