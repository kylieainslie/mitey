x <- 0:25

y <- unlist(lapply(x, f_norm,
                   weights = result6$wts,
                   mu = result6$mean,
                   sigma = result6$sd))

df <- data.frame(x = x, y = y)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_line()

p
