# analysis of weekly scabies data from NIVEL
# data courtesy of Mariette Hooiveld

# load required packages
library(readxl)
library(tidyverse)
library(astsa)
install.packages('EpiEstim', repos = c('https://mrc-ide.r-universe.dev',
                                       'https://cloud.r-project.org'))
library(EpiEstim) # it needs to be version 2.4

# read in data -----------------------------------------------------------------
setwd("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies/data")
nivel_wkly_data <- read_xlsx("scabies_data_weekly.xlsx") %>%
  # fix/translate variable names
  rename(diagnosis_code = `Diagnose (ICPC)`,
         year = `ISO-jaar`,
         week_num = `ISO-weeknr (ma-zo)`,
         pop_size = `Aantal populatie`,
         cases = `Aantal prevalente cases`,
         prev_per_100000 = `Prevalentie per 100.000`) %>%
  # drop diagnosis var
  select(-diagnosis_code) %>%
  # create new var that combines year and week
  mutate(yr_wk = paste(year, week_num, sep = "_"),
         year = as.factor(year))

# plot frequencies over time ---------------------------------------------------
yrs <- unique(nivel_wkly_data$year)
x_breaks <- paste(2011:2024, c(rep(1,14)), sep = "_")
freq_plot <- ggplot(data = nivel_wkly_data,
                    aes(x = yr_wk, y = prev_per_100000, fill = year)) +
  geom_bar(stat = "identity") +
  labs(y = "Cases per 100,000", x = "Year") +
  scale_fill_viridis_d() +
  scale_x_discrete(breaks = x_breaks, labels = yrs) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    panel.background = element_blank(),
    axis.text.x = element_text(hjust = 0),
    #axis.ticks = element_blank()
  )
freq_plot

# lets plot as a line plot
line_plot <- ggplot(data = nivel_wkly_data %>%
                      mutate(time = seq(from = 1, to = nrow(nivel_wkly_data),
                                        by = 1)),
                    aes(x = time, y = prev_per_100000)) +
  geom_line()
line_plot


# spectral analysis ------------------------------------------------------------

# time series
y <- nivel_wkly_data$prev_per_100000

# raw periodogram
mvspec(y, log = "no")

# smoothed periodogram using a Daniell kernel with m = 2
mvspec(y, span = 5, log = "no")

# smoothed periodogram using two passes of a Daniell kernel with m = 2 on each pass
scabies_spec <- mvspec(y, spans = c(5,5), log = "no")
# find where the peaks are located
details_df <- as.data.frame(scabies_spec$details)
# I think we can ignore this because it is the first point, so the value of freq,
# period, and spectrum at the next highest value is below (respectively)
# [14,]    0.0194  51.4286  557.0444
# so the peak occurs at 51.4 weeks
# another peak occurs at 26.67 weeks
# [27,]    0.0375  26.6667   58.9336

# add lines at the peak
abline(v=1/51.4, lty="dotted")
abline(v=1/26.6667, lty = "dotted")

# estimate Rt ------------------------------------------------------------------
# we will first approximate daily incidence using EpiEstim from our weekly
# aggregated incidence data
# then Rt will be estimated using the approximated daily incidence
rt <- estimate_R(incid = nivel_wkly_data$cases,
                 dt = 7L,
                 dt_out = 7L,
                 recon_opt = "naive",
                 iter = 10L,
                 tol = 1e-6,
                 grid = list(precision = 0.001, min = -1, max = 1),
                 method = "parametric_si",
                 config = make_config(list(
                   mean_si = 95.57,
                   std_si = 15.17
                   )))

# plot output
plot(rt) # full output

plot(rt, "incid") # Reconstructed daily incidence only
plot(rt, "R") # Rt estimates only
plot(rt, "SI") # SI estimates only

# better plot
r_dat <- rt$R %>%
  mutate(date = as.Date(t_start, origin = "2011-01-01"))

# Breaks for background rectangles
start_breaks <- c(as.Date("2011-01-01"),
                  seq(as.Date("2011-03-01"), tail(r_dat$date, 1), by = "quarter"))
end_breaks <- c(start_breaks[-1] - 1, tail(r_dat$date,1))
rects <- data.frame(xstart = start_breaks,
                    xend = end_breaks,
                    season = factor(c(rep(c("winter", "spring", "summer", "autumn"),
                                   length.out = length(start_breaks))),
                                   levels = c("winter", "spring", "summer", "autumn"))
                    )

fancy_plot <- ggplot() +
  geom_rect(data = rects, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,
                              fill = season), alpha = 0.4) +
  scale_fill_viridis_d() +
  geom_line(data = r_dat, aes(x = date, y = `Mean(R)`)) +
  # geom_ribbon(data = r_dat, aes(x = date, ymin = `Quantile.0.025(R)`,
  # ymax = `Quantile.0.975(R)`)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme(
    panel.background = element_blank() #,
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )
fancy_plot
# Evaluate growth rate ---------------------------------------------------------
nivel_wkly_data2 <- nivel_wkly_data %>%
  mutate(week_tot = row_number())

# Define the exponential growth model
exponential_model <- function(x, a, b) {
  a * exp(b * x)
}

# get starting values
# Fit linear regression to log-transformed data
lm_fit <- lm(log(cases) ~ week_tot, data = nivel_wkly_data2)

# Extract coefficients
a_init <- exp(coef(lm_fit)[1])
b_init <- coef(lm_fit)[2]

# Fitting the model with initial values
fit <- nls(cases ~ exponential_model(week_tot, a, b),
           data = nivel_wkly_data2,
           start = list(a = a_init, b = b_init),
           control = list(minFactor = 0.00001, maxiter = 100000))

# Summary of the model
summary(fit)

# Plot growth rate
plot(nivel_wkly_data2$week_tot, nivel_wkly_data2$cases,
     ylim = c(0, max(nivel_wkly_data2$cases) * 1.1),
     xlim = c(min(nivel_wkly_data2$week_tot), max(nivel_wkly_data2$week_tot)),
     xlab = "Year", ylab = "Cases per 1000", main = "Scabies Exponential Growth Model Fit")
lines(nivel_wkly_data2$week_tot, predict(fit), col = "red")
legend("topleft", legend = c("Data", "Model Fit"), col = c("black", "red"), lty = 1)

# add line to fancy plot from above
df_gr_plot <- nivel_wkly_data2 %>%
  select(year, cases, week_tot) %>%
  mutate(model_fit = predict(fit),
         date = as.Date(7*(week_tot - 1), origin = "2011-01-01"))

# gr_plot <- ggplot(data = df_gr_plot, aes(x = date)) +
#   geom_point(aes(y = cases, alpha = 0.7)) +
#   geom_line(aes(y = model_fit), col="red") +
#   theme(
#     panel.background = element_blank()
#   )
# gr_plot
r_dat2 <- r_dat %>%
  mutate(growth_rate = (`Mean(R)` - 1) / 98.4,
         year = lubridate::year(date)) %>%
  group_by(year) %>%
  summarise(mean_gr = mean(growth_rate))

ggplot() +
  geom_point(data = r_dat2, aes(x = year, y = mean_gr))

# changepoint analysis of growth rate
library(Rbeast)

out = beast(r_dat2$mean_gr, season="no")
plot(out)
print(out)



