# Calculate reproduction number from scabies time series

# 1. We have reported cases of scabies reported weekly. We will first randomly
#    assign each date of symptom onset to a day of the reported week. This will
#    give us a daily time series instead of weekly.
#
# 2. We will use the method proposed by Wallinga and Teunis 2004 to estimate the
#    time-varying reproduction number
#    a. we will first assume a Normal serial interval distribution
#    b. next, we will assume a Gamma serial interval distribution

# load required packages -------------------------------------------------------
library(dplyr)
library(tidyr)
library(lubridate)
library(ISOweek)
library(readxl)
library(devtools)
load_all()

# read in data -----------------------------------------------------------------
nivel_wkly_data <-
  read_xlsx("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies/data/scabies_data_weekly.xlsx") %>%
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

# 1. randomly assign a date within reporting date for symptom onset date -------
nivel_daily_data <- nivel_wkly_data %>%
  uncount(cases) %>% # Repeat rows based on the number of cases
  mutate(
    iso_week = paste0(year, "-W", sprintf("%02d", as.numeric(week_num))),
    first_day = ISOweek2date(paste0(iso_week, "-1")),
    random_day = sample(0:6, n(), replace = TRUE),
    onset_date = first_day + days(random_day)
  ) %>%
  select(-iso_week, -first_day, -random_day)

# 2. We will implement the method from Wallinga and Teunis 2004 ----------------
# To begin, we will try a smaller example and only use data from the first week.
week1 <- nivel_daily_data %>%
  filter(year == 2011, week_num == 1) %>%
  # order by onset date
  arrange(onset_date)

# Determine likelihood of transmission pair for every possible pair
# we assume the serial interval distribution is N(95.57, (15.17)^2)
# we will loop over each observation and calculate the likelihood of the
# transmission pair.
for (i in 1:dim(week1)[1]){

}
