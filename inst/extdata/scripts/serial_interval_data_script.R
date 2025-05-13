# In this script we read in and wrangle all of the data sources that we'll use
# to estimate serial interval. We will concetenate all of the data sources into
# a single data frame si_data and save it as an rds file in
# inst/extdata/data/si_data.rds

# load required packages
library(dplyr)
library(tidyr)
library(lubridate)
library(readxl)

# Kaburi et al. ----------------------------------------------------------------
# This paper describes an outbreak of scabies in a preschool in Ghana
# Kaburi et al. 2019 BMC Public Health https://doi.org/10.1186/s12889-019-7085-6

# read in epidemic curve data file
kaburi_df <- read_xlsx("./inst/extdata/data/Kaburi_et_al_data_scabies.xlsx") %>%
  # select the number of days since symptom onset
  select(`number of days since onset`) %>%
  # rename variables for ease
  rename(no_days_since_onset = `number of days since onset`) %>%
  # make an identification variable for whether or not an individual is an index case
  mutate(index_case = if_else(no_days_since_onset == max(no_days_since_onset), 1, 0),
         id = row_number(),
         # calculate ICC intervals
         icc_interval = abs(no_days_since_onset - max(no_days_since_onset)),
         study = "Kaburi et al.",
         date_onset = (as.Date("05-12-2016", format = "%d-%m-%Y") + max(no_days_since_onset) - no_days_since_onset)
  ) %>%
  select(id, date_onset, index_case, icc_interval, study)


# Ariza et al. -----------------------------------------------------------------
# these data are from a scabies outbreak at a pre-school in Germany
ariza_df <- data.frame(
  id = c(seq(1, 16, by = 1)),
  date_onset = as.Date(c("01/12/2011", "25/02/2012", "01/03/2012", "02/03/2012",
                         "05/03/2012", "08/03/2012", "10/03/2012", "10/03/2012",
                         "13/03/2012", "13/03/2012", "14/03/2012", rep("16/03/2012",3),
                         "17/02/2012", "18/03/2012"), format = "%d/%m/%Y"),
  index_case = c(1, rep(0,15))
) %>%
# create icc intervals from date of symptom onset
  mutate(icc_interval = as.integer(date_onset - min(date_onset)),
         study = "Ariza et al.")

# Akunzirwe et al 2023. --------------------------------------------------------
# outbreak in a fishing community in Uganda in 2022.
# Data are provided weekly, so all infections for each week are attributed to
# the first day of the week
# https://uniph.go.ug/an-outbreak-of-scabies-in-a-fishing-community-in-hoima-district-uganda-february%E2%88%92june-2022/

akunzirwe_df <- data.frame(
  date_onset = as.Date(c("01/01/2022", "08/01/2022", "15/01/2022", "29/01/2022",
                         "05/02/2022", "12/02/2022", "19/02/2022", "26/02/2022",
                         "05/03/2022", "12/03/2022", "19/03/2022", "26/03/2022",
                         "02/04/2022", "09/04/2022", "16/04/2022", "23/04/2022",
                         "30/04/2022", "07/05/2022", "14/05/2022", "21/05/2022",
                         "28/05/2022", "11/06/2022", "18/06/2022", "25/06/2022",
                         "02/07/2022"), format = "%d/%m/%Y"),
  num_cases = c(3, 4, 2, 5, 3, 2, 3, 6, 6, 4, 3, 2, 6, 20, 7, 4, 14, 8, 18, 2,
                10, 5, 10, 5, 3)
) %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         index_case = if_else(date_onset == min(date_onset), 1, 0),
         icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
         study = "Akunzirwe et al."
  ) %>%
  select(id, date_onset, index_case, icc_interval, study)

# Tjon-Kon-Fat et al. ----------------------------------------------------------
# outbreak in a care home in the Netherlands
# https://doi.org/10.1371/journal.pntd.0009485
# cases identified by week, date of onset is attributed to first day of the week

tkf_df <- data.frame(
  week_onset = c(24, 28, 31, 33, 36, 37, 38, 40, 41, 42),
  num_cases = c(1, 1, 1, 2, 1, 2, 1, 1, 9, 2)
)

# Define the year
year <- 2018

# Convert week number to the first day of each week using piping
tkf_df <- tkf_df %>%
  mutate(date_onset = ymd(paste0(year, "-01-01")) + weeks(week_onset - 1)) %>%
  mutate(date_onset = floor_date(date_onset, unit = "week", week_start = 1)) %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         index_case = if_else(date_onset == min(date_onset), 1, 0),
         icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
         study = "Tjon-Kon-Fat et al."
  ) %>%
  select(id, date_onset, index_case, icc_interval, study)

# Wochebo et al. 2019
# outbreak in Southern Ethiopia
# https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-019-4317-x

# wochebo_df <- data.frame(
#   date_onset = as.Date(c("2017-05-19", "2017-06-06", "2017-06-08", "2017-11-10",
#                    "2017-06-11", "2017-06-12", "2017-06-13", "2017-06-14",
#                    "2017-06-15", "2017-06-16", "2017-06-17", "2017-06-18",
#                    "2017-06-19", "2017-06-20", "2017-06-21", "2017-06-22",
#                    "2017-06-23")),
#   cases = c(1, 6, 4, 28, 33, 29, 33, 33, 23, 19, 10, 13, 11, 12, 9, 7, 8)
# ) %>%
#   uncount(cases) %>% # Repeat rows based on the number of cases
#   mutate(id = row_number(), # Add an id variable
#          index_case = if_else(date_onset == min(date_onset), 1, 0),
#          icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
#          study = "Wochebo et al."
#   )
#
# # Tsutsumi et al. 2005
# # outbreak in patients with dementia in a long-term care facility
# # https://doi.org/10.1186/1471-2334-5-85
#
# tsutsumi_df <- data.frame(
#   date_onset = as.Date(c("1989-05-06", "1989-06-03", "1989-06-17", "1989-07-15", "1989-08-26",
#            "1989-09-21", "1989-09-28", "1989-10-05", "1989-10-19", "1989-10-26",
#            "1989-11-02", "1989-11-09", "1989-11-30", "1989-12-07")),
#   cases = c(1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 1, 2, 1, 2)
# ) %>%
#   uncount(cases) %>% # Repeat rows based on the number of cases
#   mutate(id = row_number(), # Add an id variable
#          index_case = if_else(date_onset == min(date_onset), 1, 0),
#          icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
#          study = "Tsutsumi et al."
#   )


# Combine all studies together and output as RDS file
si_data <- bind_rows(kaburi_df, ariza_df, akunzirwe_df, tkf_df
                     #, wochebo_df,tsutsumi_df
                     )
saveRDS(si_data, "inst/extdata/data/si_data.rds")

