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

# Division of Public and Behavioral Health (DPBH) ------------------------------
# Scabies outbreak among residents of a long term care facility in Clark County,
# Nevada, USA in 2015
dpbh_df <- data.frame(
  date_onset = as.Date(c("07/04/2015", "08/04/2015", "16/04/2015", "17/04/2015",
                         "22/04/2015", "01/05/2015", "14/05/2015", "16/05/2015",
                         "18/05/2015"), format = "%d/%m/%Y"),
  num_cases = c(1, 1, 2, 1, 1, 1, 1, 1, 2)) %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         index_case = if_else(date_onset == min(date_onset), 1, 0),
         icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
         study = "Division of Public and Behavioural Health"
         ) %>%
  select(id, date_onset, index_case, icc_interval, study)

# Larrosa et al. ---------------------------------------------------------------
# Outbreak in a hospital in Spain in 2002
# https://doi.org/10.2807/esm.08.10.00429-en

# Define start date
start_date <- as.Date("2002-11-05")

larrosa_df <- data.frame(
  day_onset = c(1, 2, 5, 9, 15, 16, 19, 21, 29, 44, 51, 61, 62),
  num_cases = c(1, 1, 2, 1, 1, 3, 1, 1, 1, 2, 1, 1, 1)
) %>%
  mutate(date_onset = start_date + (day_onset - 1)) %>%
  uncount(num_cases) %>% # Repeat rows based on the number of cases
  mutate(id = row_number(), # Add an id variable
         index_case = if_else(date_onset == min(date_onset), 1, 0),
         icc_interval = as.integer(date_onset - min(date_onset)), # calculate ICC interval
         study = "Larrosa et al."
  ) %>%
  select(id, date_onset, index_case, icc_interval, study)


# Combine all studies together and output as RDS file
si_data <- bind_rows(kaburi_df, ariza_df, akunzirwe_df, tkf_df, dpbh_df, larrosa_df)
saveRDS(si_data, "inst/extdata/data/si_data.rds")

