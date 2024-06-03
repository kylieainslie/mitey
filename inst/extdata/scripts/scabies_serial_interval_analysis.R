# Determining serial interval of scabies using epidemic curves from different
# sources
# We will start with data from an outbreak of scabies in a preschool from
# Kaburi et al. 2019 BMC Public Health https://doi.org/10.1186/s12889-019-7085-6

# load required packages
library(tidyverse)
library(readxl)

# set working directory
setwd("~/Dropbox/Kylie/Projects/RIVM/Projects/scabies")

# Kaburi et al. ----------------------------------------------------------------
# read in epidemic curve data file
kaburi_df <- read_xlsx("./data/Kaburi_et_al_data_scabies.xlsx")

# we will calculate the index case-to-case (ICC) interval for each person by class
# the person with the greatest value for number of days since symptom onset will
# be considered the index case. The rest of the class members will have an ICC
# interval calculated as the number of days between their symptom onset and the
# index case.

icc_df <- kaburi_df %>%
  # for now, we only need the class ID and the number of days since symptom onset
  select(`Class (0=Creche; 1=Nursery 1; 2=Nursery 2; 3=KG1; 4=KG2)`,
         `number of days since onset`) %>%
  # rename variables for ease
  rename(class = `Class (0=Creche; 1=Nursery 1; 2=Nursery 2; 3=KG1; 4=KG2)`,
         no_days_since_onset = `number of days since onset`) %>%
  # make an identification variable for whether or not an individual is an index case
  group_by(class) %>%
  mutate(index_case = if_else(no_days_since_onset == max(no_days_since_onset), 1, 0),
         # calculate ICC intervals
         icc_interval = abs(no_days_since_onset - max(no_days_since_onset))
         )

# use method from Vink et al. to estimate SI
data <- icc_df$icc_interval
source("./code/SI_estimation_method_from_Vink_et_al.R")
# [1] mean = 81.36419, sd = 62.89626

# plot
hist(data, freq = FALSE, main = "Scabies ICC intervals", xlab = "ICC Interval (days)")
x <- seq(min(icc_df$icc_interval), max(icc_df$icc_interval), length.out = 100)
y <- dnorm(x, mean = 81.36419, sd = 62.89626)
lines(x, y, col = "blue", lwd = 2) 
# ------------------------------------------------------------------------------

# Ariza et al. -----------------------------------------------------------------
ariza_df <- data.frame(
  id = c(seq(1, 16, by = 1)),
  date_onset = as.Date(c("01/12/2011", "25/02/2012", "01/03/2012", "02/03/2012", 
                 "05/03/2012", "08/03/2012", "10/03/2012", "10/03/2012",
                 "13/03/2012", "13/03/2012", "14/03/2012", rep("16/03/2012",3),
                 "17/02/2012", "18/03/2012"), format = "%d/%m/%Y"),
  index_case = c(1, rep(0,15))
)

# calculate icc intervals from date of symptom onset
icc_df2 <- ariza_df %>%
  mutate(icc_interval = as.integer(date_onset - min(date_onset)))

# use method from Vink et al. to estimate SI
data <- icc_df2$icc_interval
source("./code/SI_estimation_method_from_Vink_et_al.R")
# [1] mean = 98.4, sd = 8.542332

# plot
hist(data, freq = FALSE, ylim = c(0, 0.05), xlim = c(0, 130),
     main = "Scabies ICC intervals", xlab = "ICC Interval (days)")
x <- seq(min(icc_df2$icc_interval), max(icc_df2$icc_interval), length.out = 100)
y <- dnorm(x, mean = 98.4, sd = 8.542332)
lines(x, y, col = "blue", lwd = 2) 

# ------------------------------------------------------------------------------
# plot both estimated serial interval distributions
x <- seq(0,150, by = 1)
y1 <- dnorm(x, mean = 81.36419, sd = 62.89626)
y2 <- dnorm(x, mean = 98.4, sd = 8.542332)
plot(y2~x, type = "l", lwd = 2, ylab = "Density", xlab = "Time (days)",
     main = "Estimated Serial Interval Distributions")
lines(x, y1, col = "blue", lwd = 2)
legend(1, 0.045, legend=c("Ariza et al.", "Kaburi et al."),
       col=c("black", "blue"), lty=c(1,1), cex=1)
