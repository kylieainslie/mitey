library(httr2)
library(jsonlite)
library(dplyr)
library(tibble)
library(lubridate)

lat <- 48.8566
lon <- 2.3522

start_date <- "2025-01-01"
end_date   <- "2025-12-31"

hourly_vars <- c(
  "temperature_2m"
)

resp <- request("https://archive-api.open-meteo.com/v1/archive") |>
  req_url_query(
    latitude = lat,
    longitude = lon,
    start_date = start_date,
    end_date = end_date,
    hourly = paste(hourly_vars, collapse = ","),
    timezone = "Europe/Paris"
  ) |>
  req_perform()

data_raw <- fromJSON(resp_body_string(resp), flatten = TRUE)

# Vérifie ce qu'il y a dans hourly
names(data_raw$hourly)

meteo_paris_2025 <- as_tibble(data_raw$hourly) |>
  mutate(
    time = ymd_hm(time, tz = "Europe/Paris")
  )

glimpse(meteo_paris_2025)
head(meteo_paris_2025)

write.csv(meteo_paris_2025, "meteo_paris_2025.csv", row.names = FALSE)
