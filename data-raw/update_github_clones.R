# Daily GitHub clones for mitey, accumulated over time.
# GitHub's traffic API only keeps the last 14 days, so this must run at least
# every 2 weeks. New data is merged into the existing CSV; where dates overlap
# the newest values win (the most recent day is partial until the next run).
# Needs the `gh` CLI, authenticated with a token that has admin read access
# to the repo (in GitHub Actions: GH_TOKEN set from the TRAFFIC_TOKEN secret).
# Run from the package root.

library(jsonlite)

repo <- "kylieainslie/mitey"
file <- "data-raw/github_clones_daily.csv"

res <- system2("gh", c("api", paste0("repos/", repo, "/traffic/clones")),
               stdout = TRUE)
if (!is.null(attr(res, "status"))) stop("gh api call failed")

clones <- fromJSON(paste(res, collapse = "\n"))$clones
if (length(clones) == 0 || nrow(clones) == 0) {
  message("No clone data in the last 14 days; nothing to update")
  quit(save = "no")
}

new <- data.frame(
  date    = as.Date(substr(clones$timestamp, 1, 10)),
  clones  = clones$count,
  uniques = clones$uniques
)

if (file.exists(file)) {
  old <- read.csv(file, colClasses = c(date = "Date"))
  old <- old[!old$date %in% new$date, ]
  new <- rbind(old, new)
}
new <- new[order(new$date), ]

write.csv(new, file, row.names = FALSE)
message("Wrote ", nrow(new), " days (", sum(new$clones), " clones)")
