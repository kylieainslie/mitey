# Monthly CRAN downloads for mitey (Posit / cloud.r-project.org mirror only).
# Rebuilds the full history each run, so it is safe to re-run at any time.
# Only complete months are written. Run from the package root.

library(cranlogs)

pkg   <- "mitey"
start <- as.Date("2025-09-01")                           # first CRAN release month
end   <- as.Date(format(Sys.Date(), "%Y-%m-01")) - 1     # end of last complete month

daily <- cran_downloads(packages = pkg, from = start, to = end)
if (nrow(daily) == 0) stop("No data returned from cranlogs")

daily$month <- format(daily$date, "%Y-%m")
monthly <- aggregate(count ~ package + month, data = daily, FUN = sum)
monthly <- monthly[order(monthly$month), ]

write.csv(monthly, "data-raw/cran_downloads_monthly.csv", row.names = FALSE)
message("Wrote ", nrow(monthly), " months (", sum(monthly$count), " downloads)")
