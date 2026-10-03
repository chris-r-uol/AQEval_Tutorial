#!/usr/bin/env Rscript
# R's half of the parity check: what does R AQEval return for each case?
#
# For every row of cases.csv this averages the series with openair, runs
# findBreakPoints and quantBreakSegments, and writes out the exact input it
# analysed and everything it got back. compare.py then gives Python the same
# input and compares the two, value for value.
#
#   Rscript pipeline/parity/r_side.R
#
# Needs the series cached by the main pipeline (pipeline/cache/iso_r). Cases
# with more than three breaks are searched but not measured, to keep this to a
# few minutes.

suppressMessages(library(openair))
if (requireNamespace("AQEval", quietly = TRUE)) {
  suppressMessages(library(AQEval))
} else {
  # See isolate.R: AQEval's source, where the package cannot be installed.
  source_dir <- Sys.getenv("AQEVAL_SOURCE")
  if (!nzchar(source_dir)) stop("AQEval is not installed and AQEVAL_SOURCE is not set")
  suppressMessages({
    library(strucchange); library(mgcv); library(segmented); library(dplyr); library(ggplot2)
  })
  for (file in list.files(source_dir, full.names = TRUE)) {
    tryCatch(sys.source(file, envir = globalenv()), error = function(e) NULL)
  }
}

cache <- "pipeline/cache/iso_r"
out <- "pipeline/cache/parity"
dir.create(out, showWarnings = FALSE, recursive = TRUE)
cases <- read.csv("pipeline/parity/cases.csv", stringsAsFactors = FALSE)

# Every digit of a double, so that Python reads back exactly the same number.
exact <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
save <- function(frame, case, what) {
  write.csv(frame, file.path(out, paste0(case, "_", what, ".csv")), row.names = FALSE, quote = FALSE)
}

for (i in seq_len(nrow(cases))) {
  case <- cases[i, ]
  if (case$source == "raw") {
    series <- readRDS(file.path(cache, paste0("raw_", case$site, ".rds")))[, c("date", case$pollutant)]
  } else {
    series <- read.csv(file.path(cache, paste0(case$key, ".csv")))
    series$date <- as.POSIXct(series$date, tz = "UTC")
  }
  names(series)[2] <- "value"
  averaged <- as.data.frame(timeAverage(series, avg.time = case$averaging))
  save(data.frame(epoch = sprintf("%.0f", as.numeric(averaged$date)), value = exact(averaged$value)), case$id, "input")

  breaks <- findBreakPoints(averaged, "value", h = case$h)
  found <- if (is.null(breaks)) 0 else nrow(breaks)
  save(if (found) breaks else data.frame(lower = integer(), bpt = integer(), upper = integer()), case$id, "breaks")

  if (found <= 3) {
    fitted <- suppressWarnings(suppressMessages(quantBreakSegments(averaged, "value", breaks, show = "none")))
    save(data.frame(pred = exact(fitted$data2$pred)), case$id, "trend")
    if (!is.null(fitted$segments)) save(fitted$segments, case$id, "segments")
  }
  cat(sprintf("%s  %s %s %s %s h=%s  breaks=%d\n", case$id, case$source, case$site, case$pollutant,
              case$averaging, case$h, found))
}
