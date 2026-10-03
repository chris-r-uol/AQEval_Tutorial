#!/usr/bin/env Rscript
# Signal isolation in R, for the R half of the website's results.
#
# R's AQEval fits its isolation model with mgcv; the Python port uses pygam.
# The two do not fit identically, and the difference is enough to move where
# the changes are found. So the website carries two sets of results, and shows
# the one for the language selected. This script produces the isolated series
# for the R set; `build_data.py` calls it and does everything else.
#
# The break analysis itself is not repeated in R: the Python port of that part
# reproduces R exactly, and is hundreds of times faster.
#
#   Rscript pipeline/isolate.R pipeline/cache/iso_r
#
# reads `jobs.csv` from that folder and writes one CSV per job beside it.
#
# AQEval needs Java to install (through its 'loa' dependency). Where it is not
# installed, set AQEVAL_SOURCE to the R/ folder of the AQEval source and the
# two files this script needs are read from there instead.

suppressMessages({
  library(openair)
  library(dplyr)
  library(mgcv)
})

if (requireNamespace("AQEval", quietly = TRUE)) {
  suppressMessages(library(AQEval))
  aqeval_version <- as.character(packageVersion("AQEval"))
} else {
  source_dir <- Sys.getenv("AQEVAL_SOURCE")
  if (!nzchar(source_dir)) stop("AQEval is not installed and AQEVAL_SOURCE is not set")
  for (file in c("isolate.signal.R", "aqe.misc.R")) {
    sys.source(file.path(source_dir, file), envir = globalenv())
  }
  aqeval_version <- read.dcf(file.path(source_dir, "..", "DESCRIPTION"))[1, "Version"]
}

out_dir <- commandArgs(trailingOnly = TRUE)[1]
jobs <- read.csv(file.path(out_dir, "jobs.csv"), stringsAsFactors = FALSE, na.strings = "")
# Python writes its booleans as "True" and "False".
jobs$deseason <- toupper(jobs$deseason) == "TRUE"
jobs$deweather <- toupper(jobs$deweather) == "TRUE"
years <- min(jobs$first_year):max(jobs$last_year)

# Download each site once, before any work is shared out.
sites <- unique(c(jobs$site, jobs$control[!is.na(jobs$control)]))
data <- setNames(lapply(sites, function(site) {
  cache <- file.path(out_dir, paste0("raw_", site, ".rds"))
  if (!file.exists(cache)) saveRDS(importAURN(site = tolower(site), year = years), cache)
  readRDS(cache)
}), sites)

isolate <- function(job) {
  out <- file.path(out_dir, paste0(job$key, ".csv"))
  if (file.exists(out)) return(job$key)

  site_data <- data[[job$site]]
  background <- if (is.na(job$background)) NULL else job$background
  if (!is.na(job$control)) {
    # As in the tutorial: the control site's pollutant, joined on by date.
    control <- data[[job$control]][, c("date", job$pollutant)]
    names(control)[2] <- "background"
    site_data <- left_join(site_data, control, by = "date")
    background <- "background"
  }

  # AQEval reports the model it fits as a message.
  fitted <- character()
  isolated <- withCallingHandlers(
    isolateContribution(site_data, job$pollutant,
      background = background, deseason = job$deseason, deweather = job$deweather
    ),
    message = function(m) {
      fitted <<- c(fitted, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    }
  )

  # Written to a temporary name first, so an interrupted run never leaves a
  # half-written file that looks finished.
  partial <- paste0(out, ".part")
  write.csv(
    data.frame(
      date = format(site_data$date, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
      isolated = isolated
    ),
    partial,
    row.names = FALSE
  )
  writeLines(paste(fitted, collapse = " "), file.path(out_dir, paste0(job$key, ".txt")))
  file.rename(partial, out)
  job$key
}

workers <- as.integer(Sys.getenv("WORKERS", "4"))
done <- parallel::mclapply(split(jobs, seq_len(nrow(jobs))), isolate, mc.cores = workers, mc.preschedule = FALSE)

failed <- vapply(done, inherits, logical(1), what = "try-error")
if (any(failed)) stop(sum(failed), " isolation(s) failed: ", done[[which(failed)[1]]])

writeLines(
  sprintf("AQEval %s, mgcv %s, %s", aqeval_version, packageVersion("mgcv"), R.version.string),
  file.path(out_dir, "versions.txt")
)
cat("isolated", length(done), "series in R\n")
