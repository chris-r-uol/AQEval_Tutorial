#!/usr/bin/env Rscript
# The break analysis in R, for the R half of the website's results.
#
# The Python port does not always fit the same segments as R, and averages
# weeks from a different day (see pipeline/parity/README.md), so the R results
# on the website are averaged and measured by R itself: openair's timeAverage
# and AQEval's quantBreakSegments, as the tutorial script runs them.
#
# One step is borrowed. R's search for break points, findBreakPoints, takes a
# quarter of an hour on an 8-hour series at the more sensitive settings, and
# there are thousands of searches to do. The Python port of that step returns
# the same break points and confidence intervals as R in every case checked,
# in under a second. So this script runs in two stages, with Python finding
# the breaks in between:
#
#   STAGE=average   average every series and write the result out
#   (build_data.py finds the break points in those averaged series)
#   STAGE=quantify  measure the segments for each set of breaks
#
#   STAGE=average Rscript pipeline/breaks.R pipeline/cache/iso_r pipeline/cache/breaks_r
#
# The first folder holds the site data and the isolated series written by
# isolate.R. The second holds `jobs.csv` and receives everything else.

suppressMessages({
  library(openair)
  library(jsonlite)
})

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

args <- commandArgs(trailingOnly = TRUE)
in_dir <- args[1]
out_dir <- args[2]
stage <- Sys.getenv("STAGE")

jobs <- read.csv(file.path(out_dir, "jobs.csv"), stringsAsFactors = FALSE, na.strings = "")
averaging <- strsplit(Sys.getenv("AVERAGING"), ";")[[1]]

quietly <- function(expr) suppressWarnings(suppressMessages(expr))

# Written under a temporary name first, so an interrupted run never leaves a
# half-written file that looks finished.
finish <- function(partial, path) file.rename(partial, path)

# --- Stage 1: average ----------------------------------------------------------

average <- function(job) {
  out <- file.path(out_dir, paste0(job$key, ".averaged.csv"))
  if (file.exists(out)) return(job$key)

  if (job$isolating) {
    series <- read.csv(file.path(in_dir, paste0(job$key, ".csv")))
    series$date <- as.POSIXct(series$date, tz = "UTC")
  } else {
    series <- readRDS(file.path(in_dir, paste0("raw_", job$site, ".rds")))[, c("date", job$pollutant)]
  }
  names(series)[2] <- "value"

  rows <- do.call(rbind, lapply(averaging, function(period) {
    averaged <- as.data.frame(timeAverage(series, avg.time = period))
    data.frame(
      period = period,
      epoch = sprintf("%.0f", as.numeric(averaged$date)),
      # Every digit of the double, so that Python reads back exactly this number.
      value = ifelse(is.na(averaged$value), "NA", sprintf("%.17g", averaged$value))
    )
  }))
  partial <- paste0(out, ".part")
  write.csv(rows, partial, row.names = FALSE, quote = FALSE)
  finish(partial, out)
  job$key
}

# --- Stage 2: quantify ---------------------------------------------------------

quantify <- function(job) {
  out <- file.path(out_dir, paste0(job$key, ".json"))
  if (file.exists(out)) return(job$key)

  averaged_all <- read.csv(file.path(out_dir, paste0(job$key, ".averaged.csv")), stringsAsFactors = FALSE)
  # One row per break; a search that found none has a single row with no numbers.
  breaks_all <- read.csv(file.path(out_dir, paste0(job$key, ".breaks.csv")), stringsAsFactors = FALSE)

  result <- list()
  for (period in averaging) {
    one <- averaged_all[averaged_all$period == period, ]
    averaged <- data.frame(date = as.POSIXct(one$epoch, origin = "1970-01-01", tz = "UTC"), value = one$value)
    searches <- breaks_all[breaks_all$period == period, ]

    fits <- list()        # one per distinct set of breaks
    signatures <- character()
    runs <- list()

    for (h in unique(searches$h)) {
      found <- searches[searches$h == h, ]
      status <- found$status[1]
      # Handed NULL, quantBreakSegments fits a single straight line.
      breaks <- if (status == "no_breaks") NULL else found[, c("lower", "bpt", "upper")]
      if (!is.null(breaks)) rownames(breaks) <- NULL
      run <- list(breaks = if (is.null(breaks)) list() else unname(as.matrix(breaks)), fit = NULL, status = status)

      if (status != "too_many_breaks") {
        signature <- paste(unlist(breaks), collapse = ",")
        index <- match(signature, signatures)
        if (is.na(index)) {
          fitted <- tryCatch(
            quietly(quantBreakSegments(averaged, "value", breaks, show = "none")),
            error = function(e) e
          )
          signatures <- c(signatures, signature)
          index <- length(signatures)
          if (inherits(fitted, "error")) {
            fits[[index]] <- list(error = substr(conditionMessage(fitted), 1, 200))
          } else {
            report <- fitted$report
            fits[[index]] <- list(
              pred = fitted$data2$pred,
              err = fitted$data2$err,
              report = if (is.null(report)) NULL else list(
                date1 = as.numeric(report$s1.date1),
                date2 = as.numeric(report$s1.date2),
                days = as.numeric(report$s1.date.delta, units = "days"),
                c0 = report$s1.c0,
                c1 = report$s1.c1,
                change = report$s1.c.delta,
                percent = report$s1.per.delta
              )
            )
          }
        }
        run$fit <- index
      }
      runs[[as.character(h)]] <- run
    }

    result[[period]] <- list(epoch = one$epoch, values = one$value, fits = fits, runs = runs)
  }

  partial <- paste0(out, ".part")
  write_json(list(key = job$key, averaging = result), partial, auto_unbox = TRUE, na = "null", digits = NA, null = "null")
  finish(partial, out)
  job$key
}

# --- Run -------------------------------------------------------------------------

task <- switch(stage, average = average, quantify = quantify, stop("STAGE must be 'average' or 'quantify'"))
workers <- as.integer(Sys.getenv("WORKERS", "4"))
done <- parallel::mclapply(split(jobs, seq_len(nrow(jobs))), task, mc.cores = workers, mc.preschedule = FALSE)

failed <- vapply(done, inherits, logical(1), what = "try-error")
if (any(failed)) stop(sum(failed), " job(s) failed: ", done[[which(failed)[1]]])
cat(stage, "done for", length(done), "series\n")
