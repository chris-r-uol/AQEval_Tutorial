# Open the tutorial script when RStudio starts in this folder.
#
# Only in RStudio, and only for a new session, so it does not reopen the file
# every time R restarts or get in the way of running R from a terminal.
setHook("rstudio.sessionInit", function(newSession) {
  script <- "R/aqeval_tutorial.R"
  if (newSession && file.exists(script) && requireNamespace("rstudioapi", quietly = TRUE)) {
    rstudioapi::navigateToFile(script)
  }
}, action = "append")
