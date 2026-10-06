# Finding changes in Bradford's air quality with AQEval
#
# This is the R version of the interactive tutorial at
# https://chris-r-uol.github.io/AQEval_Tutorial/
#
# How to use it
# Click on the first line of code, then press Ctrl+Enter (Cmd+Enter on a Mac)
# to run it and move to the next. Plots appear in the Plots pane and tables in
# the Console.
#
# This file is generated from site/src/lib/code.js, so that it always matches
# the website. Change it as much as you like while you work.

# 1. Load the packages -----------------------------------------------------
# openair downloads UK air quality data and AQEval finds the changes in it.
library(openair)
library(AQEval)

# 2. Get the data ----------------------------------------------------------
# Download every hourly measurement from 2018 to 2025 for one monitoring site
# on the national network (AURN). Each site has a short code. The download
# also brings the wind speed, wind direction and air temperature for each
# hour.
data <- importAURN(site = "bdma", year = 2018:2025)

# 3. Look at the data ------------------------------------------------------
# Always plot the measurements before analysing them. Look for gaps, and for
# anything that seems out of place.
timePlot(data, pollutant = "no2")

# 4. Isolate the signal ----------------------------------------------------
# Pollution rises and falls with the weather, the time of day and the time of
# year. This step fits a model of those patterns and takes them away, leaving
# the part of the signal that a change in emissions could explain.
data$isolated <- isolateContribution(
  data, "no2",
  background = "air_temp",
  deseason = TRUE,
  deweather = TRUE
)

# 5. Average the data ------------------------------------------------------
# Hourly values are noisy. Averaging them into longer blocks smooths out the
# short-lived ups and downs, and gives the next step fewer points to work
# through.
data_avg <- timeAverage(data, avg.time = "8 hour")

# 6. Find the break points -------------------------------------------------
# A break point is a moment where the average level of the series shifts. h
# sets the shortest stretch allowed between two breaks, as a fraction of the
# whole series: a smaller h can find more breaks, closer together, and takes
# longer. 0.3 is a quick first look. For 8-hour data we would normally use
# 0.12.
breaks <- findBreakPoints(data_avg, "isolated", h = 0.3)
breaks

# 7. Measure the changes ---------------------------------------------------
# Real changes rarely happen in a single day. This step fits a line through
# each stretch of the series and reports when each change started, when it
# finished, and how big it was.
result <- quantBreakSegments(data_avg, "isolated", breaks)
result$report

# Things to try ------------------------------------------------------------
# 1. Change "no2" to "nox". Do the changes fall on the same dates?
# 2. Change the site from "bdma" to "led6" (Leeds Headingley Kerbside), a
#    roadside site with no Clean Air Zone.
# 3. Make the analysis more sensitive: change h from 0.3 to 0.2, then to
#    0.12, the value normally used for 8-hour data. Smaller values find more
#    changes. In R the search then takes tens of minutes on 8-hour data, so
#    change "8 hour" to "day" first: on daily data it takes a minute or two.
# 4. Use a background site as a control: the website shows the extra lines
#    when you choose one under "Also account for".
