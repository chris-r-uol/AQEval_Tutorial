# Air quality trend detection with AQEval

A tutorial on finding changes in air quality data with [AQEval](https://karlropkins.github.io/AQEval/), using monitoring data from Bradford. It asks one question: did Bradford's Clean Air Zone, which started on 26 September 2022, change the air?

There are three ways in, and they all teach the same analysis.

| | |
| --- | --- |
| **[Interactive tutorial](https://chris-r-uol.github.io/AQEval_Tutorial/)** | Run the experiment with dials instead of code. Choose a site, remove the weather, set the sensitivity, and watch the figure, the findings and the code change together. Nothing to install. |
| **[Run the code](https://codespaces.new/chris-r-uol/AQEval_Tutorial?quickstart=1)** | A ready-made workspace in your browser, with RStudio for the R version and a notebook for the Python version. Every package is already installed. Needs a free GitHub account. |
| **[Step-by-step walkthrough](walkthrough.md)** | The original written practical, in R, with every step explained and its output shown. |

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/chris-r-uol/AQEval_Tutorial?quickstart=1)

## The code

The whole analysis is a dozen lines. In R:

```r
library(openair)
library(AQEval)

data <- importAURN(site = "bdma", year = 2018:2025)

data$isolated <- isolateContribution(
  data, "no2",
  background = "air_temp",
  deseason = TRUE,
  deweather = TRUE
)

data_avg <- timeAverage(data, avg.time = "8 hour")

breaks <- findBreakPoints(data_avg, "isolated", h = 0.3)
result <- quantBreakSegments(data_avg, "isolated", breaks)
result$report
```

And in Python, with [`aqeval`](https://pypi.org/project/aqeval/), a port of the R package:

```python
import aqeval

data = aqeval.download_aurn_data("bdma", 2018, 2025, source="aurn")

data["isolated"] = aqeval.isolate_contribution(
    data, "no2",
    background="air_temp",
    deseason=True,
    deweather=True,
)

data_avg = aqeval.time_average(data, "8 hour")

breaks = aqeval.find_break_points(data_avg, "isolated", h=0.3)
result = aqeval.quant_break_segments(data_avg, "isolated", breaks=breaks)
result["report"]
```

The full versions, with a comment explaining each step, are [`R/aqeval_tutorial.R`](R/aqeval_tutorial.R) and [`python/aqeval_tutorial.ipynb`](python/aqeval_tutorial.ipynb).

## R and Python do not give identical answers

Both languages run the same method, and from the same series they give the same break points and the same changes. What differs is the series: the model that removes the weather and the seasons is fitted with `mgcv` in R and `pygam` in Python. The two fits are close (about 1 µg/m³ apart) but not identical, and that is enough to move the date of a change.

The interactive tutorial therefore holds two complete sets of results, one for each language, and shows the set for the language selected in its code panel. What a student sees on the site is what their own code will print. (The R set is isolated, averaged and measured in R. Its break points are found with the Python port, which returns the same ones as R in a fraction of the time.)

This needs `aqeval` 0.7.3 or later. Earlier versions differed from R in other ways too; [`pipeline/parity`](pipeline/parity/README.md) records what was wrong and holds the checks that now pass.

## What is in this repository

| Path | What it is |
| --- | --- |
| `R/aqeval_tutorial.R` | The tutorial as an R script, for RStudio. |
| `python/aqeval_tutorial.ipynb` | The tutorial as a Python notebook. |
| `python/requirements.txt` | The Python packages, pinned. |
| `walkthrough.md`, `walkthrough.qmd` | The original step-by-step practical and its Quarto source. |
| `site/` | The interactive tutorial: a SvelteKit site published to GitHub Pages. |
| `site/static/data/` | The precomputed results the site shows. |
| `pipeline/` | The scripts that compute those results. |
| `pipeline/parity/` | A check of Python `aqeval` against R AQEval, and what it found. |
| `.devcontainer/` | The definition of the Codespace. |
| `START_HERE.md` | What a student sees when the Codespace opens. |

## Maintaining it

### Changing the code the tutorial teaches

The steps are defined once, in [`site/src/lib/code.js`](site/src/lib/code.js). The website's code panel, the R script and the notebook are all generated from it. After changing a step:

```bash
cd site
npm install
npm run examples
```

The publishing workflow fails if the script or the notebook is out of date.

### Working on the site

```bash
cd site
npm install
npm run dev
```

### Recomputing the results

The site is static, so it cannot run AQEval itself. `pipeline/build_data.py` runs the real analysis for every combination of dials, once in Python and once in R, and writes the answers to `site/static/data`. It caches each stage, so it can be stopped and resumed. Allow three to four hours on a laptop, most of it R measuring the segments.

```bash
python -m venv .venv
.venv/bin/pip install -r pipeline/requirements.txt
.venv/bin/python pipeline/build_data.py
```

The R half needs `openair`, `dplyr`, `jsonlite` and AQEval. AQEval needs Java to install; where that is a problem, download the AQEval source from CRAN and set `AQEVAL_SOURCE` to its `R/` folder. The Codespace has everything.

To add a site, a pollutant or a value of `h`, change the lists at the top of `pipeline/build_data.py` and run it again. The dials read their options from the results.

The Python results are tied to the versions in `python/requirements.txt`, because that is what students run. After changing a version there, delete `pipeline/cache/runs/python` and run the pipeline again.

### Publishing

Pushing to `main` publishes the site through `.github/workflows/pages.yml`. This needs one setting, once: **Settings > Pages > Source: GitHub Actions**.

### Making the Codespace start quickly

Building the workspace from nothing takes several minutes, mostly installing R packages. A prebuild does that ahead of time, so a student's workspace opens in under a minute: **Settings > Codespaces > Set up prebuild**, for the `main` branch. It is worth doing before a class.

## Data and credits

- Measurements: the Automatic Urban and Rural Network, © Crown copyright Defra, via [uk-air.defra.gov.uk](https://uk-air.defra.gov.uk), licensed under the Open Government Licence.
- Clean Air Zone boundary: City of Bradford Metropolitan District Council, Open Government Licence v3.0.
- Method: Ropkins, Walker and Tate, AQEval. [doi:10.21105/joss.08839](https://doi.org/10.21105/joss.08839)
- Base map: © OpenStreetMap contributors.
