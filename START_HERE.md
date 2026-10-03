# Start here

This workspace has everything installed. Pick a language.

## R, in RStudio

1. Click the **Ports** tab in the panel at the bottom of this window.
2. Find the row labelled **RStudio** and click the globe icon beside its address. RStudio opens in a new browser tab.
3. The tutorial script, `R/aqeval_tutorial.R`, opens by itself. If it does not, click the `R` folder in the **Files** pane, then the script. Click on its first line of code, then press **Ctrl+Enter** (**Cmd+Enter** on a Mac) to run one line at a time.

RStudio can take a minute to appear the first time.

Plots appear in the **Plots** pane and tables in the **Console**.

## Python, in a notebook

1. In the file list on the left, open `python/aqeval_tutorial.ipynb`.
2. Click on the first grey code box and press **Shift+Enter** to run it and move to the next.
3. If you are asked to choose a kernel, choose **Jupyter Kernel**, then **Python 3 (AQEval tutorial)**.

Run the boxes in order, from the top.

## How long it takes

Most lines finish in a second or two. A few do not, and it is the step, not you. Roughly:

| What you run | R | Python |
| --- | --- | --- |
| The tutorial as it comes: 8-hour averages, `h = 0.3` | about a minute in all | about three minutes, nearly all of it in "Isolate the signal" |
| A more sensitive search on 8-hour data: `h = 0.2` or less | **tens of minutes**, in "Find the break points" | **several minutes**, in "Measure the changes" |
| A more sensitive search on daily data: `h = 0.12` | two or three minutes | under a minute more |

So in either language, change the averaging to `"day"` before you lower `h`.

## If something goes wrong

- **No RStudio row in the Ports tab.** Open the **Terminal** tab, type `rserver` and press Enter, then look again.
- **A step seems stuck.** In RStudio a red stop sign shows at the top of the Console while R is working. In the notebook the box shows a spinning timer. Check the table above before giving up.
- **You want to start again.** In RStudio choose **Session > Restart R**. In the notebook choose **Restart** at the top, then run the boxes from the top.

## When you have finished

Close the browser tab. The workspace stops by itself after half an hour, and your changes are kept. You can find it again, or delete it, at [github.com/codespaces](https://github.com/codespaces).

## The interactive version

The website shows the same analysis with dials instead of code, and shows this code changing as the dials move: <https://chris-r-uol.github.io/AQEval_Tutorial/>
