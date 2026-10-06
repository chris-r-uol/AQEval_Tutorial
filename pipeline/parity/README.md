# Why Python `aqeval` and R AQEval give different results

Checked on 2 October 2026: Python `aqeval` 0.7.2 (PyPI, identical to GitHub `main` at `248545a`) with pandas 3.0.6 and numpy 2.5.3, against R AQEval 0.6.11 with openair 3.1.0 and mgcv 1.9-4. Data: AURN hourly NO₂ and NOx for BDMA, LED6, LEED and DYAG, 2018 to 2025.

> **Fixed in `aqeval` 0.7.3**, released the same day. Against it, `compare.py` gives 20 of 20 for break points and segments in all three modes. `compare_grid.py` gives Python every series R measured for the website, 1,332 fits in all: Python chooses the same fit as R in every one, with the fitted trends agreeing to within 1e-6 µg/m³ in all but three, which differ by about that much. For cause 2 the fix ports R's own sequence of fits rather than applying `segmented.patch`, which also settles the question about tolerances left open below. Cause 5, the isolation model, is not a defect and remains. What follows describes 0.7.2 and is kept as the record of what was wrong.

## Summary

Given the same hourly data, the two packages agree on the break points and disagree on almost everything after them. There are five separate causes. Four are defects in the port; the fifth is the known difference between the GAM libraries.

| # | Cause | Affects | Size of effect |
| --- | --- | --- | --- |
| 1 | Dates are assumed to be nanosecond-resolution. Under pandas 3 they are not. | `test_break_points`, `quant_break_points`, `quant_break_segments`, `spectral_frequency` | Large. The break-point test rejects every break. The segment fit picks a different model in most cases. |
| 2 | A trial fit that is rank-deficient is scored as a failure, where R scores it by its residual sum of squares. A rank-deficient final fit is not detected at all. | `quant_break_segments` | Occasional, but large when it happens: a different model wins. |
| 3 | Multi-day averaging periods start on a different day from openair 3.1. | `time_average(x, "7 day")` and everything downstream of it | Every weekly value differs. |
| 4 | `breaks=None` is treated as "not given". R treats `NULL` as "no breaks". | `quant_break_segments`, `quant_break_points` | Only when no breaks are found. |
| 5 | pygam is not mgcv. | `isolate_contribution` | Small and expected: r ≈ 0.998, about 1 µg/m³. |

What does agree exactly: the downloaded data (every value, all four sites), `time_average` for 8-hour and daily periods, and `find_break_points` (20 of 20 cases, including the confidence intervals).

Fixing 1 and 2 was tested here, without changing the package, and brings `quant_break_segments` from 8 of 20 cases matching R to 20 of 20. The change for 2 is in [`segmented.patch`](segmented.patch).

## The evidence

`r_side.R` runs R on 20 cases and writes out the exact averaged series it analysed and everything it returned. `compare.py` gives Python that same series and compares break table, segment table and fitted trend.

| Python given R's averaged series | Break points | Segments |
| --- | --- | --- |
| as shipped, with the second-resolution dates `download_aurn_data` and `time_average` produce | 20 of 20 | **8 of 20** |
| dates converted to nanoseconds first (`--ns`) | 20 of 20 | **18 of 20** |
| nanosecond dates and `segmented.patch` | 20 of 20 | **20 of 20** |

Where the segments differ the difference is not small. The fitted trend is up to 13 µg/m³ away from R's and the changes are reported on different dates. For the tutorial's own default set-up over 2018 to 2025 (BDMA, NO₂, isolated, 8-hour, h = 0.3) R reports a 7% rise in August 2020 and Python an 11% fall in March 2020.

## Cause 1: dates are assumed to be nanoseconds

pandas 3 no longer stores every datetime in nanoseconds. It keeps the resolution it was given:

| How the dates were made | pandas 2 | pandas 3 |
| --- | --- | --- |
| `pd.to_datetime(x, unit="s")`, which is what `download_aurn_data` does | `datetime64[ns]` | `datetime64[s]` |
| parsed from text, which is what `load_aq_data` does | `datetime64[ns]` | `datetime64[us]` |
| `pd.date_range(...)` | `datetime64[ns]` | `datetime64[us]` |

`time_average` passes the resolution through. So the package's own loaders now hand its analysis functions dates that four lines misread:

```
aqeval/find_breaks.py:51          epoch = data["date"].astype("int64").to_numpy() // 10**9
aqeval/quantify_breaks.py:282     epoch = data2["date"].astype("int64").to_numpy() / 10**9
aqeval/spectral_analysis.py:195   ts_sec = ts.asi8 / 10**9
aqeval/spectral_analysis.py:202   x_sec = df["date"].astype("int64").to_numpy()[ok] / 10**9
```

Each divides by 10⁹ to turn nanoseconds into seconds. When the integers are already seconds, the result is not seconds.

### `find_breaks.py:51`: the break-point model loses its time variable

This is the serious one. It is an integer division, so for second-resolution dates in this century `epoch` is `1` on every row:

```
dates stored as datetime64[s]    epoch = [1, 1, ..., 1]
dates stored as datetime64[ns]   epoch = [1514764800, 1514851200, ..., 1767139200]
```

`_fit_break_points_model` uses `epoch` as the regressor within each segment. With a constant there, each segment's "slope" column is a copy of its indicator, the columns sum to the intercept, and the model is a set of flat steps fitted to a singular design. On BDMA daily NO₂ with h = 0.3:

| | `datetime64[s]` | `datetime64[ns]` |
| --- | --- | --- |
| `test_break_points` | every model with a break is "not significant", adjusted R² is NaN, and it suggests **no breaks** | all models significant, suggests break 2 |
| `quant_break_points` | 41.2 → 35.5 (−14%) and 35.5 → 28.4 (−20%) | 38.2 → 38.8 (+1.5%) and 32.7 → 31.1 (−5%) |
| `quant_break_segments(data, "value", h=0.3)`, the default route that finds, tests, then measures | a single straight line, no report | 3 segments |

So anyone calling `quant_break_segments` or `quant_break_points` without supplying `breaks`, on data from `download_aurn_data`, under pandas 3, is told there are no changes. With microsecond dates (the bundled `load_aq_data`) the division gives whole kiloseconds: the model still runs but time is rounded down to the nearest 1,000 seconds.

`dt_app` calls exactly this route (`aqeval.quant_break_segments(isolated, "isolated", h=h, show=())`). Its `requirements.txt` pins pandas 2.2.3, so the deployed service is on nanoseconds and unaffected, but its local `.venv` has pandas 3.0.5.

### `quantify_breaks.py:282`: a last-digit error that changes the answer

Here the division is a float division and the result is then rescaled to run from 0 to 1, so the wrong unit cancels out mathematically. It does not cancel exactly. Dividing by 10⁹ first introduces rounding, and the rescaled time `d_prop` ends up a few units in the last place away from R's `(t - min(t)) / max(t - min(t))`:

```
row 131, dates as [s]:   0.31175059952038375
row 131, dates as [us]:  0.3117505995203836
```

That would not matter for most algorithms. It matters here because of how AQEval starts the segmented fit. The starting breakpoints are built from `d_prop` at the break rows, offset by whole multiples of the distance to the confidence limits, so they sit exactly on, or one rounding error away from, data points. The fit begins by asking of every observation whether `z > psi`. For the observation that psi is sitting on, the answer is decided by the last digit. A different answer puts that point on the other side of the breakpoint, the first Muggeo step moves somewhere else, and after 9^k candidates a different one has the best adjusted R².

With `d_prop` computed in exact seconds it is bit-for-bit what R computes, the comparisons come out the same, and the fits agree to about 1e-9. As shipped, the answer depends on how the dates happen to be stored: of the 20 cases, the segments match R in 8 with the dates as `[s]`, 6 as `[us]` and 18 as `[ns]`.

### `spectral_analysis.py:195` and `:202`

Not tested here. The two lines scale two sets of dates that can have different resolutions (a `date_range` and the data's own column). If they differ, `np.interp` is interpolating between axes in different units.

### The fix

Compute seconds without assuming a resolution, in one helper used by all four:

```python
def _epoch_seconds(dates) -> np.ndarray:
    """Seconds since 1970 as floats, whatever resolution the dates are stored in."""
    return pd.to_datetime(dates).astype("datetime64[ns]").astype("int64").to_numpy() / 10**9
```

and use `np.floor(...)` or `// 1` where whole seconds are wanted. Whole seconds are exactly representable, so this reproduces R's arithmetic. A test that runs the same series as `[s]`, `[us]` and `[ns]` and requires identical output would have caught all of this.

## Cause 2: rank-deficient fits in the segmented step

`aqeval/_segmented.py`, `fit_segmented_lm`. Two differences from R's `local_seg.lm.fit` and `local_segmented`, both about what happens when a design matrix loses rank.

### 2a. A trial fit that cannot be fully estimated

After the Muggeo update `psi = psi_old + gamma / beta`, R fits `y ~ [1, z, U]` at the new psi to see whether the step improved things:

```r
obj1 <- try(mylm(cbind(XREG, U1), y, w, offs), silent = TRUE)
if (class(obj1)[1] == "try-error")
  obj1 <- try(lm.wfit(cbind(XREG, U1), y, w, offs), silent = TRUE)
L1 <- if (class(obj1)[1] == "try-error") L0 + 10 else sum(obj1$residuals^2 * w)
```

`mylm` fails on a singular design, but `lm.wfit` does not: it drops the columns it cannot estimate and returns a fit, with a perfectly good residual sum of squares. So `L1` is a real number, usually below `L0`, and R **accepts the step**.

The port does this:

```python
try:
    obj1 = _lstsq_fit(np.column_stack((XREG, U1)), y)
    L1 = obj1[3]
except SegmentedError:
    obj1 = None
    L1 = L0 + 10
```

`_lstsq_fit` raises on rank deficiency, so `L1 = L0 + 10`, the step is judged to have made things worse, and the port **halves it** until the design is full rank again.

This happens when the update throws a breakpoint outside the data. One candidate from BDMA weekly NO₂, h = 0.15:

```
start psi         0.110   0.183   0.595   0.693
psi after update  0.054  -0.217   1.035   1.425      three of four outside [0, 1]
```

The `U` columns for those three are all zero or copies of `z`, so the design has rank 3 of 6.

* R accepts the step, then refits at the final psi, gets `NA` for the `U` coefficients, builds `Vxb = V %*% diag(beta)` full of `NA`, and the closing `lm()` stops with `0 (non-NA) cases`. AQEval's grid search catches the error and discards the candidate.
* Python halves the step back to `0.083, 0.096, 0.705, 0.876`, gets a respectable fit, and keeps it. It has the best adjusted R² of the 81, so it wins. R never considered it.

In that case R fitted 59 of the 81 candidates and Python all 81.

### 2b. A final design that is singular

Before returning, R fits the full model `y ~ 1 + z + U + Vxb` and refuses it if any coefficient is `NA`:

```r
if (isNAcoef) stop("at least one coef is NA: breakpoint(s) at the boundary? (possibly with many x-values replicated)")
```

The port builds the same matrix and inverts `X'X`, relying on `np.linalg.inv` to raise. It does not raise on a matrix that is singular only to rounding error; it returns enormous numbers. One candidate from DYAG daily NO₂, h = 0.2, ends with two breakpoints at 0.1853 and 0.1980 with **no observations between them** (DYAG has a long gap there). The two `Vxb` columns are then multiples of each other, the design has rank 4 of 6 (smallest singular value 1e-20 of the largest), R rejects it, and Python returns it as the winner.

### The fix

[`segmented.patch`](segmented.patch), three small changes:

* score a trial fit by the residual sum of squares of a least-squares fit of whatever rank (`np.linalg.lstsq`), as `lm.wfit` does, instead of `L0 + 10`;
* require the fit kept at the end to be full rank, which is where R's candidate dies;
* check the rank of the final design explicitly instead of trusting `inv`.

With it, the same grid rows succeed in both languages in all three cases logged row by row (71 of 81, 8 of 9 and 59 of 81), their adjusted R² agree to 1e-11, and the same row wins.

One thing not established: R's `lm` calls a column deficient at a relative tolerance of 1e-7, and `matrix_rank` at about 1e-13. In every case seen the deficiency was exact, so the two agree. A design that is nearly but not exactly collinear could still be judged differently.

## Cause 3: where a week starts

For hourly data starting at 00:00 on 1 January 2018:

| | first bin starts |
| --- | --- |
| openair 3.1.0 `timeAverage(avg.time = "7 day")` | Thursday 28 December 2017 |
| `aqeval.time_average(x, "7 day")` | Monday 1 January 2018 |

openair's bins are whole multiples of seven days counted from 1 January 1970, a Thursday. `time_average` anchors at midnight on the first day of the data, which its docstring says is to match openair; that was presumably true of an earlier openair. 8-hour and daily averages are identical, to within summation order (1e-13). Any other multi-day period, such as `"2 day"` or `"14 day"`, is likely to be affected in the same way and was not tested.

## Cause 4: `breaks=None`

`findBreakPoints` returns `NULL` when it finds no breaks, and

```r
breaks <- findBreakPoints(data, "x", h = 0.5)
result <- quantBreakSegments(data, "x", breaks)
```

then fits a single straight line, because `quantBreakSegments` only searches for itself `if (missing(breaks))` and a `NULL` that was passed is not missing. In Python `breaks=None` is the default value, so the same two lines silently run a fresh search with `h = 0.15` and the significance test, and report breaks the caller never asked for. (And, under pandas 3, cause 1 then rejects them all.) A sentinel default would separate "not given" from "none".

## Cause 5: the isolation model

Not a defect. Over the 104 isolation settings the tutorial uses, the series isolated by pygam and by mgcv correlate at 0.991 to 1.000 (median 0.9986), with a root-mean-square difference of 0.08 to 6.3 µg/m³ (median 1.0) and no difference in mean. The largest disagreements are at the extremes of a background term, where the two spline bases extrapolate differently.

It is worth saying plainly in the package's README that this difference, small as it is, is enough to move the changes that `quant_break_segments` reports. Once causes 1 to 4 are fixed it will be the only reason the two languages differ.

## Running the check again

```bash
Rscript pipeline/parity/r_side.R
python pipeline/parity/compare.py                   # as shipped
python pipeline/parity/compare.py --ns              # nanosecond dates
python pipeline/parity/compare.py --own-averaging   # Python averages the hourly data itself
```

It needs the series cached by the main pipeline in `pipeline/cache/iso_r`. Where AQEval cannot be installed, set `AQEVAL_SOURCE` to the `R/` folder of its source, as for the rest of the pipeline. The cases are in `cases.csv`; add rows to widen it. R's outputs in `pipeline/cache/parity` are small CSV files and would serve as test fixtures for the package, so that its tests do not need R.

`compare_grid.py` goes further than the twenty cases: for every series R averaged and measured for the website it gives Python the same series and the same breaks, and compares the fitted trends.

```bash
python pipeline/parity/compare_grid.py
```

Whenever the `aqeval` version changes, the tutorial needs its Python results rebuilt: raise the version in `python/requirements.txt` and `pipeline/requirements.txt`, delete `pipeline/cache/runs/python`, and run `pipeline/build_data.py`.
