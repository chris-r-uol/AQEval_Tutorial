"""Precompute every result the interactive tutorial can show.

The website is static (GitHub Pages), so it cannot run AQEval itself. Instead
this script runs the real analysis, with the ``aqeval`` Python package, for
every combination of dials the site offers and writes the answers out as JSON.
Moving a dial on the site looks up the matching answer.

The work is the tutorial's own workflow, repeated:

1. download hourly AURN data for each site;
2. isolate the signal (``isolate_contribution``) at hourly resolution;
3. average (``time_average``);
4. find break points (``find_break_points``) for each value of ``h``;
5. quantify the break segments (``quant_break_segments``).

All of it is done twice, once with Python ``aqeval`` and once with R AQEval
(``isolate.R`` and ``breaks.R``), because the two do not give the same answer
and a student must see on the site what their own code will print. The
isolation model is fitted with different libraries, pygam and mgcv, and the
small difference between their fits is enough to move the changes that are
found.

From aqeval 0.7.3 that is the only difference: given the same series, the
port's averaging, break points and segments match R's (``pipeline/parity``
checks this, including against every R result computed here). Each language's
results still come from that language, so that the site never depends on the
match holding: R isolates, averages and measures the segments for its half.
The one thing borrowed is the search for break points, where R's own takes a
quarter of an hour per series (see ``find_for_r_job``). R's segment fit is the
slow part: allow an hour or two.

Every stage is cached under ``pipeline/cache`` so an interrupted run resumes.

    python pipeline/build_data.py            # everything
    python pipeline/build_data.py --export   # rewrite the JSON from the cache
"""

from __future__ import annotations

import os

# The GAM and the segmented fits are run many at a time, one per process. Left
# alone, each process would also start a thread per core for its linear
# algebra and they would fight each other.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import contextlib
import io
import itertools
import json
import math
import pickle
import shutil
import subprocess
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "pipeline" / "cache"
OUT = ROOT / "site" / "static" / "data"

# --- What the site offers ----------------------------------------------------

YEARS = (2018, 2025)

# The Bradford site, and the nearest AURN stations to compare it with. Bradford
# has only one station with a long NOx record, so the comparison and control
# sites are its neighbours.
SITES = ["BDMA", "LED6", "LEED", "DYAG"]

POLLUTANTS = ["no2", "nox"]

# The languages the tutorial is offered in. Each has its own isolated series.
ENGINES = ["python", "r"]

# What can be handed to ``background=``: nothing, the air temperature column
# (as the original tutorial does), or the same pollutant measured at an urban
# background station.
CONTROL_SITES = ["LEED", "DYAG"]
BACKGROUNDS = [None, "air_temp", *CONTROL_SITES]

AVERAGING = ["8 hour", "day", "7 day"]

# h is the smallest allowed segment, as a fraction of the series. 0.12 is the
# value normally used for 8-hour data and is the most sensitive offered: the
# number of candidate breaks grows as h shrinks, and quantifying k breaks takes
# 9^k model fits.
H_VALUES = [0.12, 0.15, 0.2, 0.3, 0.5]

# 9^4 fits takes half a minute on 8-hour data and 9^5 several minutes, which
# is too long to spend on every set-up. Beyond this the breaks are still
# reported, but not quantified.
MAX_BREAKS_TO_QUANTIFY = 4

EVENTS = [
    {"id": "lockdown", "date": "2020-03-23", "label": "First COVID-19 lockdown", "short": "Lockdown"},
    {"id": "caz", "date": "2022-09-26", "label": "Bradford Clean Air Zone starts", "short": "Clean Air Zone"},
]

BRADFORD_CENTRE = (53.7950, -1.7594)  # City Park

CAZ_BOUNDARY_URL = (
    "https://gis.bradford.gov.uk/server/rest/services/Open_Data/"
    "Clean_Air_Zone_Boundary/MapServer/0/query"
    "?where=1%3D1&outFields=*&outSR=4326&f=geojson"
)


# --- Keys ----------------------------------------------------------------------


def iso_key(deseason: bool, deweather: bool, background: str | None) -> str:
    """Name for one signal-isolation setting, used in file names."""
    return f"ds{int(deseason)}-dw{int(deweather)}-bg{(background or 'none').lower()}"


def iso_settings(site: str):
    """Every isolation setting offered for a site, the 'no isolation' one first."""
    for background, deseason, deweather in itertools.product(BACKGROUNDS, [False, True], [False, True]):
        if background == site:
            continue  # a site cannot be its own control
        yield deseason, deweather, background


# --- Stage 1: download -------------------------------------------------------


def load_site(site: str) -> pd.DataFrame:
    path = CACHE / "raw" / f"{site}.pkl"
    if path.exists():
        return pd.read_pickle(path)
    from aqeval import download_aurn_data

    data = download_aurn_data(site.lower(), YEARS[0], YEARS[1], source="aurn")
    if data.empty:
        raise RuntimeError(f"No AURN data downloaded for {site}")
    path.parent.mkdir(parents=True, exist_ok=True)
    data.to_pickle(path)
    return data


def load_with_background(site: str, pollutant: str, background: str | None) -> tuple[pd.DataFrame, str | None]:
    """The site's data, with the control site's pollutant joined on if asked for.

    Returns the frame and the name of the column to pass as ``background=``.
    """
    data = load_site(site)
    if background in CONTROL_SITES:
        control = load_site(background)[["date", pollutant]].rename(columns={pollutant: "background"})
        joined = data.merge(control, on="date", how="left")
        assert len(joined) == len(data), f"{background} has repeated hours"
        return joined, "background"
    return data, background


# --- Stage 2: signal isolation -------------------------------------------------


def isolate_job(site: str, pollutant: str, deseason: bool, deweather: bool, background: str | None) -> str:
    key = f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}"
    path = CACHE / "iso" / f"{key}.pkl"
    if path.exists():
        return key
    from aqeval import isolate_contribution

    data, column = load_with_background(site, pollutant, background)
    printed = io.StringIO()
    with warnings.catch_warnings(), contextlib.redirect_stdout(printed):
        warnings.simplefilter("ignore")
        isolated = isolate_contribution(data, pollutant, background=column, deseason=deseason, deweather=deweather)
    # aqeval prints the model it fitted, e.g. "fitting: no2 ~ te(wd,ws) + ..."
    formula = printed.getvalue().strip().removeprefix("fitting:").strip()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as handle:
        pickle.dump({"date": data["date"].to_numpy(), "isolated": np.asarray(isolated, dtype=float), "formula": formula}, handle)
    return key


def isolate_in_r(settings: list[tuple], workers: int) -> None:
    """Run the same isolations in R. See ``isolate.R``."""
    folder = CACHE / "iso_r"
    folder.mkdir(parents=True, exist_ok=True)
    jobs = pd.DataFrame(
        {
            "key": f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}",
            "site": site,
            "pollutant": pollutant,
            "deseason": deseason,
            "deweather": deweather,
            # R is given either a column that is already there, or a site to join on.
            "background": None if background in CONTROL_SITES else background,
            "control": background if background in CONTROL_SITES else None,
            "first_year": YEARS[0],
            "last_year": YEARS[1],
        }
        for site, pollutant, deseason, deweather, background in settings
    )
    jobs.to_csv(folder / "jobs.csv", index=False)
    subprocess.run(
        ["Rscript", str(ROOT / "pipeline" / "isolate.R"), str(folder)],
        check=True,
        env={**os.environ, "WORKERS": str(workers)},
    )


# --- Stage 3: break points and segments ---------------------------------------


def iso_date(value) -> str | None:
    if value is None or pd.isna(value):
        return None
    return pd.Timestamp(value).strftime("%Y-%m-%dT%H:%M")


def number(value, digits: int = 2) -> float | None:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return round(value, digits) if math.isfinite(value) else None


def trend_knots(frame: pd.DataFrame) -> list[list]:
    """The fitted trend, reduced to the points where it bends.

    The segmented fit is piecewise linear, so its 8,000-odd predictions carry
    only a handful of numbers. A point is kept if it starts or ends a run of
    data, or if the slope changes there.
    """
    valid = frame.dropna(subset=["pred"])
    if valid.empty:
        return []
    time_ms = valid["date"].astype("datetime64[ms]").astype("int64").to_numpy()
    pred = valid["pred"].to_numpy()
    err = valid["err"].to_numpy()
    keep = np.zeros(len(valid), dtype=bool)
    keep[[0, -1]] = True
    if len(valid) > 2:
        slope = np.diff(pred) / np.diff(time_ms)
        scale = max(float(np.nanmax(np.abs(slope))), 1e-18)
        bends = np.abs(np.diff(slope)) > 1e-6 * scale
        keep[1:-1] |= bends
    # The standard error can come out as NaN where a fit's covariance is
    # singular. JSON has no NaN, so it is written as null.
    return [[int(t), number(p), number(e)] for t, p, e in zip(time_ms[keep], pred[keep], err[keep])]


def segments_from_report(report: pd.DataFrame | None, fitted: pd.DataFrame) -> list[dict]:
    """AQEval's segment report as plain records (it uses R's dotted names).

    The report gives the fitted concentration at each end of each segment by
    looking it up on that row of the data. Where the data has a gap on that
    row — very often the first row of the series — it reports nothing, and the
    size of the change goes with it. The fitted line still has a value there:
    it is straight within a segment, so it is read off the line through the
    segment's own fitted points. ``filled`` records which values came that way.
    """
    if report is None:
        return []
    time_ms = fitted["date"].astype("datetime64[ms]").astype("int64").to_numpy()
    pred = fitted["pred"].to_numpy(dtype=float)
    valid = ~np.isnan(pred)

    rows = []
    for record in report.to_dict("records"):
        duration = record.get("s1.date.delta")
        start, end = record.get("s1.date1"), record.get("s1.date2")
        c0, c1 = float(record.get("s1.c0")), float(record.get("s1.c1"))
        change, percent = record.get("s1.c.delta"), record.get("s1.per.delta")

        filled = []
        if np.isnan(c0) or np.isnan(c1):
            t0 = pd.Timestamp(start).value // 10**6
            t1 = pd.Timestamp(end).value // 10**6
            inside = np.flatnonzero(valid & (time_ms >= t0) & (time_ms <= t1))
            if len(inside) >= 2:
                a, b = inside[0], inside[-1]
                slope = (pred[b] - pred[a]) / (time_ms[b] - time_ms[a])
                if np.isnan(c0):
                    c0 = pred[a] + slope * (t0 - time_ms[a])
                    filled.append("c0")
                if np.isnan(c1):
                    c1 = pred[b] + slope * (t1 - time_ms[b])
                    filled.append("c1")
                change = c1 - c0
                percent = change / c0 * 100 if c0 else None

        rows.append(
            {
                "from": iso_date(start),
                "to": iso_date(end),
                "days": None if pd.isna(duration) else round(pd.Timedelta(duration).total_seconds() / 86400, 1),
                "c0": number(c0),
                "c1": number(c1),
                "change": number(change),
                "percent": number(percent),
                **({"filled": filled} if filled else {}),
            }
        )
    return rows


def series_block(dates: pd.Series, values: np.ndarray, runs: dict) -> dict:
    """One averaging period's results, in the shape the site reads."""
    return {
        "start": iso_date(dates.iloc[0]),
        "stepHours": round((dates.iloc[1] - dates.iloc[0]).total_seconds() / 3600),
        # The series the analysis ran on: isolated, or as measured.
        "values": [None if np.isnan(v) else round(float(v), 1) for v in values],
        "runs": runs,
    }


def write_runs(path: Path, result: dict) -> None:
    """Written under a temporary name first, so an interrupted run never leaves
    a half-written file that looks finished."""
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(".part")
    # allow_nan=False: Python would write NaN, which a browser cannot read.
    partial.write_text(json.dumps(result, separators=(",", ":"), allow_nan=False))
    partial.rename(path)


def breaks_job(site: str, pollutant: str, deseason: bool, deweather: bool, background: str | None) -> str:
    """Python's half: average one series and run the break analysis for every h."""
    key = f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}"
    path = CACHE / "runs" / "python" / f"{key}.json"
    if path.exists():
        return key
    from aqeval import find_break_points, quant_break_segments, time_average

    if deseason or deweather or background is not None:
        with open(CACHE / "iso" / f"{key}.pkl", "rb") as handle:
            stored = pickle.load(handle)
        frame = pd.DataFrame({"date": stored["date"], "value": stored["isolated"]})
        formula = stored["formula"]
    else:
        frame = load_site(site)[["date", pollutant]].rename(columns={pollutant: "value"})
        formula = None

    result = {
        "site": site,
        "pollutant": pollutant,
        "deseason": deseason,
        "deweather": deweather,
        "background": background,
        "formula": formula,
        "averaging": {},
    }

    for averaging in AVERAGING:
        averaged = time_average(frame, averaging)
        dates = averaged["date"].reset_index(drop=True)
        runs: dict[str, dict] = {}
        quantified: dict[tuple, dict] = {}  # several h often find the same breaks

        for h in H_VALUES:
            run: dict = {"breaks": [], "trend": [], "segments": [], "status": "ok"}
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                try:
                    breaks = find_break_points(averaged, "value", h=h)
                except Exception as error:  # noqa: BLE001 - a numerical failure, reported to the reader
                    runs[str(h)] = {**run, "status": "error", "message": str(error)[:200]}
                    continue

                found = 0 if breaks is None else len(breaks)
                if found:
                    # Break tables hold 1-based row numbers, as in R.
                    run["breaks"] = [
                        {
                            "date": iso_date(dates[int(row.bpt) - 1]),
                            "lower": iso_date(dates[int(row.lower) - 1]),
                            "upper": iso_date(dates[int(row.upper) - 1]),
                        }
                        for row in breaks.itertuples()
                    ]

                if found > MAX_BREAKS_TO_QUANTIFY:
                    run["status"] = "too_many_breaks"
                else:
                    if not found:
                        # Handed None, as the tutorial's code hands it, the fit
                        # is a single straight line.
                        run["status"] = "no_breaks"
                    signature = tuple(breaks.to_numpy().ravel()) if found else ()
                    if signature not in quantified:
                        try:
                            fitted = quant_break_segments(averaged, "value", breaks=breaks, show=())
                            quantified[signature] = {
                                "trend": trend_knots(fitted["data2"]),
                                "segments": segments_from_report(fitted["report"], fitted["data2"]),
                            }
                        except Exception as error:  # noqa: BLE001
                            quantified[signature] = {"status": "error", "message": str(error)[:200]}
                    run.update(quantified[signature])
            runs[str(h)] = run

        result["averaging"][averaging] = series_block(dates, averaged["value"].to_numpy(), runs)

    write_runs(path, result)
    return key


def r_stage(stage: str, settings: list[tuple], workers: int) -> None:
    """Run one stage of ``breaks.R`` over every series."""
    folder = CACHE / "breaks_r"
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "key": f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}",
            "site": site,
            "pollutant": pollutant,
            "isolating": "TRUE" if (deseason or deweather or background is not None) else "FALSE",
        }
        for site, pollutant, deseason, deweather, background in settings
    ).to_csv(folder / "jobs.csv", index=False)
    subprocess.run(
        ["Rscript", str(ROOT / "pipeline" / "breaks.R"), str(CACHE / "iso_r"), str(folder)],
        check=True,
        env={**os.environ, "STAGE": stage, "WORKERS": str(workers), "AVERAGING": ";".join(AVERAGING)},
    )


def find_for_r_job(site: str, pollutant: str, deseason: bool, deweather: bool, background: str | None) -> str:
    """Find the break points in the series R averaged, for R to measure.

    This is the one step of R's half that R does not do itself. Its own search
    takes up to a quarter of an hour per series; the port returns the same
    break points and confidence intervals (``pipeline/parity`` checks that) in
    under a second. It works from the values alone, so the defects that affect
    the port's segment fit do not touch it.
    """
    key = f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}"
    path = CACHE / "breaks_r" / f"{key}.breaks.csv"
    if path.exists():
        return key
    from aqeval import find_break_points

    # round_trip: pandas' default float parser is one unit in the last place out
    # for about one value in eight, and these must be exactly R's numbers.
    averaged = pd.read_csv(CACHE / "breaks_r" / f"{key}.averaged.csv", float_precision="round_trip")
    rows = []
    for period, series in averaged.groupby("period", sort=False):
        frame = pd.DataFrame({"date": pd.to_datetime(series["epoch"], unit="s"), "value": series["value"].to_numpy()})
        for h in H_VALUES:
            with warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
                warnings.simplefilter("ignore")
                breaks = find_break_points(frame, "value", h=h)
            if breaks is None or len(breaks) == 0:
                rows.append({"period": period, "h": h, "status": "no_breaks"})
                continue
            status = "too_many_breaks" if len(breaks) > MAX_BREAKS_TO_QUANTIFY else "ok"
            rows.extend({"period": period, "h": h, "status": status, **row} for row in breaks.to_dict("records"))

    partial = path.with_suffix(".part")
    pd.DataFrame(rows, columns=["period", "h", "status", "lower", "bpt", "upper"]).astype(
        {"lower": "Int64", "bpt": "Int64", "upper": "Int64"}
    ).to_csv(partial, index=False)
    partial.rename(path)
    return key


def floats(values: list) -> np.ndarray:
    """A JSON array with nulls, as numbers with NaN."""
    return np.array([np.nan if v is None else v for v in values], dtype=float)


def r_runs_job(site: str, pollutant: str, deseason: bool, deweather: bool, background: str | None) -> str:
    """Turn what ``breaks.R`` wrote for one series into the site's format.

    R hands over the series it averaged, the breaks it was given, and for each
    distinct set of breaks the fitted trend and the report that
    ``quantBreakSegments`` returned. From there the two languages go through
    the same code, so that a figure looks the same whichever produced it.
    """
    key = f"{site}_{pollutant}_{iso_key(deseason, deweather, background)}"
    raw = json.loads((CACHE / "breaks_r" / f"{key}.json").read_text())
    isolating = deseason or deweather or background is not None

    result = {
        "site": site,
        "pollutant": pollutant,
        "deseason": deseason,
        "deweather": deweather,
        "background": background,
        "formula": (CACHE / "iso_r" / f"{key}.txt").read_text().strip() if isolating else None,
        "averaging": {},
    }

    for averaging, block in raw["averaging"].items():
        dates = pd.Series(pd.to_datetime(block["epoch"], unit="s"))

        fits = []
        for fit in block["fits"]:
            if "error" in fit:
                fits.append({"status": "error", "message": fit["error"]})
                continue
            fitted = pd.DataFrame({"date": dates, "pred": floats(fit["pred"]), "err": floats(fit["err"])})
            report = None
            if fit.get("report"):
                r = fit["report"]
                report = pd.DataFrame(
                    {
                        "s1.date1": pd.to_datetime(r["date1"], unit="s"),
                        "s1.date2": pd.to_datetime(r["date2"], unit="s"),
                        "s1.date.delta": pd.to_timedelta(r["days"], unit="D"),
                        "s1.c0": floats(r["c0"]),
                        "s1.c1": floats(r["c1"]),
                        "s1.c.delta": floats(r["change"]),
                        "s1.per.delta": floats(r["percent"]),
                    }
                )
            fits.append({"trend": trend_knots(fitted), "segments": segments_from_report(report, fitted)})

        runs = {}
        for h, run in block["runs"].items():
            converted = {
                # Break tables hold 1-based row numbers.
                "breaks": [
                    {"date": iso_date(dates[bpt - 1]), "lower": iso_date(dates[lower - 1]), "upper": iso_date(dates[upper - 1])}
                    for lower, bpt, upper in run["breaks"]
                ],
                "trend": [],
                "segments": [],
                "status": run["status"],
            }
            if run.get("message"):
                converted["message"] = run["message"]
            if run.get("fit") is not None:
                converted.update(fits[run["fit"] - 1])  # R counts from one
            runs[h] = converted

        result["averaging"][averaging] = series_block(dates, floats(block["values"]), runs)

    write_runs(CACHE / "runs" / "r" / f"{key}.json", result)
    return key


# --- Stage 4: export -----------------------------------------------------------


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    p = math.pi / 180
    a = math.sin((lat2 - lat1) * p / 2) ** 2 + math.cos(lat1 * p) * math.cos(lat2 * p) * math.sin((lon2 - lon1) * p / 2) ** 2
    return 2 * 6371 * math.asin(math.sqrt(a))


def caz_boundary() -> dict:
    """Bradford's Clean Air Zone boundary, simplified for the web.

    Published by Bradford Council under the Open Government Licence v3.0. The
    original has 12,000 vertices; a few metres of tolerance leaves a few
    hundred and no visible difference at city scale.
    """
    import requests
    from shapely.geometry import mapping, shape

    path = CACHE / "caz_raw.geojson"
    if not path.exists():
        response = requests.get(CAZ_BOUNDARY_URL, timeout=60)
        response.raise_for_status()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(response.content)
    polygon = shape(json.loads(path.read_text())["features"][0]["geometry"])
    simplified = polygon.simplify(0.00005, preserve_topology=True)  # about 4 m
    geometry = json.loads(json.dumps(mapping(simplified)))
    geometry["coordinates"] = [[[round(x, 5), round(y, 5)] for x, y in ring] for ring in geometry["coordinates"]]
    return {
        "type": "Feature",
        "properties": {
            "name": "Bradford Clean Air Zone",
            "source": "City of Bradford Metropolitan District Council",
            "licence": "Open Government Licence v3.0",
        },
        "geometry": geometry,
    }


def export() -> None:
    """Write everything the site reads."""
    from shapely.geometry import Point, shape

    import aqeval
    from aqeval import import_aq_meta

    shutil.rmtree(OUT, ignore_errors=True)
    OUT.mkdir(parents=True)

    boundary = caz_boundary()
    (OUT / "caz.geojson").write_text(json.dumps(boundary, separators=(",", ":")))
    zone = shape(boundary["geometry"])

    meta = import_aq_meta("aurn")
    meta.columns = [str(c) for c in meta.columns]
    launch = pd.Timestamp(next(e["date"] for e in EVENTS if e["id"] == "caz"))

    sites = []
    for code in SITES:
        row = meta[meta["site_id"] == code].iloc[0]
        lat, lon = float(row["latitude"]), float(row["longitude"])
        data = load_site(code)
        year = data["date"].dt.year

        stats = {}
        for pollutant in POLLUTANTS:
            before = data.loc[(data["date"] >= launch - pd.DateOffset(years=1)) & (data["date"] < launch), pollutant]
            after = data.loc[(data["date"] >= launch) & (data["date"] < launch + pd.DateOffset(years=1)), pollutant]
            stats[pollutant] = {
                "annual": {str(y): number(v, 1) for y, v in data.groupby(year)[pollutant].mean().items()},
                # Share of hours with a valid measurement, per year.
                "capture": {str(y): number(v * 100, 0) for y, v in data.groupby(year)[pollutant].apply(lambda s: s.notna().mean()).items()},
                "before": number(before.mean(), 1),
                "after": number(after.mean(), 1),
                "mean": number(data[pollutant].mean(), 1),
            }

        sites.append(
            {
                "code": code,
                "name": str(row["site_name"]),
                "type": str(row["location_type"]),
                "authority": str(row["local_authority"]),
                "lat": lat,
                "lon": lon,
                "insideCaz": bool(zone.contains(Point(lon, lat))),
                "kmFromBradford": round(haversine_km(*BRADFORD_CENTRE, lat, lon), 1),
                "isolations": [iso_key(*setting) for setting in iso_settings(code)],
                "stats": stats,
            }
        )

    offered = {str(h) for h in H_VALUES}
    for engine in ENGINES:
        (OUT / "runs" / engine).mkdir(parents=True)
        for source in sorted((CACHE / "runs" / engine).glob("*.json")):
            result = json.loads(source.read_text())
            # The cache may hold values of h that are no longer offered.
            for block in result["averaging"].values():
                block["runs"] = {h: run for h, run in block["runs"].items() if h in offered}
            (OUT / "runs" / engine / source.name).write_text(
                json.dumps(result, separators=(",", ":"), allow_nan=False)
            )

    versions = CACHE / "iso_r" / "versions.txt"
    manifest = {
        "generated": time.strftime("%Y-%m-%d"),
        "engines": {
            "python": f"aqeval {aqeval.__version__}",
            "r": versions.read_text().strip() if versions.exists() else "AQEval",
        },
        "years": list(YEARS),
        "sites": sites,
        "pollutants": POLLUTANTS,
        "controlSites": CONTROL_SITES,
        "averaging": AVERAGING,
        "h": H_VALUES,
        "maxBreaksToQuantify": MAX_BREAKS_TO_QUANTIFY,
        "events": EVENTS,
        "centre": list(BRADFORD_CENTRE),
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=1, allow_nan=False))
    size = sum(f.stat().st_size for f in OUT.rglob("*") if f.is_file())
    print(f"wrote {sum(1 for _ in OUT.rglob('*.json'))} files, {size / 1e6:.1f} MB, to {OUT.relative_to(ROOT)}")


# --- Driver ---------------------------------------------------------------------


def run_parallel(function, jobs: list[tuple], workers: int, label: str) -> None:
    started = time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(function, *job): job for job in jobs}
        for done, future in enumerate(as_completed(futures), start=1):
            future.result()  # re-raise anything that went wrong in the worker
            print(f"{label} {done}/{len(jobs)}  {time.time() - started:,.0f}s  {futures[future]}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--export", action="store_true", help="only rewrite the JSON from the cache")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 2))
    args = parser.parse_args()

    if not args.export:
        for site in SITES:
            frame = load_site(site)
            print(f"{site}: {len(frame):,} hourly rows, {frame['date'].min():%Y-%m-%d} to {frame['date'].max():%Y-%m-%d}")

        settings = [(site, pollutant, *setting) for site in SITES for pollutant in POLLUTANTS for setting in iso_settings(site)]
        # (False, False, None) is 'no isolation': there is no model to fit.
        isolating = [s for s in settings if any(s[2:])]
        run_parallel(isolate_job, isolating, args.workers, "isolate")
        isolate_in_r(isolating, args.workers)
        run_parallel(breaks_job, settings, args.workers, "breaks (Python)")
        r_stage("average", settings, args.workers)
        run_parallel(find_for_r_job, settings, args.workers, "find breaks (for R)")
        r_stage("quantify", settings, args.workers)
        run_parallel(r_runs_job, settings, args.workers, "segments (R)")

    export()


if __name__ == "__main__":
    sys.exit(main())
