"""Does Python ``aqeval`` give the same break analysis as R AQEval?

Run ``r_side.R`` first. This then hands Python exactly the series R analysed,
case by case, and compares the break points, the segments and the fitted
trend with what R returned.

    Rscript pipeline/parity/r_side.R
    python pipeline/parity/compare.py           # dates as aqeval itself produces them
    python pipeline/parity/compare.py --ns      # dates forced to nanosecond resolution
    python pipeline/parity/compare.py --own-averaging   # Python averages the hourly data itself

With aqeval 0.7.3 everything agrees in all three modes. What it found for
0.7.2 is written up in README.md beside this file: the break points agreed
with R in every case but the segments in only 8 of 20, because the package
assumed nanosecond dates (18 of 20 with ``--ns``) and handled rank-deficient
fits differently from R.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import build_data  # noqa: E402 - the pipeline's own loaders

from aqeval import find_break_points, quant_break_segments, time_average  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = build_data.CACHE / "parity"
NO_BREAKS = pd.DataFrame({"lower": [], "bpt": [], "upper": []}, dtype=int)


def read(case: str, what: str) -> pd.DataFrame | None:
    path = OUT / f"{case}_{what}.csv"
    if not path.exists():
        return None
    # round_trip: pandas' default float parser is one unit in the last place
    # out for about one value in eight, and the point here is to hand Python
    # exactly the numbers R had. A missing value in a one-column file is an
    # empty line, which must be kept.
    return pd.read_csv(path, skip_blank_lines=False, float_precision="round_trip")


def hourly(case) -> pd.DataFrame:
    """The hourly series a case starts from, as Python loads it."""
    if case.source == "raw":
        return build_data.load_site(case.site)[["date", case.pollutant]].rename(columns={case.pollutant: "value"})
    frame = pd.read_csv(
        build_data.CACHE / "iso_r" / f"{case.key}.csv", parse_dates=["date"], float_precision="round_trip"
    )
    return frame.rename(columns={"isolated": "value"})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--ns", action="store_true", help="force nanosecond-resolution dates")
    parser.add_argument("--own-averaging", action="store_true", help="let Python average the hourly data itself")
    args = parser.parse_args()
    warnings.simplefilter("ignore")

    rows = []
    for case in pd.read_csv(HERE / "cases.csv", dtype={"key": str}).itertuples():
        given = read(case.id, "input")
        if given is None:
            sys.exit("Run pipeline/parity/r_side.R first.")
        # What aqeval's own functions hand on: second-resolution dates.
        frame = pd.DataFrame({"date": pd.to_datetime(given["epoch"], unit="s").astype("datetime64[s]"), "value": given["value"]})
        input_note = "R's"
        if args.own_averaging:
            own = time_average(hourly(case), case.averaging)[["date", "value"]].reset_index(drop=True)
            if len(own) != len(frame) or not (own["date"].to_numpy() == frame["date"].to_numpy()).all():
                input_note = "DIFFERENT BINS"
            else:
                # Means are summed in a different order, so the last digit or two can differ.
                gap = np.nanmax(np.abs(own["value"].to_numpy() - frame["value"].to_numpy()))
                input_note = "same" if gap < 1e-9 else f"DIFFERS by {gap:.2g}"
            frame = own
        if args.ns:
            frame = frame.assign(date=frame["date"].astype("datetime64[ns]"))

        with contextlib.redirect_stdout(io.StringIO()):
            breaks = find_break_points(frame, "value", h=case.h)
            breaks = NO_BREAKS if breaks is None else breaks
            fitted = quant_break_segments(frame, "value", breaks=breaks, show=()) if len(breaks) <= 3 else None

        r_breaks, r_trend, r_segments = read(case.id, "breaks"), read(case.id, "trend"), read(case.id, "segments")
        row = {
            "case": case.id,
            "series": f"{'isolated' if case.source != 'raw' else 'measured'} {case.site} {case.pollutant}, {case.averaging}, h={case.h}",
            "averaged series": input_note,
            "break points": "same" if np.array_equal(breaks.to_numpy(), r_breaks.to_numpy()) else "DIFFER",
        }
        if fitted is not None and r_trend is not None:
            segments = fitted["segments"]
            if segments is not None and r_segments is not None:
                same = segments.shape == r_segments.shape and np.array_equal(segments.to_numpy(float), r_segments.to_numpy(float), equal_nan=True)
                row["segments"] = "same" if same else "DIFFER"
            ours, theirs = fitted["data2"]["pred"].to_numpy(), r_trend["pred"].to_numpy(float)
            if len(ours) == len(theirs):
                row["largest trend difference"] = float(np.nanmax(np.abs(ours - theirs)))
        rows.append(row)

    table = pd.DataFrame(rows)
    with pd.option_context("display.width", 200, "display.max_rows", 100, "display.float_format", "{:.2g}".format):
        print(table.to_string(index=False))
    measured = table.dropna(subset=["segments"]) if "segments" in table else table.iloc[:0]
    print(f"\nbreak points agree in {(table['break points'] == 'same').sum()} of {len(table)} cases; "
          f"segments in {(measured['segments'] == 'same').sum()} of {len(measured)}")


if __name__ == "__main__":
    main()
