"""Does Python ``aqeval`` reproduce every R result the website shows?

``compare.py`` checks twenty cases in detail. This checks all of them. For
each series that R averaged and measured for the website (``breaks.R``),
Python is given the same averaged series and the same break points, and its
fitted trend is compared with R's at every point.

    python pipeline/parity/compare_grid.py

It needs the main pipeline to have run, and reads its cache. Nothing is
written. With aqeval 0.7.3, 1,329 of 1,332 fits agreed to within 1e-6 µg/m³
and the other three by about that much; a different choice of fit is out by
whole units.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import build_data  # noqa: E402 - the pipeline's own settings and paths

FOLDER = build_data.CACHE / "breaks_r"

# R's trend comes through JSON with 15 significant figures, so agreement is
# judged to a little short of that. A different choice of fit is out by whole
# units, not by the ninth decimal place.
TOLERANCE = 1e-6


def check(key: str) -> list[dict]:
    """Compare every fit for one series. One record per distinct set of breaks."""
    from aqeval import quant_break_segments

    warnings.simplefilter("ignore")
    averaged = pd.read_csv(FOLDER / f"{key}.averaged.csv", float_precision="round_trip")
    searches = pd.read_csv(FOLDER / f"{key}.breaks.csv")
    r_results = json.loads((FOLDER / f"{key}.json").read_text())["averaging"]

    records = []
    for period, series in averaged.groupby("period", sort=False):
        frame = pd.DataFrame({"date": pd.to_datetime(series["epoch"], unit="s"), "value": series["value"].to_numpy()})
        seen = set()
        for h, found in searches[searches["period"] == period].groupby("h", sort=False):
            run = r_results[period]["runs"][str(h)]
            if run["fit"] is None or run["fit"] in seen:
                continue  # too many breaks to measure, or this set of breaks is already checked
            seen.add(run["fit"])
            r_fit = r_results[period]["fits"][run["fit"] - 1]
            breaks = None if found["status"].iloc[0] == "no_breaks" else found[["lower", "bpt", "upper"]].astype(int)

            record = {"key": key, "period": period, "h": h, "breaks": 0 if breaks is None else len(breaks)}
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    fitted = quant_break_segments(frame, "value", breaks=breaks, show=())
                ours = fitted["data2"]["pred"].to_numpy()
                error = None
            except Exception as failure:  # noqa: BLE001 - R failing on the same fit counts as agreement
                ours, error = None, str(failure)[:80]

            if "error" in r_fit:
                record["outcome"] = "both fail" if error else "only R fails"
            elif error:
                record["outcome"] = "only Python fails"
            else:
                theirs = np.array([np.nan if v is None else v for v in r_fit["pred"]], dtype=float)
                same_gaps = np.array_equal(np.isnan(ours), np.isnan(theirs))
                record["gap"] = float(np.nanmax(np.abs(ours - theirs))) if same_gaps else float("inf")
                record["outcome"] = "same" if record["gap"] < TOLERANCE else "DIFFERENT"
            records.append(record)
    return records


def main() -> None:
    keys = sorted(path.name[: -len(".json")] for path in FOLDER.glob("*.json"))
    if not keys:
        sys.exit("Nothing to check: run pipeline/build_data.py first.")
    import aqeval

    with ProcessPoolExecutor(max_workers=max(1, (os.cpu_count() or 2) - 2)) as pool:
        records = [record for batch in pool.map(check, keys) for record in batch]
    table = pd.DataFrame(records)

    print(f"aqeval {aqeval.__version__}: {len(table)} fits across {len(keys)} series")
    print(table.groupby(["period", "breaks", "outcome"]).size().unstack(fill_value=0).to_string())
    agreed = table["outcome"].isin(["same", "both fail"])
    print(f"\nPython matches R in {agreed.sum()} of {len(table)} fits")
    if "gap" in table and (table["outcome"] == "same").any():
        print(f"largest difference between matching trends: {table.loc[table['outcome'] == 'same', 'gap'].max():.1e}")
    if not agreed.all():
        worst = table[~agreed].sort_values("gap", ascending=False).head(10)
        print("\nfits that differ:")
        print(worst.to_string(index=False))


if __name__ == "__main__":
    main()
