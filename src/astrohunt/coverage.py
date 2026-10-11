"""Which sectors each star has usable QLP light curves for.

Prediction uses tess-point (pointing model); the light-curve download is the
final check. Sector dates come from the official TESS orbit table.

Excluded sectors (see config/sectors.json): 99-104, whose QLP timestamps had a
barycentric-correction error that MAST had not yet redelivered as of
2026-10-10 (https://archive.stsci.edu/hlsp/qlp).
"""

from __future__ import annotations

import io
import json
import re
from multiprocessing import Pool

import numpy as np
import pandas as pd
import requests

from astrohunt import paths

ORBIT_TIMES_URL = "https://tess.mit.edu/public/files/TESS_orbit_times.csv"
QLP_ROOT = "https://archive.stsci.edu/hlsps/qlp/"
MIN_SPAN_DAYS = 365.0


def sector_config() -> dict:
    return json.loads((paths.config_dir() / "sectors.json").read_text())


def sector_times() -> pd.DataFrame:
    """Start, end and mid-time (UTC) per sector from the official orbit table."""
    r = requests.get(ORBIT_TIMES_URL, timeout=60)
    r.raise_for_status()
    lines = [ln for ln in r.text.splitlines() if ln and not ln.startswith("#")]
    df = pd.read_csv(io.StringIO("\n".join(lines)))
    df.columns = [c.strip() for c in df.columns]
    start = pd.to_datetime(df["Start of Orbit"], errors="coerce")
    end = pd.to_datetime(df["End of Orbit"], errors="coerce")
    g = pd.DataFrame({"sector": df["Sector"].astype(int), "start": start, "end": end})
    out = g.groupby("sector").agg(start=("start", "min"), end=("end", "max")).reset_index()
    out["mid"] = out["start"] + (out["end"] - out["start"]) / 2
    return out


def qlp_sectors_available() -> list[int]:
    """Sectors with a directory on the QLP HLSP server, minus excluded ones."""
    r = requests.get(QLP_ROOT, timeout=60)
    r.raise_for_status()
    found = sorted({int(s) for s in re.findall(r'href="s(\d{4})/"', r.text)})
    excluded = set(sector_config()["excluded"])
    return [s for s in found if s not in excluded]


def _stars2px_chunk(args):
    from tess_stars2px import tess_stars2px_function_entry

    ids, ras, decs = args
    out = tess_stars2px_function_entry(ids, ras, decs)
    tid, sec = out[0], out[3]
    keep = tid >= 0
    return pd.DataFrame({"ID": tid[keep].astype(np.int64), "sector": sec[keep].astype(int)})


def predict_sectors(stars: pd.DataFrame, chunk: int = 500, workers: int = 2) -> pd.DataFrame:
    """One row per (star, sector) that tess-point puts on silicon."""
    ids = stars["ID"].astype(np.int64).to_numpy()
    ras = stars["ra"].to_numpy(float)
    decs = stars["dec"].to_numpy(float)
    jobs = [(ids[i:i + chunk], ras[i:i + chunk], decs[i:i + chunk]) for i in range(0, len(ids), chunk)]
    with Pool(workers) as pool:
        parts = pool.map(_stars2px_chunk, jobs)
    return pd.concat(parts, ignore_index=True).drop_duplicates()


def coverage_table(stars: pd.DataFrame, workers: int = 2) -> pd.DataFrame:
    """Per star: usable QLP sectors, count, and time span between first and last."""
    pairs = predict_sectors(stars, workers=workers)
    usable = set(qlp_sectors_available())
    pairs = pairs[pairs["sector"].isin(usable)]
    times = sector_times().set_index("sector")
    pairs = pairs[pairs["sector"].isin(times.index)]
    pairs["mid"] = times.loc[pairs["sector"], "mid"].to_numpy()
    g = pairs.sort_values("sector").groupby("ID")
    cov = pd.DataFrame({
        "sectors": g["sector"].apply(lambda s: ",".join(map(str, s))),
        "n_sectors": g["sector"].size(),
        "span_days": (g["mid"].max() - g["mid"].min()).dt.total_seconds() / 86400.0,
    }).reset_index()
    cov["qualifies"] = (cov["n_sectors"] >= 2) & (cov["span_days"] >= MIN_SPAN_DAYS)
    return cov
