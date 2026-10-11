"""Download QLP light curves and store them as compact per-sector arrays.

File naming follows the QLP README
(https://archive.stsci.edu/hlsps/qlp/hlsp_qlp_tess_ffi_all_tess_v1_readme.txt):
hlsp_qlp_tess_ffi_s<sector:04d>-<tic:016d>_tess_v<ver>_llc.fits, under
s<sector>/<tic[0:4]>/<tic[4:8]>/<tic[8:12]>/<tic[12:16]>/.

Two column layouts exist:
- older sectors: SAP_FLUX (raw), KSPSAP_FLUX (QLP spline-detrended), KSPSAP_FLUX_ERR
- newer sectors: SAP_FLUX, SYS_RM_FLUX (systematics removed, not detrended),
  DET_FLUX (QLP-detrended), DET_FLUX_ERR

We search on the least-processed systematics-corrected flux (SYS_RM_FLUX,
else SAP_FLUX) and do our own detrending, so injection tests see the same
processing as real signals. QLP's own detrended flux is kept for cross-checks.
"""

from __future__ import annotations

import csv
import hashlib
import io
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits

from astrohunt import paths
from astrohunt.coverage import sector_config

QLP_ROOT = "https://archive.stsci.edu/hlsps/qlp/"

# Quality bits treated as bad. Rabtsevich's leo_rerun.py mask (QLP-update
# recommendations plus SPOC/TICA attitude, safe-mode, Earth-point, manual
# exclude, and the 2^29 scattered-light flag), plus 1024 (cosmic ray in
# collateral pixels) from the older-format bit list in the QLP data product
# description.
BAD_BITS = 1 | 2 | 4 | 8 | 16 | 32 | 128 | 1024 | 2048 | 4096 | (1 << 29) | (1 << 30)

MANIFEST_FIELDS = ["tic", "sector", "status", "url", "sha256", "bytes", "retrieved_utc",
                   "flux_column", "n_cadences", "n_good"]


def qlp_url(tic: int, sector: int) -> str:
    ver = "v02" if sector in set(sector_config()["v02_sectors"]) else "v01"
    t = f"{int(tic):016d}"
    return (f"{QLP_ROOT}s{sector:04d}/{t[0:4]}/{t[4:8]}/{t[8:12]}/{t[12:16]}/"
            f"hlsp_qlp_tess_ffi_s{sector:04d}-{t}_tess_{ver}_llc.fits")


@dataclass
class SectorLC:
    tic: int
    sector: int
    time: np.ndarray        # BTJD (BJD - 2457000), TDB
    flux: np.ndarray        # normalised, not detrended
    flux_err: np.ndarray
    quality: np.ndarray
    qlp_det: np.ndarray     # QLP's own detrended flux, for cross-checks
    flux_column: str

    def good(self) -> np.ndarray:
        return (np.isfinite(self.time) & np.isfinite(self.flux)
                & ((self.quality & BAD_BITS) == 0))


def parse_fits(data: bytes, tic: int, sector: int) -> SectorLC:
    with fits.open(io.BytesIO(data)) as h:
        d = h[1].data
        names = set(d.columns.names)
        if "SYS_RM_FLUX" in names:
            flux_col, err_col, det_col = "SYS_RM_FLUX", "DET_FLUX_ERR", "DET_FLUX"
        else:
            flux_col, err_col, det_col = "SAP_FLUX", "KSPSAP_FLUX_ERR", "KSPSAP_FLUX"
        t = np.asarray(d["TIME"], float)
        f = np.asarray(d[flux_col], float)
        fe = np.asarray(d[err_col], float)
        q = np.asarray(d["QUALITY"], np.int64)
        det = np.asarray(d[det_col], float)
    med = np.nanmedian(f[(q & BAD_BITS) == 0]) if np.any((q & BAD_BITS) == 0) else np.nanmedian(f)
    return SectorLC(int(tic), int(sector), t, f / med, fe / med, q, det, flux_col)


def lc_path(tic: int, sector: int) -> Path:
    return paths.sub("lc", f"{int(tic) % 1000:03d}") / f"tic{int(tic)}_s{sector:04d}.npz"


def save(lc: SectorLC) -> Path:
    p = lc_path(lc.tic, lc.sector)
    np.savez_compressed(p, time=lc.time, flux=lc.flux, flux_err=lc.flux_err,
                        quality=lc.quality, qlp_det=lc.qlp_det,
                        meta=np.array([lc.tic, lc.sector]), flux_column=np.array(lc.flux_column))
    return p


def load(tic: int, sector: int) -> SectorLC | None:
    p = lc_path(tic, sector)
    if not p.exists():
        return None
    z = np.load(p)
    return SectorLC(int(tic), int(sector), z["time"], z["flux"], z["flux_err"], z["quality"],
                    z["qlp_det"], str(z["flux_column"]))


def fetch(tic: int, sector: int, session: requests.Session | None = None,
          keep_fits: bool = False, retries: int = 3) -> dict:
    """Download one light curve; return a manifest row. Never raises on 404."""
    url = qlp_url(tic, sector)
    row = {"tic": int(tic), "sector": int(sector), "url": url, "status": "",
           "sha256": "", "bytes": 0, "retrieved_utc": "", "flux_column": "",
           "n_cadences": 0, "n_good": 0}
    if lc_path(tic, sector).exists():
        row["status"] = "cached"
        return row
    s = session or requests.Session()
    for attempt in range(retries):
        try:
            r = s.get(url, timeout=120)
            if r.status_code == 404:
                row["status"] = "missing"
                return row
            r.raise_for_status()
            break
        except requests.RequestException as exc:
            if attempt == retries - 1:
                row["status"] = f"error:{type(exc).__name__}"
                return row
            time.sleep(5 * (attempt + 1))
    data = r.content
    row.update(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data),
               retrieved_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"))
    try:
        lc = parse_fits(data, tic, sector)
    except Exception as exc:
        row["status"] = f"corrupt:{type(exc).__name__}"
        return row
    save(lc)
    if keep_fits:
        (paths.sub("fits") / url.rsplit("/", 1)[-1]).write_bytes(data)
    row.update(status="ok", flux_column=lc.flux_column, n_cadences=int(lc.time.size),
               n_good=int(lc.good().sum()))
    return row


def append_manifest(rows: list[dict]) -> None:
    p = paths.sub("lc") / "manifest.csv"
    new = not p.exists()
    with p.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=MANIFEST_FIELDS)
        if new:
            w.writeheader()
        for row in rows:
            if row["status"] != "cached":
                w.writerow(row)


def load_spoc(tic: int) -> tuple[list[SectorLC], list[dict]]:
    """SPOC 2-min PDCSAP light curves for a TIC from MAST (cached under data/spoc/<tic>).

    Independent photometry for cross-checks; not used by the search. PDCSAP is already
    corrected for crowding (CROWDSAP) and aperture losses (FLFRCSAP). Flux errors are left
    NaN, as for QLP, so search.prepare estimates noise from the scatter. Quality flags are
    reduced to good (0) / bad (1).
    Returns (light curves, per-sector metadata).
    """
    import glob

    from astroquery.mast import Observations

    obs = Observations.query_criteria(target_name=str(tic), obs_collection="TESS",
                                      dataproduct_type="timeseries", provenance_name="SPOC")
    if not len(obs):
        return [], []
    prods = Observations.get_product_list(obs)
    prods = prods[[str(f).endswith("-s_lc.fits") for f in prods["productFilename"]]]
    ddir = paths.sub("spoc", str(tic))
    Observations.download_products(prods, download_dir=str(ddir), cache=True)
    lcs, meta = [], []
    for f in sorted(glob.glob(str(ddir / "**" / f"*{tic:016d}*-s_lc.fits"), recursive=True)):
        with fits.open(f) as h:
            d, sec = h[1].data, int(h[0].header["SECTOR"])
            hdr = h[1].header
            t = np.asarray(d["TIME"], float)
            fl = np.asarray(d["PDCSAP_FLUX"], float)
            q = np.asarray(d["QUALITY"], np.int64)
            meta.append({"sector": sec, "CROWDSAP": hdr.get("CROWDSAP"), "FLFRCSAP": hdr.get("FLFRCSAP"),
                         "file": Path(f).name})
        ok = np.isfinite(t) & np.isfinite(fl)
        med = np.nanmedian(fl[ok & (q == 0)])
        lcs.append(SectorLC(tic, sec, t, fl / med, np.full(t.size, np.nan), np.where(q == 0, 0, 1),
                            fl / med, "PDCSAP_FLUX"))
    return lcs, meta
