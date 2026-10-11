"""Run a frozen pre-registered transit test on new TESS sectors. Written before the data exist.

Usage: uv run python scripts/prereg_test.py PREREG_JSON SECTORS [--dry-run]
       e.g. uv run python scripts/prereg_test.py docs/preregistration/toi2457_p3.2911.json 109,110

Everything is fixed by PREREG_JSON (from scripts/preregister.py): template shape, ephemeris,
prior depth, detrending masks. Nothing is refitted except where a rule says so. Rules are in
the Markdown file next to the JSON; this script implements them:

Data:  SPOC 2-min PDCSAP for SECTORS if MAST has it; otherwise QLP (same loader as the search).
       Detrending: search.prepare with the frozen masks; transits of the other masked signals removed.
Coverage: a predicted transit counts when >= 50% of its expected in-transit cadences exist
       (quality-good, after detrending). Fewer than 8 counted transits -> INCONCLUSIVE.
T1 detection:   depth d at the frozen ephemeris; e = max(white, 200-phase empirical error);
                pass if d / e >= 3.
T2 consistency: |d - prior| / hypot(e, prior_err) <= 3.
T3 cycle count: in the new data alone, the frozen alias k=0 has a lower (better) chi-square
                improvement over a flat line than k=-1 and k=+1, and that improvement is < 0.
Verdict: CONFIRMED if T1, T2, T3 pass; NOT CONFIRMED if any fails with >= 8 counted transits.
Reported, not judged: free-epoch timing offset (template refit on a grid, +-0.5 duration).
--dry-run: allow sectors already used in the prior (pipeline check only; verdict not meaningful).
"""

import json
import sys
from pathlib import Path

import numpy as np

from astrohunt import gauntlet as G
from astrohunt import lightcurves as L
from astrohunt import search as S

args = [a for a in sys.argv[1:] if not a.startswith("--")]
dry = "--dry-run" in sys.argv
pre = json.loads(Path(args[0]).read_text())
sectors = [int(x) for x in args[1].split(",")]
used = set(pre["inputs"]["spoc_sectors"]) | set(pre["inputs"]["qlp_sectors"])
if used & set(sectors) and not dry:
    sys.exit(f"sectors {sorted(used & set(sectors))} were used to build the prediction; not a test")
tic = pre["tic"]

tm = G.template(pre["template"])
per, tp = pre["ephemeris"]["period_d"], pre["ephemeris"]["t_mid_btjd"]
dur = pre["template"]["transit_dur"]
masks = [S.Signal(m["period"], m["t0"], m["duration"], 0, 0, 0, 0, 0, 0, 0) for m in pre["inputs"]["detrend_masks"]]
others = [s for s, m in zip(masks, pre["inputs"]["detrend_masks"]) if not m["is_candidate"]]

slcs, meta = L.load_spoc(tic)
slcs = [lc for lc in slcs if lc.sector in sectors]
source = "SPOC 2-min PDCSAP"
missing = sorted(set(sectors) - {lc.sector for lc in slcs})
if missing:
    for s in missing:
        L.fetch(tic, s)
    q = [lc for lc in (L.load(tic, s) for s in missing) if lc is not None]
    slcs += q
    source += f"; QLP for {[lc.sector for lc in q]}" if q else ""
res = {"prereg": args[0], "sectors": sectors, "data_source": source, "dry_run": dry,
       "sectors_with_data": sorted({lc.sector for lc in slcs})}
if not slcs:
    res["verdict"] = "INCONCLUSIVE"
    res["why"] = "no data for the requested sectors yet"
    print(json.dumps(res, indent=1))
    sys.exit(0)

p = S.prepare(slcs, masks=masks)
keep = np.ones(p.t.size, bool)
for o in others:
    keep &= ~S.in_transit(p.t, o, 1.5)

# coverage
cad = np.median(np.diff(np.sort(p.t)))
n_lo, n_hi = int(np.floor((p.t.min() - tp) / per)), int(np.ceil((p.t.max() - tp) / per))
counted = []
for n in range(n_lo, n_hi + 1):
    tc = tp + n * per
    got = int(np.sum(keep & (np.abs(p.t - tc) < 0.5 * dur)))
    if got >= 0.5 * dur / cad:
        counted.append(round(tc, 4))
res["cadence_min"] = cad * 1440
res["n_transits_counted"] = len(counted)
res["transits_counted_btjd"] = counted

# T1, T2
d, e_w, n_in = G.fit_scale(p, keep, tm, per, tp)
e_emp = G.empirical_err(p, keep, tm, per, tp, n=200)
e = float(np.nanmax([e_w, e_emp]))
prior, prior_e = pre["prior_depth_ppm"] * 1e-6, pre["prior_depth_err_ppm"] * 1e-6
res.update({"depth_ppm": d * 1e6, "err_ppm": e * 1e6, "white_err_ppm": e_w * 1e6, "empirical_err_ppm": e_emp * 1e6,
            "snr": d / e, "T1_pass": bool(d / e >= 3)})
z = abs(d - prior) / np.hypot(e, prior_e)
res.update({"depth_z_vs_prior": float(z), "T2_pass": bool(z <= 3)})


# T3: aliases k=-1,0,+1 from the frozen alias table, chi-square improvement over flat on the new data
def improvement(per_k, tp_k):
    near = keep & (np.abs(((p.t - tp_k + 0.5 * per_k) % per_k) - 0.5 * per_k) < 3 * dur)
    mod = G._model(tm, p.t[near], per_k, tp_k)
    w = 1 / p.e[near] ** 2
    flat = np.sum(w * p.f[near]) / np.sum(w)
    return float(np.sum(((p.f[near] - mod) / p.e[near]) ** 2) - np.sum(((p.f[near] - flat) / p.e[near]) ** 2))


al = {a["k"]: a for a in pre["aliases"]}
imp = {k: improvement(al[k]["period"], al[k]["tp"]) for k in (-1, 0, 1)}
res["alias_chi2_improvement"] = imp
res["T3_pass"] = bool(imp[0] < 0 and imp[0] < imp[-1] and imp[0] < imp[1])

# reported only: free-epoch timing offset
dts = np.linspace(-0.5 * dur, 0.5 * dur, 121)
near = keep & (np.abs(((p.t - tp + 0.5 * per) % per) - 0.5 * per) < 3 * dur)
c2 = [np.sum(((p.f[near] - G._model(tm, p.t[near], per, tp + dt)) / p.e[near]) ** 2) for dt in dts]
i = int(np.argmin(c2))
inside = dts[np.array(c2) - c2[i] < 1]
res["timing_offset_min"] = float(dts[i] * 1440)
res["timing_offset_range_min"] = [float(inside.min() * 1440), float(inside.max() * 1440)]

if len(counted) < 8:
    res["verdict"] = "INCONCLUSIVE"
    res["why"] = f"only {len(counted)} predicted transits covered (< 8)"
elif res["T1_pass"] and res["T2_pass"] and res["T3_pass"]:
    res["verdict"] = "CONFIRMED"
else:
    res["verdict"] = "NOT CONFIRMED"
if dry:
    res["verdict"] += " (dry run: sectors already in the prior; not a test)"
out = Path(args[0]).with_name(Path(args[0]).stem + f"_result_s{'-'.join(map(str, sectors))}{'_dryrun' if dry else ''}.json")
out.write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
print(out)
