"""Per-candidate vetting gauntlet (plan: "Vetting gauntlet").

Implemented here:
- test 4, hide-a-season: fit the ephemeris without one segment, predict that
  segment's transits, measure them there.
- test 5, depth consistency: per-sector depths, and per photometry version
  (QLP switched to TGLC photometry from Sector 94).
- test 7, Gaia companion indicators: RUWE and RV-scatter flags (Gaia DR3).

Tests 1-3 come from LEO-Vetter (odd/even, secondary, shape); test 6 (pixel
localisation) runs LEO's difference-image check. TRICERATOPS (8) is separate. Pass rules are the
plan's draft values, recalibrated where noted; current: config/gauntlet_v2.json.
"""

from __future__ import annotations

import io
import json

import numpy as np
import pandas as pd
import requests

from astrohunt import paths
from astrohunt.search import Prepared, Signal, coherent_refine

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def rules() -> dict:
    return json.loads((paths.config_dir() / "gauntlet_v2.json").read_text())


def template(metrics: dict):
    """LEO's fitted transit model (batman, quadratic limb darkening) as m(t)."""
    from leo_vetter.models import TransitModel

    tm = TransitModel(metrics["transit_per"], metrics["transit_epo"], metrics["transit_RpRs"],
                      metrics["transit_aRs"], metrics["transit_b"], metrics["transit_u1"],
                      metrics["transit_u2"], 1.0)
    return tm


def _model(tm, t, per=None, epo=None):
    params = tm.params.copy()
    if per is not None:
        params["per"].value = per
    if epo is not None:
        params["epo"].value = epo
    return tm.model(params, t)


def fit_scale(p: Prepared, sel: np.ndarray, tm, per=None, epo=None) -> tuple[float, float, int]:
    """Least-squares depth scale of the template in the selected cadences.

    Returns (depth in relative flux, 1-sigma error, n in-transit cadences).
    The template's own depth sets the unit; a = 1 means 'as deep as the global fit'.
    """
    m = _model(tm, p.t[sel], per, epo)
    s = 1.0 - m                          # template dip, >= 0
    intr = s > 0
    n = int(intr.sum())
    if n < 3:
        return np.nan, np.nan, n
    w = 1.0 / p.e[sel] ** 2
    resid = 1.0 - p.f[sel]
    denom = float(np.sum(w * s * s))
    a = float(np.sum(w * resid * s) / denom)
    depth0 = float(s.max())
    return a * depth0, depth0 / np.sqrt(denom), n


def empirical_err(p: Prepared, sel: np.ndarray, tm, per=None, epo=None, n: int = 60,
                  seed: int = 0) -> float:
    """Depth error from the template fitted at shifted phases away from real transits.

    Captures red noise that the white-noise error misses. Shifts are drawn in
    [0.15, 0.85] of the period, which keeps the shifted windows off the real
    transits for any duration below 0.15 P.
    """
    per = tm.params["per"].value if per is None else per
    epo = tm.params["epo"].value if epo is None else epo
    rng = np.random.default_rng(seed)
    vals = []
    for f in rng.uniform(0.15, 0.85, n):
        d, _, k = fit_scale(p, sel, tm, per, epo + f * per)
        if k >= 3 and np.isfinite(d):
            vals.append(d)
    if len(vals) < 10:
        return np.nan
    v = np.asarray(vals)
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def refit_ephemeris(p: Prepared, sel: np.ndarray, tm, per0: float, epo0: float, dur: float,
                    n_grid: int = 121) -> dict:
    """Chi-square scan in (epoch, period), template shape fixed.

    The epoch is pivoted to the transit nearest the middle of the selected
    data, so epoch and period are nearly uncorrelated. Returns the best
    values and every grid point inside delta-chi2 < 1, from which any
    predicted transit time gets its uncertainty.
    """
    near = sel & (np.abs(((p.t - epo0 + 0.5 * per0) % per0) - 0.5 * per0) < 3 * dur)
    t, f, e = p.t[near], p.f[near], p.e[near]
    tp0 = epo0 + round((np.median(t) - epo0) / per0) * per0
    span = float(t.max() - t.min())
    dP = dur / max(span / per0, 1.0)                 # period change that shifts the ends by one duration
    dts = np.linspace(-0.5 * dur, 0.5 * dur, n_grid)
    dps = np.linspace(-dP, dP, n_grid)
    chi2 = np.empty((n_grid, n_grid))
    for i, dt in enumerate(dts):
        for j, dp in enumerate(dps):
            m = _model(tm, t, per0 + dp, tp0 + dt)
            chi2[i, j] = np.sum(((f - m) / e) ** 2)
    i, j = np.unravel_index(np.argmin(chi2), chi2.shape)
    ok = chi2 - chi2[i, j] < 1.0
    ii, jj = np.nonzero(ok)
    edge = bool(i in (0, n_grid - 1) or j in (0, n_grid - 1))
    return {"per": float(per0 + dps[j]), "tp": float(tp0 + dts[i]), "edge": edge,
            "dt_step": float(dts[1] - dts[0]), "dp_step": float(dps[1] - dps[0]),
            "region": [(float(tp0 + dts[a]), float(per0 + dps[b])) for a, b in zip(ii, jj)]}


def predicted_sigma(fit: dict, target: float) -> float:
    """Half-range of predicted transit times near `target` over the delta-chi2 < 1 region."""
    n = round((target - fit["tp"]) / fit["per"])
    times = [tp + n * per for tp, per in fit["region"]]
    # grid resolution at the target is a floor (one grid cell may hold the whole
    # delta-chi2 < 1 region when the surface is steep)
    floor = 0.5 * float(np.hypot(fit["dt_step"], n * fit["dp_step"]))
    return 0.5 * (max(times) - min(times)) + floor


def hide_a_season(p: Prepared, sig: Signal, metrics: dict) -> list[dict]:
    """Hold out each segment; ephemeris from the rest; measure the held-out depth.

    'insufficient' when: the rest leaves the cycle count ambiguous (BLS alias
    ratio > max_alias_ratio); the predicted time is uncertain by more than a
    quarter of the duration; or the held-out segment has no predicted in-transit
    data.
    """
    r = rules()["hide_a_season"]
    tm = template(metrics)
    dur = sig.duration
    out = []
    for k, idx in enumerate(p.segments):
        held = np.zeros(p.t.size, bool)
        held[idx] = True
        rest = ~held
        row = {"segment": k, "sectors": sorted({int(s) for s in np.unique(p.sector[idx])})}
        bls = coherent_refine(p, rest, sig.period, sig.duration)
        if bls is None:
            out.append({**row, "result": "insufficient", "why": "no BLS fit without this segment"})
            continue
        row["alias_ratio"] = bls.alias_ratio
        if bls.alias_ratio > r["max_alias_ratio"]:
            out.append({**row, "result": "insufficient", "why": "cycle count ambiguous without this segment"})
            continue
        fit = refit_ephemeris(p, rest, tm, bls.period, metrics["transit_epo"], dur)
        per, epo = fit["per"], fit["tp"]
        tmid = float(np.median(p.t[idx]))
        sig_t = predicted_sigma(fit, tmid)
        row.update({"P_rest": per, "sigma_t_pred_min": sig_t * 1440, "grid_edge": fit["edge"]})
        if fit["edge"] or not np.isfinite(sig_t) or sig_t > 0.25 * dur:
            out.append({**row, "result": "insufficient", "why": "predicted transit time too uncertain"})
            continue
        d_h, e_h, n_h = fit_scale(p, held, tm, per, epo)
        d_r, e_r, _ = fit_scale(p, rest, tm, per, epo)
        e_h = float(np.nanmax([e_h, empirical_err(p, held, tm, per, epo)]))
        e_r = float(np.nanmax([e_r, empirical_err(p, rest, tm, per, epo)]))
        # propagate the prediction's own uncertainty: spread of the held-out depth
        # over ephemerides inside the delta-chi2 < 1 region (up to 50 of them)
        reg = fit["region"]
        if len(reg) > 1 and np.isfinite(d_h):
            pick = np.linspace(0, len(reg) - 1, min(50, len(reg))).astype(int)
            spread = np.nanstd([fit_scale(p, held, tm, reg[q][1], reg[q][0])[0] for q in pick])
            e_h = float(np.hypot(e_h, spread))
            row["ephemeris_spread_ppm"] = spread * 1e6
        row.update({"n_in_held": n_h, "depth_held_ppm": d_h * 1e6, "err_held_ppm": e_h * 1e6,
                    "depth_rest_ppm": d_r * 1e6, "err_rest_ppm": e_r * 1e6})
        if n_h < 3 or not np.isfinite(d_h):
            out.append({**row, "result": "insufficient", "why": "no predicted in-transit data"})
            continue
        z = abs(d_h - d_r) / np.hypot(e_h, e_r)
        row.update(snr_held=d_h / e_h, depth_z=z)
        row["result"] = "pass" if (d_h / e_h >= r["min_snr"] and z <= r["max_depth_z"]) else "fail"
        out.append(row)
    return out


def depth_consistency(p: Prepared, metrics: dict, tglc_from: int = 94) -> dict:
    """Per-sector template depths; chi-square within each photometry version."""
    from scipy.stats import chi2 as chi2d

    r = rules()["depth_consistency"]
    tm = template(metrics)
    rows = []
    for s in np.unique(p.sector):
        d, e, n = fit_scale(p, p.sector == s, tm)
        e = float(np.nanmax([e, empirical_err(p, p.sector == s, tm)]))
        if n >= 3 and np.isfinite(d):
            rows.append({"sector": int(s), "version": "TGLC" if s >= tglc_from else "pre_TGLC",
                         "depth_ppm": d * 1e6, "err_ppm": e * 1e6, "n_in": n})
    df = pd.DataFrame(rows)
    res = {"per_sector": rows, "by_version": {}}
    verdicts = []
    for v, g in df.groupby("version") if len(df) else []:
        w = 1 / g["err_ppm"] ** 2
        mean = float((w * g["depth_ppm"]).sum() / w.sum())
        err = float(1 / np.sqrt(w.sum()))
        entry = {"mean_ppm": mean, "err_ppm": err, "n_sectors": len(g)}
        if len(g) >= 2:
            chi2 = float((((g["depth_ppm"] - mean) / g["err_ppm"]) ** 2).sum())
            pval = float(chi2d.sf(chi2, len(g) - 1))
            entry.update(chi2=chi2, p_value=pval)
            verdicts.append(pval >= r["min_p_value"])
        res["by_version"][v] = entry
    bv = res["by_version"]
    if "pre_TGLC" in bv and "TGLC" in bv:
        a, b = bv["pre_TGLC"], bv["TGLC"]
        res["cross_version_z"] = abs(a["mean_ppm"] - b["mean_ppm"]) / np.hypot(a["err_ppm"], b["err_ppm"])
    res["result"] = "insufficient" if not verdicts else ("pass" if all(verdicts) else "fail")
    return res


def gaia_flags(dr2_source_id: int) -> dict:
    """Gaia DR3 indicators of an unresolved companion, via the DR2->DR3 crossmatch."""
    r = rules()["gaia"]
    q = ("SELECT TOP 1 g.source_id, g.ruwe, g.ipd_frac_multi_peak, g.phot_bp_rp_excess_factor, "
         "g.radial_velocity_error, g.rv_nb_transits, g.rv_chisq_pvalue, g.rv_renormalised_gof, "
         "g.non_single_star FROM gaiadr3.dr2_neighbourhood AS n JOIN gaiadr3.gaia_source AS g "
         f"ON g.source_id = n.dr3_source_id WHERE n.dr2_source_id = {int(dr2_source_id)} "
         "ORDER BY n.angular_distance ASC")
    resp = requests.get(GAIA_TAP, params={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q},
                        timeout=120)
    resp.raise_for_status()
    df = pd.read_csv(io.StringIO(resp.text))
    if df.empty:
        return {"result": "insufficient", "why": "no DR3 match"}
    g = df.iloc[0].to_dict()
    ruwe_flag = bool(g["ruwe"] > r["ruwe_max"])
    rv_flag = bool(pd.notna(g["rv_chisq_pvalue"]) and g["rv_chisq_pvalue"] < r["rv_pvalue_max"]
                   and pd.notna(g["rv_renormalised_gof"]) and g["rv_renormalised_gof"] > r["rv_gof_min"]
                   and (g["rv_nb_transits"] or 0) >= r["rv_min_transits"])
    nss_flag = bool(g.get("non_single_star", 0) and g["non_single_star"] > 0)
    flags = [n for n, f in (("RUWE", ruwe_flag), ("RV_variable", rv_flag), ("non_single_star", nss_flag)) if f]
    return {**{k: (None if pd.isna(v) else v) for k, v in g.items()}, "flags": flags,
            "result": "flag" if flags else "pass"}


def run_candidates(run: str, log=print) -> pd.DataFrame:
    """Gauntlet tests 4, 5 and 7 for every row of data/runs/<run>/shortlist_new.csv."""
    from astrohunt import lightcurves as L
    from astrohunt import search as S
    from astrohunt import vetting as V

    d = paths.sub("runs", run)
    cands = pd.read_csv(d / "shortlist_new.csv")
    stars = pd.read_parquet(paths.sub("catalogs") / "tic_dwarfs.parquet")
    stars["ID"] = stars["ID"].astype(np.int64)
    stars = stars.set_index("ID")
    out_dir = paths.sub("runs", run, "gauntlet")
    summary = []
    for _, c in cands.iterrows():
        tic = int(c["tic"])
        st = stars.loc[tic]
        lcs = [lc for lc in (L.load(tic, int(s)) for s in str(c["sectors"]).split(",")) if lc is not None]
        p = S.prepare(lcs)
        sigs = S.search(p)
        match = [s for s in sigs if abs(s.period / c["period"] - 1) < 1e-3]
        if not match:
            summary.append({"tic": tic, "period": c["period"], "result": "not reproduced"})
            continue
        sig = match[0]
        p = S.prepare(lcs, masks=sigs)      # transits masked in the detrending
        m = V.leo_vet(p, sig, V.star_dict(st), 1, others=[o for o in sigs if o is not sig])
        res = {"tic": tic, "period": sig.period, "leo": {k: m.get(k) for k in
               ("MES", "FA_fail", "FP_fail", "transit_per", "transit_epo", "transit_RpRs", "transit_b")},
               "hide_a_season": hide_a_season(p, sig, m),
               "depth_consistency": depth_consistency(p, m)}
        # test 6: pixel localisation (difference images from TESScut, LEO PRF centroid)
        try:
            secs = sorted({int(s) for s in np.unique(p.sector[S.in_transit(p.t, sig, 1.0)])})
            res["pixel"] = V.pixel_vet(p, sig, V.star_dict(st), float(st["ra"]), float(st["dec"]), secs,
                                       str(paths.sub("tdi")), others=[o for o in sigs if o is not sig])
            res["pixel"]["result"] = "fail" if res["pixel"]["FP_fail"] else (
                "pass" if res["pixel"]["sectors_used"] else "insufficient")
        except Exception as exc:
            res["pixel"] = {"result": "error", "error": f"{type(exc).__name__}: {exc}"}
        try:
            res["gaia"] = gaia_flags(int(st["GAIA"])) if pd.notna(st["GAIA"]) else {"result": "no Gaia id"}
        except Exception as exc:
            res["gaia"] = {"result": "error", "error": str(exc)}
        (out_dir / f"tic{tic}_p{sig.period:.5f}.json").write_text(json.dumps(res, indent=1, default=float))
        hs = [h["result"] for h in res["hide_a_season"]]
        summary.append({"tic": tic, "period": sig.period, "hide_pass": hs.count("pass"),
                        "hide_fail": hs.count("fail"), "hide_insufficient": hs.count("insufficient"),
                        "depth": res["depth_consistency"]["result"], "pixel": res["pixel"].get("result"),
                        "offset_arcsec": res["pixel"].get("offset_qual_arcsec"), "gaia": res["gaia"].get("result"),
                        "gaia_flags": ",".join(res["gaia"].get("flags", []))})
        log(f"TIC {tic} P={sig.period:.5f}: hide {hs}, depth {summary[-1]['depth']}, gaia {summary[-1]['gaia']}")
    df = pd.DataFrame(summary)
    df.to_csv(d / "gauntlet_summary.csv", index=False)
    return df
