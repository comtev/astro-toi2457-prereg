"""Multi-sector transit search.

Why this shape: one coherent BLS over a multi-year baseline needs a period
grid so fine (frequency step ~ duration / baseline^2) that it is impractical
on CPUs. Instead:

1. Detrend each sector (wotan biweight) and split the data into segments of
   at most two consecutive sectors.
2. Run BLS (log-likelihood objective) on each segment over one shared period
   grid, sized so a transit cannot drift by more than EPS of the shortest
   duration across the longest segment.
3. Sum the segment log-likelihoods: an incoherent stack that lets each
   segment have its own epoch.
4. Refine the strongest peaks coherently over the full baseline on a narrow,
   fine grid. This also resolves the cycle-count ambiguity between years.
5. Keep the best signal, mask its transits, repeat (up to MAX_SIGNALS).

Numbers below are the pipeline-v0 settings; the injection tests decide
whether they change before the freeze.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.timeseries import BoxLeastSquares
from wotan import flatten

from astrohunt.lightcurves import SectorLC

P_MIN, P_MAX = 0.5, 15.0                     # days
# Period chunks with the transit durations physically plausible for K/M
# dwarfs (stellar density ~1-20 solar; T ~ 13 h (P/365 d)^(1/3) rho^(-1/3)
# for a central transit), five durations each, ~1.45x apart.
CHUNKS = [
    (0.5, 2.0, np.array([0.02, 0.03, 0.045, 0.065, 0.10])),
    (2.0, 6.0, np.array([0.03, 0.045, 0.065, 0.09, 0.14])),
    (6.0, 15.0, np.array([0.04, 0.06, 0.085, 0.12, 0.19])),
]
DURATIONS = np.unique(np.concatenate([c[2] for c in CHUNKS]))
OVERSAMPLE = 5                               # BLS phase bins per shortest duration
POLISH_OVERSAMPLE = 30                       # final ephemeris pass (~3-4 min phase bins)
EPS = 0.25                                   # max transit drift across a segment, in units of min duration
WINDOW = 1.0                                 # wotan biweight window, days
CLIP_SIGMA = 4.0                             # upper-only clipping (flares, cosmics)
SEGMENT_GAP = 10.0                           # days; a larger gap starts a new segment
SEGMENT_MAX_SECTORS = 2
N_PEAKS = 5                                  # incoherent peaks refined per iteration
MAX_SIGNALS = 3
KEEP_SNR = 7.0                               # coherent BLS depth S/N to keep a signal for vetting
BIN_DAYS = 10.0 / 1440                       # bin faster cadences to 10 min (QLP's Year 3-4 cadence)


def bin_to(t: np.ndarray, f: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray]:
    """Mean-bin to `width` if the native cadence is finer; empty bins dropped."""
    if t.size < 2 or np.median(np.diff(t)) > 0.9 * width:
        return t, f
    idx = np.floor((t - t[0]) / width).astype(np.int64)
    n = np.bincount(idx)
    keep = n > 0
    tb = np.bincount(idx, weights=t)[keep] / n[keep]
    fb = np.bincount(idx, weights=f)[keep] / n[keep]
    return tb, fb


@dataclass
class Prepared:
    tic: int
    t: np.ndarray
    f: np.ndarray           # detrended, clipped
    e: np.ndarray
    raw: np.ndarray         # normalised, undetrended (for LEO's 'raw')
    sector: np.ndarray
    segments: list[np.ndarray] = field(default_factory=list)

    @property
    def baseline(self) -> float:
        return float(self.t.max() - self.t.min())


def detrend_sector(lc: SectorLC, window: float = WINDOW, model=None, masks=None):
    """Return t, detrended flux, raw flux, per-point error for one sector.

    `model`, if given, is a callable m(t) multiplied into the raw flux before
    detrending: that is how injections see the full processing chain.
    """
    g = lc.good()
    t, f = lc.time[g], lc.flux[g].copy()
    if model is not None:
        f = f * model(t)
    if t.size < 100:
        return None
    t, f = bin_to(t, f, BIN_DAYS)
    tmask = np.zeros(t.size, bool)
    for sig in masks or []:
        # 1.5 durations: covers ingress/egress plus ephemeris error
        tmask |= in_transit(t, sig, width=1.5)
    # Split only at gaps > 0.5 d (orbit/sector breaks) and cut 0.25 d at those
    # edges. Rabtsevich's 0.1 d break tolerance split at every flagged block and
    # cut 0.25 d at each, deleting ~4 days of QLP Sector 6 around TIC 4206066's
    # transits; no edge cut at all let a ramp at a Sector 98 gap become the
    # strongest "signal". Chosen on TIC 4206066 (2026-10-11), see commit message.
    trend = flatten(t, f, method="biweight", window_length=window, break_tolerance=0.5,
                    edge_cutoff=0.25, return_trend=True, mask=tmask if tmask.any() else None)[1]
    fd = f / trend
    ok = np.isfinite(fd)
    t, fd, f = t[ok], fd[ok], f[ok]
    if t.size < 100:
        return None
    sigma = 1.4826 * np.median(np.abs(fd - np.median(fd)))
    if tmask.any():
        sigma = 1.4826 * np.median(np.abs(fd[~tmask[ok]] - np.median(fd[~tmask[ok]])))
    keep = fd < 1 + CLIP_SIGMA * sigma
    t, fd, f = t[keep], fd[keep], f[keep]
    return t, fd, f, np.full_like(fd, sigma)


def make_segments(t: np.ndarray, sector: np.ndarray) -> list[np.ndarray]:
    secs = sorted(np.unique(sector), key=lambda s: t[sector == s].min())
    groups, cur, last_end = [], [], None
    for s in secs:
        ts = t[sector == s]
        if cur and (ts.min() - last_end > SEGMENT_GAP or len(cur) >= SEGMENT_MAX_SECTORS):
            groups.append(cur)
            cur = []
        cur.append(s)
        last_end = ts.max()
    if cur:
        groups.append(cur)
    return [np.flatnonzero(np.isin(sector, g)) for g in groups]


def prepare(lcs: list[SectorLC], window: float = WINDOW, model=None, masks=None) -> Prepared | None:
    """Detrend and stack sectors. `masks`: signals whose transits the trend must ignore
    (used after the search, so vetting and depths see undistorted transits)."""
    parts = []
    for lc in lcs:
        out = detrend_sector(lc, window, model, masks)
        if out is not None:
            parts.append((*out, np.full(out[0].size, lc.sector)))
    if not parts:
        return None
    t, f, raw, e, s = (np.concatenate(x) for x in zip(*parts))
    o = np.argsort(t)
    p = Prepared(int(lcs[0].tic), t[o], f[o], e[o], raw[o], s[o])
    p.segments = make_segments(p.t, p.sector)
    return p


def durations_for(period: float) -> np.ndarray:
    for lo, hi, durs in CHUNKS:
        if lo <= period < hi:
            return durs
    return CHUNKS[-1][2]


def period_grid(t: np.ndarray, segments: list[np.ndarray]) -> list[np.ndarray]:
    """One ascending period array per chunk; fractional step EPS*dmin/T_segment."""
    tseg = max(t[idx].max() - t[idx].min() for idx in segments)
    grids = []
    for lo, hi, durs in CHUNKS:
        step = EPS * durs.min() / tseg
        n = int(np.ceil(np.log(hi / lo) / np.log1p(step)))
        freqs = (1 / hi) * (1 + step) ** np.arange(n + 1)
        grids.append(np.sort(1 / freqs[(freqs <= 1 / lo) & (freqs > 1 / hi * 0.999999)]))
    return grids


def incoherent_periodogram(p: Prepared, mask: np.ndarray, grids: list[np.ndarray]):
    """Summed segment log-likelihood per period, chunk by chunk (each chunk its own durations)."""
    out = []
    for (lo, hi, durs), periods in zip(CHUNKS, grids):
        total = np.zeros(periods.size)
        for idx in p.segments:
            idx = idx[mask[idx]]
            if idx.size < 200:
                continue
            bls = BoxLeastSquares(p.t[idx], p.f[idx], p.e[idx])
            r = bls.power(periods, durs, objective="likelihood", method="fast", oversample=OVERSAMPLE)
            total += np.nan_to_num(r.power, nan=0.0)
        out.append(total)
    return out


def peak_score(power: np.ndarray, width: int = 1001) -> np.ndarray:
    """Log-likelihood above a running median (removes the slow red-noise rise).

    Kept in log-likelihood units on purpose: dividing by a local scatter
    (SDE-style) blew up at chunk edges in testing (TIC 4206066, 2026-10-10).
    """
    from scipy.ndimage import median_filter

    w = max(3, min(width, (power.size // 2) * 2 - 1))
    return power - median_filter(power, size=w, mode="reflect")


def distinct_peaks(periods: np.ndarray, score: np.ndarray, k: int, tol: float = 0.01) -> list[int]:
    picked: list[int] = []
    for i in np.argsort(np.nan_to_num(score, nan=-np.inf))[::-1]:
        if all(abs(periods[i] - periods[j]) / periods[j] > tol for j in picked):
            picked.append(int(i))
        if len(picked) == k:
            break
    return picked


@dataclass
class Signal:
    period: float
    t0: float
    duration: float
    depth: float
    depth_snr: float
    loglike: float
    peak_score: float       # incoherent log-likelihood above running median
    n_transits: int
    n_sectors_with_transit: int
    alias_ratio: float      # second-best / best coherent log-likelihood near this period


def transit_epochs(p: Prepared, mask: np.ndarray, period: float, t0: float, dur: float):
    ph = ((p.t - t0 + 0.5 * period) % period) - 0.5 * period
    intr = (np.abs(ph) < 0.5 * dur) & mask
    epochs = np.round((p.t[intr] - t0) / period)
    return np.unique(epochs), np.unique(p.sector[intr])


def coherent_refine(p: Prepared, mask: np.ndarray, period0: float, dur0: float) -> Signal | None:
    tseg = np.median([p.t[i].max() - p.t[i].min() for i in p.segments])
    half = 3.0 * period0 * dur0 / tseg            # ~3x the single-segment period uncertainty
    half = min(half, 0.2 * period0)               # short segments: keep the window positive
                                                  # (TIC 300642809 inverted run, pilot_v0)
    step = EPS * dur0 / p.baseline
    n = int(min(20000, np.ceil(2 * half / (period0 * step))))
    periods = np.linspace(period0 - half, period0 + half, max(n, 11))
    durs = np.unique(np.clip(dur0 * np.array([0.7, 0.85, 1.0, 1.2, 1.45]), 0.02, 0.5 * periods.min()))
    bls = BoxLeastSquares(p.t[mask], p.f[mask], p.e[mask])
    r = bls.power(periods, durs, objective="likelihood", method="fast", oversample=OVERSAMPLE)
    pw = np.nan_to_num(r.power, nan=-np.inf)
    i = int(np.argmax(pw))
    if not np.isfinite(pw[i]) or pw[i] <= 0:
        return None
    # alias check: best power outside +-1 cycle-count step of the winner
    cyc = periods[i] ** 2 / p.baseline
    far = np.abs(periods - periods[i]) > 0.5 * cyc
    alias = float(pw[far].max() / pw[i]) if far.any() and np.isfinite(pw[far].max()) else 0.0
    # polish: finer period steps and phase bins around the winner (BLS 'fast'
    # bins phase at duration/OVERSAMPLE, ~20 min, too coarse for per-sector depths)
    P1 = float(periods[i])
    dstep = (periods[1] - periods[0]) if periods.size > 1 else P1 * step
    fine = np.linspace(P1 - 2 * dstep, P1 + 2 * dstep, 41)
    fdurs = np.unique(np.clip(float(r.duration[i]) * np.array([0.85, 1.0, 1.15]), 0.02, 0.5 * fine.min()))
    r2 = bls.power(fine, fdurs, objective="likelihood", method="fast", oversample=POLISH_OVERSAMPLE)
    pw2 = np.nan_to_num(r2.power, nan=-np.inf)
    j = int(np.argmax(pw2))
    if pw2[j] >= pw[i]:
        r, pw, periods, i = r2, pw2, fine, j
    P, T0, D = float(periods[i]), float(r.transit_time[i]), float(r.duration[i])
    ep, secs = transit_epochs(p, mask, P, T0, D)
    return Signal(P, T0, D, float(r.depth[i]), float(r.depth_snr[i]), float(pw[i]), np.nan,
                  int(ep.size), int(secs.size), alias)


def search(p: Prepared, max_signals: int = MAX_SIGNALS, grids: list[np.ndarray] | None = None,
           keep_snr: float = KEEP_SNR) -> list[Signal]:
    grids = period_grid(p.t, p.segments) if grids is None else grids
    mask = np.ones(p.t.size, bool)
    found: list[Signal] = []
    for _ in range(max_signals):
        powers = incoherent_periodogram(p, mask, grids)
        # SDE within each chunk (different durations -> different baselines), then pooled
        periods = np.concatenate(grids)
        score = np.concatenate([peak_score(pw) for pw in powers])
        best: Signal | None = None
        for i in distinct_peaks(periods, score, N_PEAKS):
            # refine the peak and its 2x and 1/2x harmonics; keep the best likelihood
            for h in (1.0, 2.0, 0.5):
                ph_ = periods[i] * h
                if not (P_MIN <= ph_ <= P_MAX):
                    continue
                dur0 = _best_duration(p, mask, ph_)
                s = coherent_refine(p, mask, ph_, dur0)
                if s is not None:
                    s.peak_score = float(score[i])
                    if best is None or s.loglike > best.loglike:
                        best = s
        if best is None or best.depth_snr < keep_snr or best.depth <= 0:
            break
        found.append(best)
        mask &= ~in_transit(p.t, best, width=2.0)
    return found


def in_transit(t: np.ndarray, sig: "Signal", width: float = 1.0) -> np.ndarray:
    """Cadences within +-width/2 durations of a transit centre."""
    ph = ((t - sig.t0 + 0.5 * sig.period) % sig.period) - 0.5 * sig.period
    return np.abs(ph) < 0.5 * width * sig.duration


def _best_duration(p: Prepared, mask: np.ndarray, period: float) -> float:
    durs = durations_for(period)
    total = np.zeros(durs.size)
    for idx in p.segments:
        idx = idx[mask[idx]]
        if idx.size < 200:
            continue
        bls = BoxLeastSquares(p.t[idx], p.f[idx], p.e[idx])
        # duration-resolved power at a single period: one call per duration
        for k, d in enumerate(durs):
            r = bls.power([period], [d], objective="likelihood", method="fast", oversample=OVERSAMPLE)
            total[k] += np.nan_to_num(r.power[0], nan=0.0)
    return float(durs[int(np.argmax(total))])
