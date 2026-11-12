"""take a hydrograph and identify peaks, working in the opposite direction of RREDI, then set a thresh from pre-fire flow """

import numpy as np 
import pandas as pd 
from hashlib import sha1
import pandas as pd
from scipy.signal import find_peaks, argrelextrema
from scipy.stats import pearson3 
from bisect import bisect_left

# some config parameters (tunable)

MIN_SEP_H = 12          #two peaks are independent only if they are >= MIN_SEP_H hours apart
TROUGH_RATIO = 0.75      # AND flow between this drops below TROUGH_RATIO * smaller peak
POT_RATE = 5            # threshold = flow independent peaks exceed POT_RATE times/yr
TROUGH_ORDER_H = 6      # a trough must be the minimum within +/- this many hrs 
MIN_YEARS = 5.0          # shorter pre-fire records get no threshold 
HOURS_PER_YEAR = 8766.0  # account for leap years! 
SQKM_PER_SQMI = 2.58999

def utc(t) -> pd.Timestamp: 
    t = pd.Timestamp(t)
    return t.tz_localize("UTC") if t.tz is None else t.tz_convert("UTC")

# define the pre-fire period over which a flood frequency analysis is performed to select a q threshold to filter post-fire flows by
def pre_fire_record(q: pd.Series, ig_date, prior_fires = (), recovery_years: int = 3) -> pd.Series: 
    """Flow before 'Ig_date' with the post-fire period from previous fires in the same watershed masked out, such that 
    previous fire's flood do not inflate the threshold"""
    q = q.loc[:pd.Timestamp(ig_date, tz="UTC")].copy() # index to everything up to ig date
    for f in prior_fires: 
        t0 = utc(f)
        q.loc[t0:t0 + pd.DateOffset(years = recovery_years)] = np.nan # all flows 3 yrs after a previous fire are rewritten to nan
    return q

def record_years(q: pd.Series) -> float: 
    """record length of the whole pre_fire period sampled in yrs, not calendar gaps so that the masked periods dont count"""
    q = q.sort_index()
    return float(q.resample("1h").mean().notna().sum() / HOURS_PER_YEAR) 

def _range_min(v: np.ndarray):
    """Sparse table for O(1) minimum over v[a:b+1], so the trough test stays fast on decades of hourly flow"""
    table = [v]
    k = 1
    while 2 * k <= len(v):
        prev = table[-1]
        table.append(np.minimum(prev[:-k], prev[k:]))
        k *= 2
    def query(a: int, b: int) -> float:
        lvl = int(np.log2(b - a + 1))
        return float(min(table[lvl][a], table[lvl][b - (1 << lvl) + 1]))
    return query

def independent_peaks(q: pd.Series, sep: pd.Timedelta, trough_ratio: float = 0.75, max_keep: int | None = None) -> pd.Series:
    """Independent peaks, greedy for largest first, following USWRC (1976) Bulletin 17: a candidate is kept only if it is at least
    'sep' from its nearest kept neighbours on both sides (see independence_sep()) AND flow between it and each of those neighbours
    drops below 'trough_ratio' * the smaller of the two peaks. Returned largest first."""
    q = q.dropna().sort_index()
    if q.empty:
        return q.rename("peak")
    v = q.to_numpy(dtype=float)
    t = ((q.index - q.index[0]) / pd.Timedelta(seconds=1)).to_numpy() # seconds since start, convert timestamps to handle irregulr intervals
    sep_s = sep / pd.Timedelta(seconds=1)
    trough = _range_min(v)
    idx, _ = find_peaks(v) # every local maximum is a candidate; independence is decided below
    kept: list[int] = []      # positions, in the order kept (largest first)
    kept_pos: list[int] = []  # the same positions, kept sorted by time
    for p in idx[np.argsort(v[idx], kind="stable")[::-1]]:  # sort then slice to reverse list
        j = bisect_left(kept_pos, p) # get the index where to insert candidate peak in kept positions
        neighbours = ([kept_pos[j - 1]] if j > 0 else []) + ([kept_pos[j]] if j < len(kept_pos) else [])
        if any(abs(t[p] - t[nb]) < sep_s or trough(min(p, nb), max(p, nb)) >= trough_ratio * min(v[p], v[nb])
               for nb in neighbours):
            continue
        kept_pos.insert(j, p)
        kept.append(p)
        if max_keep and len(kept) >= max_keep:
            break
    return q.iloc[kept].rename("peak")
        
    
def pot_threshold(qpre: pd.Series, rate: float = POT_RATE, min_sep_h = MIN_SEP_H, 
                  trough_ratio: float = TROUGH_RATIO, min_years: int = MIN_YEARS) -> dict: 
    """ peaks over threshod for nominal rate (rate), pass in entire pre-fire  series
    """
    y = record_years(qpre)
    m = y * rate
    pk = independent_peaks(qpre, pd.Timedelta(hours = min_sep_h), trough_ratio, max_keep=int(np.ceil(m)) + 2)
    x = pk[pk > 0 ].to_numpy()
    if y < min_years or len(x) < np.ceil(m): 
        val = np.nan
    else: 
        lo = max(int(np.floor(m)), 1)
        hi = min(lo + 1, len(x))
        val = float(np.interp(m, [lo, hi], [x[lo-1], x[hi -1]]))
    return dict(qthr = val, rate = rate, record_years = y, short_record = y < min_years) # populates short record boolean if short


def melt_season(q: pd.Series) -> dict: 
    """ 2 harmonic fit of daily mean log flow on DOY (pre-fire record), returns melt_doys = the run of days where 
    around the fitted maximum where the fit is in the top half of its range (the second harmonic can raise monsoon so we try and avoid that)"""
    d = q.dropna().resample("1D").mean().dropna()
    if len(d) < 365:                    # the record is less than 1 yr, nothing to fit 
        return dict(h_r2 = np.nan, h_peak_hoy = np.nan, melt_doys = set())   # get the empty dict 
    y = np.log(d.clip(lower = max))
