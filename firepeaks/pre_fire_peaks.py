"""take a hydrograph and identify peaks, working in the opposite direction of RREDI """

import numpy as np 
import pandas as pd 
from hashlib import sha1
import pandas as pd
from scipy.signal import find_peaks, argrelextrema
from scipy.stats import pearson3 
from bisect import bisect_left

# some config parameters (tunable)

MIN_SEP_H = 12
REL_SLOPE = 0.05
ABS_SLOPE_FRAC = 0.01   # replaces RREDIS fixed 1 cfs/h, using a fraction of the peak instead (allows for tiny basin response flows)
TROUGH_ORDER_H = 6      # a trough must be the minimum within +/- this many hrs 
MAX_LOOKBACK_H = 120 
HOURS_PER_YEAR = 8766.0  # account for leap years! 
SQKM_PER_SQMI = 2.58999

# define the pre-fire period over which a flood frequency analysis is performed to select a q threshold to filter post-fire flows by
def pre_fire_record(q: pd.Series, ig_date, prior_fires = (), recovery_years: float = 3.0) -> pd.Series: 
    """Flow before 'Ig_date' with the post-fire period from previous fires in the same watershed masked out, such that 
    previous fire's flood do not inflate the threshold"""
    q = q.loc[:pd.Timestamp(ig_date, tz="UTC")].copy() # index to everything up to ig date
    for f in prior_fires: 
        t0 = pd.Timestamp(f, tz = "UTC")
        q.loc[t0:t0 + pd.DateOffset(years = recovery_years)] = np.nan # all flows 3 yrs after a previous fire are rewritten to nan
    return q

def record_years(q: pd.Series) -> float: 
    """record length of the whole pre_fire period sampled in yrs, not calendar gaps so that the masked periods dont count"""
    q = q.sort_index()
    return float(q.resample("1h").mean().notna().sum() / HOURS_PER_YEAR) 

def independence_sep(area_km2: float) -> pd.Timedelta: 
    """Using the standard separation method for declustering: 5 days + ln(area in mi^2)"""
    return pd.Timedelta(days = max(1.0, 5.0 + np.log(area_km2 / SQKM_PER_SQMI))) # fall back to 1 day if lower

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
        neighbours = [kept_pos[j - 1]] if j > 0 else []
        if j < len(kept_pos):
            neighbours.append(kept_pos[j])
        dependent = False
        for nb in neighbours:
            a, b = min(p, nb), max(p, nb)
            if abs(t[p] - t[nb]) < sep_s or trough(a, b) >= trough_ratio * min(v[p], v[nb]):
                dependent = True
                break
        if dependent:
            continue
        kept_pos.insert(j, p)
        kept.append(p)
        if max_keep and len(kept) >= max_keep:
            break
    return q.iloc[kept].rename("peak")
    
def pds_flow(q: pd.Series, area_km2: float, return_period_yr: float = 1.0, min_years: float = 5.0, trough_ratio: float = 0.75) -> dict: 
    """partial duration flow for a return period: the level independent peaks exceed on acverage every 'return_period_yr' years. ie the mth largest peak with m = y/T, interpolated between ranks, some qc info returned w/ peak"""
    yrs = record_years(q)
    m = yrs / return_period_yr
    pk = independent_peaks(q, independence_sep(area_km2), trough_ratio, max_keep=int(np.ceil(m)) + 2)
    x = pk.values
    if yrs < min_years or len(x) < np.ceil(m): 
        val = np.nan
    else: 
        lo = max(int(np.floor(m)), 1)
        hi = min(lo + 1, len(x))
        val = float(np.interp(m, [lo, hi], [x[lo-1], x[hi -1]]))
    return dict(q_pds = val, return_period_yr = return_period_yr, record_years = yrs, n_peaks = len(x), short_record = yrs < min_years) 


# getting the thresholds per gage will look like: 
# q15 = ...
# qh = hourlymax(q15)
# pre = pre_fire_record(qh, ig_date, prior_fires = earlier_ig_dates)
# thr = pds_flow(pre, area_km2)
# 