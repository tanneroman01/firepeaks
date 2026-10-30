import numpy as np 
import pandas as pd 
from scipy.signal import find_peaks, argrelextrema

# some config parameters (tunable)
MIN_SEP_H = 12
REL_SLOPE = 0.05
ABS_SLOPE_FRAC = 0.01   # replaces RREDIS fixed 1 cfs/h, using a fraction of the peak instead (allows for tiny basin response flows)
TROUGH_ORDER_H = 6      # a trough must be the minimum within +/- this many hrs
MAX_LOOKBACK_H = 120




# define the post fire period over which we search the hydrograph 
def post_fire_period(ig_date, years: float = 3.0): 
    t0 = pd.Timestamp(ig_date, tz="UTC") 
    return t0, t0 + pd.DateOffset(years=years)


# pull the 
def hourly(q15: pd.Series) -> pd.Series: 
    """hourly mean of 15-min flow -- because AORC is hourly, we need to locate the rise start on hourly
    """
    return q15.resample("1h").mean()

def hourly_max(q15: pd.Series) -> pd.Series: 
    """Hourly MAXIMUM of 15 minute flow, the hourly() function shaves flashy peaks so we use it for rise detection only"""
    return q15.resample("1h").max()



def find_events(q15: pd.Series, threshold: float, start, end, min_sep_h: float = MIN_SEP_H): 
    """Independent 15 min peaks above 'threshold' within [start, end]"""
    q = q15.loc[start:end].dropna()
    step_h = q.index.to_series().diff().median() / pd.Timedelta(hours=1)
    idx, _ = find_peaks(q.values, height=threshold, distance=max(1, int(min_sep_h / step_h)))
    return pd.DataFrame({"peak_time": q.index[idx], "peak_cfs": q.values[idx]})

def last_trough(qh: pd.Series, peak_time, max_lookback = MAX_LOOKBACK_H, order_h = TROUGH_ORDER_H): 
    """define the lower bound for the rise seach, the last local minumum before the peak"""
    pre = qh.loc[peak_time - pd.Timedelta(hours = max_lookback):peak_time].interpolate(method = "time", limit_area = "inside").dropna() # define the window before the peak, pre holds the flows from the start of the lookback window up to and including the peak, currently lookback window is 120 hrs
    mins = argrelextrema(pre.values, np.less_equal, order = order_h)[0]   # calculated the relative extrema of the troughs, with the comparator being less than or equal to (not np.less which would fail on a flat line), then keep only the last minima (where the rise starts)
    mins = mins[pre.index[mins] < peak_time] # gaurd to not count the peak itself (convert to timestamp, then boolean mask)
    if len(mins):           
        return pre.index[mins[-1]], False
    return pre.index[0], True           # this means we hit the lookback limit 

def rise_start(qh: pd.Series, peak_time, peak_mag: float):
    """implementing RREDI's rule: from the minimum, the first hour whose rise exceeds the rate threshold for two consecutive hours, or the slope once
    thhen the start is the hour before, and falls back to the minimum if failed"""
    lo, at_limit = last_trough(qh, peak_time)
    seg = qh.loc[lo:peak_time].dropna()
    seg = seg.loc[seg.idxmin():]
    s, base = seg.diff(), seg.iloc[0]
    abs_thr, rel_thr = ABS_SLOPE_FRAC * peak_mag, REL_SLOPE * (peak_mag - base)
    for k in range(1, len(seg)):
        two = k + 1 < len(seg) and s.iloc[k] > abs_thr and s.iloc[k + 1] > abs_thr
        if two or s.iloc[k] > rel_thr:
            return seg.index[k - 1], float(seg.iloc[k - 1]), at_limit
    return seg.index[0], float(base), at_limit


