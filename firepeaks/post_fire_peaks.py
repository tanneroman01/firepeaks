import warnings
import numpy as np 
import pandas as pd 
from scipy.signal import find_peaks, argrelextrema, peak_prominences
from pre_fire_peaks import (MIN_SEP_H, TROUGH_RATIO, independent_peaks, melt_season, utc) 


# some config parameters (tunable)

REL_SLOPE = 0.05
ABS_SLOPE_FRAC = 0.01   # replaces RREDIS fixed 1 cfs/h, using a fraction of the peak instead (allows for tiny basin response flows)
TROUGH_ORDER_H = 6      # a trough must be the minimum within +/- this many hrs
MAX_LOOKBACK_H = 120

# event filters 
SNOW_PCT_MIN = 15.0 
AC24_MIN = 0.3 
DIURNAL_DAYS = 5
CREST_WIN_H = 48 
PROM_MIN = 0.5 
PROM_WIN_H = 72 # window for the prominence calc 

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
    if len(q) < 3 or not np.isfinite(threshold):
        return pd.DataFrame({"peak_time": pd.DatetimeIndex([], tz = "UTC"), "peak_cfs": pd.Series(dtype = float)})
        pk = independent_peaks(q, pd.Timedelta(hours = min_sep_h), trough_ratio).sort_index()
        pk = pk[pk > threshold]
        return pd.DataFrame({"peak_time": pk.index, "peak_cfs": pk.values})


#--------------event filters run after find_events-----------------

def ac24(qh: pd.Series, pk, days: int = DIURNAL_DAYS) -> float: 
    """lag 24 h auto correlation of hourly flow minum its 24 h centred mean over 'days' since the event"""
    h = qh.loc[pk - pd.Timedelta(days= days):pk -
               pd.Timedelta(hours = 1)].interpolate(limit=3)
    if h.notna().sum() < 72: 
        return np.nan 
    r = (h - h.rolling(24, center=True, min_periods= 18).mean()).dropna()
    return float(r.autocorr(24)) if len(r) > 48 else np.nan 


# quality check on the flood peak, give it peak time pk, and the peaks flow qpk. Look at the hourly flow series qh around the time and answer: 
#   dom(dominance) "is it really the top of the flood?" if nearby is higher, peak detected isnt the highest point of event
#   prom(prominence) "does it stand out in nearby peaks?" measure drop on either side as frac of peak before reaching higher flow 

def crest(qh: pd.Series, pk, qpk: float, win_h: float = CREST_WIN_H, prom_win_h: float = PROM_WIN_H):
    """dom = hourly max flow within +- win_h / peak (>1 = a larger crest nearby), 
    prom = topographic prominence of the peak within +- prom_win_h window/peak"""
    win = qh.loc[pk - pd.Timedelta(hours= win_h): pk + pd.Timedelta(hours=win_h)] 
    dom = float(win.max() / qpk) if len(win) else np.nan 
    w = qh.loc[pk - pd.Timedelta(hours= prom_win_h):pk + pd.Timedelta(hours = prom_win_h)].interpolate(method = "linear", limit= 6).dropna() # fill dont drop so that readings either side of gap dont become neighbors 
    if len(w) < 3: 
        return dom, np.nan  # crest needs atleast 3 points (1 on either side, so nan for prom)
    i = int(np.argmin(np.abs(w.index-pk).total_seconds()))
    i = max(range(max(0, i-1), min(len(w), i+2)), key= lambda k: w.iloc[k]) # the hourly crest at the peak
    if not 0 < i < len(w) - 1:   
        return dom, np.nan   # give up if crest is at the edge of the window (if not reversest the chained comparison)
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore")
        prom = peak_prominences(w.to_numpy(), [i])[0][0]   # scipy gives a bunch of warnings for this func, e.g. flat top is prom of 0, but still meaningful in our case 
    return dom, float(prom / qpk)


class EventFilter: 
    """per - gage-fire filter, built once from pre-fire record, gives melt season, basins snow share, and reason: 
    keep, snowmelt, or minor crest"""
    def __init__(self, qpre: pd.Series, snow_pct: float, snow_pct_min: float = SNOW_PCT_MIN, ac24min: float = AC24_MIN, 
                 prom_min: float = PROM_MIN): 
        self.snow_pct, self.ac24_min, self.prom_min = snow_pct, ac24_min, prom_min
        self.snowy = bool(np.isfinite(snow_pct) and snow_pct >= snow_pct_min)
        self.season = melt_season(qpre) if self.snowy else None 

    def features(self, qh: pd.Series, pk, qpk: float) -> dict: 
        in_melt = bool(self.snowy and pk.dayofyear in self.season["melt_doys"])
        a = ac24(qh, pk) if in_melt else np.nan 
        dom, prom = crest(qh, pk, qpk)
        return dict(in_melt = in_melt, ac24 = a, dom48 = dom, prom_ratio = prom)

    def reason(self, qh: pd.Series, pk, qpk: float) -> str: 
        f = self.features(qh, pk, qpk)
        if f["in_melt"]  and f["ac24"] >= self.ac24_min: 
            return "snowmelt"
        if f["dom48"] > 1.0 and f["prom_ratio"] < self.prom_min: 
            return "minor_crest"
        return ""


def detect(q15: pd.Series, threshold: float, start, end, ef: EventFilter | None = None, 
           qh: pd.Series | None = None, min_sep_h = MIN_SEP_H, trough_ration: float = TROUGH_RATIO) -> pd.DataFrame: 
    """find events with the filters, every peak above threshold with 'excluded' = "" (counted), "snowmelt" or "minor_crest" Excluded peaks 
    are kept in the table not deletedm qh is the hourly max of same record, pass it in when calling repeatedly on one gage else None """
    ev = find_events(q15, threshold, start, end, min_sep_h)
    if ef is None or ev.empty: 
        return ev.assign(excluded = "")
    qh = hourly_max(q15) if qh is None else qh
    return ev.assign(excluded = [ef.reason(qh, t, v) for t, c in zip(ev.peak_time, ev.peak_cfs)])


#---------------------------- When is the rise start?-------------------------------------------

def last_trough(qh: pd.Series, peak_time, max_lookback: int = MAX_LOOKBACK_H, order_h = TROUGH_ORDER_H): 
    """lower bound for the rise search, the last local min before the peak"""