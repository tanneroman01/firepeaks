"""Functions for 02_peaks.py, pull the entire NWIS instantaneous record needed
 for each gage/fire pair and store them on disk"""

import pandas as pd 
import numpy as np 
import io, time
import requests

#URLsssss
IV_URL   = "https://nwis.waterservices.usgs.gov/nwis/iv/"
DV_URL   = "https://waterservices.usgs.gov/nwis/dv/"
SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
UA = {"User-Agent": "firepeaks"}

def _get(url: str, params: dict, retries: int = 3, timeout: int = 180) -> str: 
    """GET with backoff, returns '' when NWIS has no data (answers 404 for that)"""
    for k in range(retries): 
        try: 
            r = requests.get(url, params=params, headers=UA, timeout=timeout)
            if r.status_code == 404: 
                return ""
            r.raise_for_status()
            return r.text
        except requests.RequestException: 
            if k == retries - 1: 
                raise
            time.sleep(2 ** k)


def gage_span(fires: pd.DataFrame, gage_id: str, record_start=None )  -> tuple[pd.Timestamp, pd.Timestamp] | None: 
    f = fires[(fires["gage_id"] == gage_id) & fires["searched"]]
    if f.empty: 
        return None
    start = pd.Timestamp(record_start or "1980-01-01", tz = "UTC")
    end = min(f["search_end"].max(), pd.Timestamp.now(tz="UTC"))
    return start, end


#def fetch_iv(gage_id: str, start, end, chuck_years: int = 5) -> pd.Series: 
 #   """NWIS instantenough discharge (parameterCd is 00060 for IV) as a cfs Series with a UTC index.
 #   Requests in chucks to avoid timeouts, convert tz_cd to UTC here"""
 #   parts, t, end = [], pd.Timestamp(start), pd.Timestamp(end) 
 #   while t < end: 
 #       t1 = min(t + pd.DateOffset(years= chuck_years), end)
 #       df = parse_rdb(_get(IV_URL, {"sites": gage_id, "parameterCd": "00060", "format": "rdb", 
 #                                    "startDT": t.strftime("%Y-%m-%d"), "endDT": t1.strftime("%Y-%m-%d")}))
 #       if len(df): 
 #           parts.append(df)
 #       t = t1
 #       time.sleep(0.5)  #tryna be polite lol 
 #   if not parts: 
 #       return pd.DataFrame(columns= ["q_cfs", "qual"])
 #   df = pd.concat(parts, ignore_index=True)
 #   vcols = [c for c in df.columns if c.endswith("_00060")]
 #   vcol = max()
#def featch_dv(gage_id: str, start, end) -> pd.Series: 
#    """daily mean discharge for coverage check of IV series"""

#def clean(q: pd.Series, qualifiers: pd.Series | None = None) -> pd.Series: 
#    """ drop duplicate timestamps, sort, set errors to Nan and keep provisional data but flag it"""
#    q = 

#def load_flow(gage_id: str, start, end, cache_dir: Path, refresh: bool = False) -> pd.Series: 
#    """02_peaks.py will call this, read {cache_dire} / {gage_id}.parquet if it convers [start,end]
#    otherwise fetch_iv, write the cache, and return the slice"""

#def coverage(q: pd.Series, start, end, freq: str = "15min") -> dict: 
#    """QC per window: fractopm of expected timestams present, first/last timestamp, longest gap. 
#    Stops a window with a gap over a flood from reading as no peak"""

