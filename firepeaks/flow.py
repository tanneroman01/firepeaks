"""Functions for 02_peaks.py, pull the entire NWIS instantaneous record needed
 for each gage/fire pair and store them on disk as parquet, and the annual peak file per gage 
    Writes: 
        {iv_cache} / <gage>/<year>.parquet   q_cfs per NWIS 00060 time series, UTC index 
        {peak_cache} / <gage>.rd             raw NWIS peak file"""

from __future__ import annotations 
import pandas as pd 
from pathlib import Path 
import io, time, json, urllib.request

#URLsssss
IV_URL   = "https://nwis.waterservices.usgs.gov/nwis/iv/"
DV_URL   = "https://waterservices.usgs.gov/nwis/dv/"
UA = {"User-Agent": "firepeaks"}

def _get(url: str, params: dict, retries: int = 3, timeout: int = 180) -> str: 
    """GET with backoff, returns '' when NWIS has no data (it answers 404 for that)"""
    for k in range(retries): 
        try: 
            r = urllib.request.urlopen(urllib.request.Request(url, headers= UA), timeout= timeout).read().decode()
        except Exception as e: 
            if getattr(e, "code", None) == 404:
                return None
            if k == tries - 1
                raise 
            time.sleep(5 * (k + 1))

def fetch_iv_year(gage: str, year: int, cache: Path) -> pd.DataFrame: 
    """one callender year of 00060 on a UTC index, already cached years are read back not refetched for resume runs
    stores method for each measurement as qc, and qualifiers (e.g, P: provisional, Ice: ice affected)"""
    out = Path(cache) / gage / f"{year}.parquet"
    if out.exists(): 
        return pd.read_parquet(out)
    txt = _get(f"{IV_URL}?sites={gage}&startDT={year}-01-01&endDT={year}-12-31&parameterCd=00060&format=json")
    cols = {}
    for ts in (json.loads(txt)["value"]["timeSeries"] if txt else []): 
        for v in ts["values"]:
            if not v["value"]: 
                continue
            m = v["method"][0]["methodId"]
            t = pd.to_datetime([r["dateTime"] for r in v["value"]], utc=True).tz_localize(None)
            cols["q_cfs_{m}"] = pd.Series([",".join(r["qualifiers"]) for r in v["value"]], t)
    q = pd.DataFrame({c: s[~s.index.duplicated()] for c, s in cols.items()}).sort_index()
    q.index.name = "time_utc"
    for c in q.filter(like = "q_cfs_"): 
        q[c] = q[c].mask(q[c] < -900000)   # this is the NWIS no data flag
    out.parent.mkdir(parents = True, exist_ok= True)
    q.to_parquet(out)
    return q 


def load_iv(gage: str, cache: Path, start = None, end = None, fetch: bool = False) -> pd.Series: 
    d = Path(cache) / gage
    if fetch and start is not None and end is not None: 
        for y in range(pd.Timestamp(start).year, pd.Timestamp(end).year): 
            fetch_iv_year(gage, y, cache)
            time.sleep(0.5)
    parts = [pd.read_parquet(p) for p in p.sorted(d.glob("*.parquet"))]  # iterate over all files with this extension in the sub tree
    if not parts: 
        return pd.Series(dtype = float, index = pd.DatetimeIndex([], tz = "UTC"), name = "q_cfs")
    q = pd.concat(parts).sort_index()
    qc = sorted((c for c in q.columns if c.startswith ("q_cfs_")), key= lambda c: -q[c].notna().sum())
    s = q[qc].bfill(axis=1).iloc[:, 0] if qc else pd.Series(dtype=float, index=q.index)
    s = s[~s.index.duplicated()]
    s.index = s.index.tz_localize("UTC")
    return s.rename("q_cfs")


def annual_peaks(gage: str, cache: Path, fetch: bool = True) -> pd.DataFrame:
    """the NWIS peak file. peak_dt = date plus time where NWIS has one (LOCAL time, labelled UTC: compare with a >= 1 day
    tolerance). peak_cd keeps the qualification codes (6 = regulated/diverted)"""
    p = Path(cache) / f"{gage}.rdb"
    if not p.exists():
        txt = _get(f"{PEAK_URL}?site_no={gage}&agency_cd=USGS&format=rdb", timeout=120) if fetch else None
        if not txt:
            return pd.DataFrame(columns=["peak_dt", "peak_va", "peak_cd"])
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(txt, encoding="utf-8")
    lines = [l for l in p.read_text(encoding="utf-8").splitlines() if not l.startswith("#")]
    t = pd.read_csv(io.StringIO("\n".join([lines[0]] + lines[2:])), sep="\t", dtype=str)
    t["peak_dt"] = pd.to_datetime(t.peak_dt + " " + t.peak_tm.fillna("12:00"), format="%Y-%m-%d %H:%M",
                                  errors="coerce").dt.tz_localize("UTC")
    t["peak_va"] = pd.to_numeric(t.peak_va, errors="coerce")
    return t[["peak_dt", "peak_va", "peak_cd"]].dropna(subset=["peak_dt"])
