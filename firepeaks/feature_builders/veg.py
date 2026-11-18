"""vegetation recovery features for every event in the table, build in one call(ie. python veg.py)
        ratio_LAI_near MODIS MCD15A3H LAI over the event window / mean of the same calendar window in the 
         5 years pre-fire
        
        ratio_NDVI  landsat C2 L2, same windows and baselines as above 
         
        delta_NDVI_near event_NDVI - baseline NDVI 
         
        n_base_lai, n_base_ndvi   baselines with data (0-5) MODIS starts 2002-07, so pre-2003 fires get 0
         these are explicit availability indicators so that drees do not learn NaN as era proxy
          
           
    windows (indentical for both sensors): 
     event    W = [max(e-90d, g + 15d), e]  where e is the event data and g is the ignition date 
     baseline W - k yr for k - k0 .. k0 + 4, k0 the smallest k >= 1 with e - k yr <= g (all end pre-fire)
     ratio    flood/ mean (available baselines); Nan if flood missing, no baseline, or mean <0.1 
     
Landsat gap filling: a pixel with no clear scene in the 90 day window falls back to its median over [max(e-180d, g + 15d), e] and baselines are 
widened the same way 

Region is basin burn area overlap, LAI is reduced at 500m, NDVI at 30m"""

# can we use the hi-q lai collection? 
# apparently terra and aqua are drifting off their orbits and scheduled to shut down Fep 2027, something to keep in mind
from __future__ import annotations 
import hashlib 
import time 
import pandas as pd 
import numpy as np 
import geopandas as gpd 
import ee
from pathlib import Path 
import json  
from concurrent.futures import ThreadPoolExecutor


NEAR_DAYS, WIDE_DAYS, FLOOR_DAYS, N_BASE = 90, 180, 15, 5 
MIN_BASE = 0.1 # to prevent ratio's from blowing up 
LAI_COL = "MODIS/061/MCD15A3H"
# FparLai_QC mask to get rid of clouds (bits 3-4 clear)
LAI_QC1, LAI_QC2 = 3, 4
LANDSAT = {
      "LANDSAT/LT05/C02/T1_L2": ("SR_B4", "SR_B3"),
    "LANDSAT/LE07/C02/T1_L2": ("SR_B4", "SR_B3"),
    "LANDSAT/LC08/C02/T1_L2": ("SR_B5", "SR_B4"),
    "LANDSAT/LC09/C02/T1_L2": ("SR_B5", "SR_B4"),
}

# Landsat QA
QA_BAD = 0b10111111  # QA_PIXEL fill, dilated cloud, cirrus, cloud, shadow, snow, water 
SCALE = {"lai": 500, "ndvi": 30}
CHUNK = {"lai": 80, "ndvi": 20}  # bands per reduceRegion, landsat is heavier imagery so gets a smaller batch so that a single request doesnt run out of memory
COLS = ["ratio_LAI_near", "n_base_lai", "ratio_NDVI_near", "delta_NDVI_near", "n_base_ndvi"]

def init(project: str = "dbflow-480621"): 
    ee.Initialize(project = project)


# defining windows (need to add some leap year handling in here? )

def windows(ig: pd.Timestamp, ev: pd.Timestamp): 
    """Event window and its 5 pre-fire baselines windows, each has (start, end_exclusive, wide_start)"""
    floor = ig + pd.Timedelta(days  = FLOOR_DAYS)
    lo = max(ev - pd.Timedelta(d = NEAR_DAYS), floor)
    if lo >= ev:   # the event is within 15 days of ignition, thus no good window 
        return None, []
    wide = max(ev - pd.Timedelta(days = WIDE_DAYS), floor)
    k0 = next(k for k in range(0, 100) if ev - pd.DateOffset(years = k) <= ig)  # find how many years we need to step back from ev to get on or before the ignition date 
    yr = [pd.DateOffset(years = k) fork in range(k0, k0 + N_BASE)]  # take the same window from just before the fire and shift it back into 5 years before the fire 
    return (lo, ev, wide), [(lo- y, ev - y, wide - y) for y in yr]  


# earth engine boiler plat e 
def _d(t): 
    return t.strftime("%Y-%m-%d")

def wkey(sensor, w): 
    """turns each window into a text key built from its dates """
    s, e, sw, = w
    return f"{s:%Y%m%d}_{e:%Y%m%d}" if sensor == "lai" else f"{sw:%Y%m%d}_{s:%Y%m%d}_{e:%Y%m%d}"

def _median(col):
    """median composites, fully masked when the window has no data (e.g. pre-sensor era)""" 
    empty = ee.Image.constant(0).toFloat().rename("v").updateMask(0)
    return col.merge(ee.ImageCollection([empty])).median()

def _lai(geom, w): 
    s, e, _ = w

    def f(img): 
        ok = img.select("FparLai_QC").bitwiseAnd(LAI_QC1, LAI_QC2).eq(0) # select pixels where bits 3 and 4 are 0 
        return img.select("Lai").multiply(0.1).toFloat().updateMask(ok).rename("v")
    return _median(ee.ImageCollection(LAI_COL).filterBounds(geom).filterDate(_d(s), _d(e)).map(f))

def _ndvi_col(geom, s, e): 
    def prep(nir, red): 
        def f(img): 
            sr = img.select([nir, red], ["nir", "red"]).multiply(0.0000275).add(-0.2) # usgs converts to Int 16 so we need to rescale back to true reflectance
            ok = img.select("QA_PIXEL").bitwiseAnd(QA_BAD).eq(0)
            return sr.NormalizedDifference(["nir"], ["red"]).rename("v").toFloat().updateMask(ok)
        return f 
    cols = [ee.ImageCollection(c).filterBounds(geom).filterDate(_d(s), _d(e)).map(prep(*b))
            for c, b in LANDSAT.items()]
    out = cols[0]
    for c in cols[1:]: 
        out = out.merge(c)
    return out 

def _ndvi(geom, w): 
    s, e, sw = w
    near = _median(_ndvi_col(geom, s, e))
    return near.masked(_median(_ndvi_col(geom, s, e))) if sw < s else near

IMG = {"lai": _lai, "ndvi": _ndvi}

def _get_em(obj, tries = 6): 
    """get with backoff for rate errors"""
    for k in range(tries): 
        try:
            return obj.getInfo()
        except Exception as e: 
            if k == tries - 1 or not any(t in str(e) for t in ("concurrent", "Too many", "rate", "timed out")): 
                raise
            time.sleep(5*2 ** k) # exponential backoff, worst case here is 2.5 min

def _reduce(geom, sensor, wins): 
    """ one reduce region per chunk instead of per window basin"""
    out = {}  # empty dict for results
    for i in range(0, len(wins), CHUNK[sensor]): 
        part = wins[i:i + CHUNK[sensor]] # define what windows go in each batch
        names = ["w" + wkey(sensor, w) for w in part] #dont let band names start with a digit (poor mans static typing)
        img = ee.Image.cat([IMG[sensor](geom, w).rename(n) for w, n in zip(part, names)]) # build one image with multi band (one per window)
        res = _get_em(img.reduceRegion(reducer=ee.Reducer.mean(), geometry=geom, scale=SCALE[sensor],
                                    bestEffort=True, maxPixels=1e9, tileScale=4)) # average over basin split into smaller tiles to avoid memory errors
        out.update({n[1:]: res.get(n) for n in names})
    return out

# Getting the burn/basin intersections and caching them 
# perhaps should re-factor such that burn/basin intersections are cached earlier in the run for other feature calc's too? 

def pair_regions(events, basins, perims, simplify_deg = 0.0003): 
    b = basins.assign(GAGE_ID = basins.GAGE_ID.astype(str).str.zfill(8))
    bg = dict(zip(b.GAGE_ID, b.geometry))
    p = perims.to_crs(b.crs)
    pg = dict(zip(p.Event_ID, p.geometry))
    out = {}
    for g, f in events[["GAGE_ID", "Event_ID"]].dropna().drop_duplicates().itertuples(index = False): 
        if g not in bg or f not in pg: 
            print(f"veg: no geometry for {g} / {f}")
            continue
        inter = bg[g].intersection(pg[f])
        if inter.is_empty: 
            continue
        out[(g, f)] = gpd.GeoSeries([inter], crs = b.crs).to_crs(4326).iloc[0].simplify(simplify_deg)

# fetch or read every window each gage-fire pairs events need, return dict of features 
def pair_values(gage, fire, region, ev_pairs, cache_dir: Path):
    """Fetch (or read) every window value this pair's events need; returns {'lai': {...}, 'ndvi': {...}}."""
    h = hashlib.sha1(region.wkb).hexdigest()[:10]
    path = cache_dir / f"{gage}_{fire}_{h}.json"
    cache = json.loads(path.read_text()) if path.exists() else {"lai": {}, "ndvi": {}}
    need = {"lai": set(), "ndvi": set()}
    for ig, ev in ev_pairs:
        flood, base = windows(ig, ev)
        for w in ([flood] if flood else []) + base:
            for s in need:
                if wkey(s, w) not in cache[s]:
                    need[s].add(w)
    if any(need.values()):
        geom = ee.Geometry(region.__geo_interface__)
        for s, wins in need.items():
            if wins:
                cache[s].update(_reduce(geom, s, sorted(wins)))
        path.write_text(json.dumps(cache))
    return cache

# calculate the ratio's and deltas of the lai/ndvi for the basins 
def _ratio(c, flood, base):
    f = c.get(flood)
    b = [v for v in (c.get(k) for k in base) if v is not None]
    if f is None or not b or abs(np.mean(b)) < MIN_BASE:
        return np.nan, np.nan, len(b)
    return f / np.mean(b), f - np.mean(b), len(b)

# global builder function 

def build(events: pd.DataFrame, basins: gpd.GeoDataFrame, perims: gpd.GeoDataFrame,
          cache_dir: str | Path, workers: int = 4) -> pd.DataFrame:
    """Return `events` (same rows, same order) with COLS appended.

    events needs GAGE_ID, Fire_ID, Ig_Date, event_date. basins needs GAGE_ID + geometry
    (gagesii_basins.gpkg); perims needs Event_ID + geometry (MTBS). Call init() first.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    ev = events.copy()
    ev["GAGE_ID"] = ev.GAGE_ID.astype(str).str.zfill(8)
    ev["Ig_Date"], ev["event_date"] = pd.to_datetime(ev.Ig_Date), pd.to_datetime(ev.event_date)

    regions = pair_regions(ev, basins, perims)
    groups = {k: list(zip(g.Ig_Date, g.event_date)) for k, g in ev.groupby(["GAGE_ID", "Fire_ID"])
              if k in regions}
    print(f"veg: {len(ev)} events, {len(groups)} (gage, fire) pairs", flush=True)

    def run(k):
        try:
            return k, pair_values(*k, regions[k], groups[k], cache_dir)
        except Exception as e:  # noqa: BLE001 - report and continue; a rerun fetches only what is missing
            print("veg FAIL", *k, str(e)[:200], flush=True)
            return k, None
    with ThreadPoolExecutor(workers) as ex:
        caches = dict(ex.map(run, groups))

    rows = []
    for r in ev.itertuples():
        o = dict.fromkeys(COLS, np.nan)
        c = caches.get((r.GAGE_ID, r.Fire_ID))
        flood, base = windows(r.Ig_Date, r.event_date) if pd.notna(r.Ig_Date) else (None, [])
        if c is not None and flood is not None:
            o["ratio_LAI_near"], _, o["n_base_lai"] = _ratio(
                c["lai"], wkey("lai", flood), [wkey("lai", w) for w in base])
            o["ratio_NDVI_near"], o["delta_NDVI_near"], o["n_base_ndvi"] = _ratio(
                c["ndvi"], wkey("ndvi", flood), [wkey("ndvi", w) for w in base])
        rows.append(o)
    out = pd.concat([events.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    print("veg coverage:", out[COLS].notna().sum().to_dict(), flush=True)
    return out


# main so that we can call this script directly 
if __name__ == "__main__":
    import yaml
    repo = Path(__file__).resolve().parents[2]
    cfg = yaml.safe_load((repo / "config" / "default.yaml").read_text())["paths"]
    init()
    table = pd.read_csv(repo / cfg["event_table"], dtype={"GAGE_ID": str})
    res = build(table, gpd.read_file(repo / cfg["basins"]), gpd.read_file(repo / cfg["mtbs"]),
                repo / "cache" / "veg")
    res.to_csv(repo / "outputs" / "veg_features.csv", index=False)  # new file; never overwrites the input table

    