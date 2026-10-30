"""
Helper functions for sites_01.py. 
  
read_sites(df)       -> standard events table (timestamps) 
get_basins(events)   -> Geodataframe, one row per site with geometry, area, souce, geom_hash
search_windows()      -> one row per event: basin + window_start, window_end, window_src
build_windows()        -> calls the previous three functions in sequence



"""

from __future__ import annotations 
import hashlib, json, urllib.request
from pathlib import Path 
from typing import Callable
import numpy as np 
import pandas as pd 
import geopandas as gpd 
from shapely import set_precision
from shapely.geometry import Point

EQUAL_AREA = 5070
MTBS_COLS = ["Event_ID", "Incid_Name", "Incid_Type", "Ig_Date", "BurnBndAc"]
CACHE = Path(__file__).resolve().parents[1] / "cache" / "streamstats"

# guarding functions to help prevent date and id issues 

def normalize_gage_id(g) -> str | None:
    """takes gage ids and converts to strings, removes decimals, etc. pandas will read in gage id's as float"""
    if g is None or (isinstance(g, float) and np.isnan(g)):  #dont pass nans 
        return None
    s = str(g).strip().split(".")[0] # converts to string and splits, then indexes to the first item (removing decimals)
    return s.zfill(8) if len(s) < 8 else s


def load_basins(path, sites: pd.DataFrame, area_tolerance: float = 0.10) -> gpd.GeoDataFrame:   
    """ensure every site in sites.csv has a basin polygon, return flag if calc'd area is > area_tolerance * gagesii area in either direction"""
    g = gpd.read_file(path)
    g["gage_id"] = g["GAGE_ID"].map(normalize_gage_id)
    g = g[["gage_id", "geometry"]].merge(sites, on = "gage_id", how = "inner", validate="one_to_one")
    missing = set(sites.gage_id) - set(g.gage_id) 
    if missing: 
        raise ValueError(f"{len(missing)} sites have no polygon: {sorted(missing)[:10]}")
    g["area_km2"] = g.to_crs(EQUAL_AREA).area / 1e6 
    g["area_ratio"] = g["area_km2"] / g["drain_sqkm_gii"]
    g["area_flag"] = (g["area_ratio"] - 1).abs() > area_tolerance
    g["geom_hash"] = g.to_crs(4326).geometry.map(geom_hash)
    return g 

def fire_pairs(basins: gpd.GeoDataFrame, mtbs_path, wildfire_only: bool = True) -> pd.DataFrame: 
    """Every MTBS fire intersecting every basin, with burned share of the basin polygon"""
    import shapely
    b = basins.to_crs(EQUAL_AREA)[["gage_id", "area_km2", "geometry"]]
    mask = basins.to_crs(4269).union_all()               # read only parameters near a basin (the file is big!), MTBS is is NAD83, so mask in 4269
    f = gpd.read_file(mtbs_path, columns=MTBS_COLS, mask = mask).to_crs(EQUAL_AREA)
    if wildfire_only: 
        f = f[f["Incid_Type"].eq("Wildfire")]
    f["geometry"] = f.geometry.make_valid()
    j = gpd.sjoin(b, f, predicate="intersects", how = "inner")
    burned = shapely.intersection(j.geometry.values, f.geometry.loc[j["index_right"]].values) 
    return pd.DataFrame({
        "gage_id": j["gage_id"].values, 
        "fire_id": j["Event_ID"].values,
        "incid_name": j["Incid_Name"].values, 
        "ig_date": pd.to_datetime(j["Ig_Date"].values).tz_localize("UTC"),
        "burned_km2": shapely.area(burned) / 1e6,
        "burned_pct": 100 * shapely.area(burned) / 1e6 / j["area_km2"].values,
    })

def search_periods(fires: pd.DataFrame, min_burned_pct: float, years: float) -> pd.DataFrame: 
    """Mark which fires get a post-fire search, and flag overlapping post-fire periods"""
    f = fires.sort_values(by = ["gage_id", "ig_date"]).reset_index(drop = True)
    f["searched"] = f["burned_pct"] >= min_burned_pct  # create a boolean col
    f["search_start"] = f["ig_date"]
    f["search_end"] = f["ig_date"] + pd.DateOffset(years = int(years))
    nxt = f[f.searched].groupby("gage_id")["ig_date"].shift(-1) # keep searched (index bool = keep true) then find the next ig date for each gage 
    f["next_fire_overlaps"] = False # create a default fist, so not Nan for unsearched
    f.loc[f.searched, "next_fire_overlaps"] = nxt < f.loc[f.searched, "search_end"]  # use loc when selecting boolean mask, not wanting to drop False
    return f 

def check_against_events(fires: pd.DataFrame, event_table, gage_ids, min_burned_pct: float) -> pd.DataFrame: 
    """every gage fir paire in the existing (haleys) event table should be found, with a similar burned %
    This is just confirming the code replicates previous methods well"""
    ev = pd.read_csv(event_table, usecols = ["GAGE_ID", "Fire_ID", "burnedarea_per"], low_memory= False) # who cares haha
    ev = ev.drop_duplicates(["GAGE_ID", "Fire_ID"]).assign(gage_id = lambda d: d.GAGE_ID.map(normalize_gage_id))
    ev = ev[ev.gage_id.isin(set(gage_ids))]
    m = ev.merge(fires, left_on = ["gage_id", "Fire_ID"], right_on = ["gage_id", "fire_id"], how = "left")
    m["found"] = m["fire_id"].notna()
    m["pct_diff"] = m["burned_pct"] - m["burnedarea_per"]
    extra = fires[fires.searched].merge(ev[["gage_id", "Fire_ID"]], left_on=["gage_id", "fire_id"],
                                        right_on = ["gage_id", "Fire_ID"], how = "left", indicator=True)
    extra = extra[extra["_merge"] == "left_only"].drop(columns=["Fire_ID", "_merge"])  # searched fires not in haley's table
    return m, extra

def _to_utc(col: pd.Series, tz: str) -> pd.Series:
    """take a column of times and return them as time zone aware timestamps in UTC, times with no timezones are assumed to be 
    local time in the zone you pass as tz (should be IANA name)"""
    t = pd.to_datetime(col, errors='coerce', utc=True)
    if t.dt.tz is None: 
        t = t.dt.tz_localize(tz, ambiguous= "NaT", nonexistent="shift_forward")
    return t.dt.tz_convert("UTC") 


# events table 
def read_sites(df: pd.DataFrame, site_col: str = "site_id", gage_col: str | None = "gage_id",
               lat_col: str | None = "lat", lon_col: str | None = "lon",
               peak_col: str | None = "peak_time", start_col: str | None = None,
               end_col: str | None = None, tz: str = "UTC") -> pd.DataFrame:
    """Standardise an events table. Each row needs a gage ID or lat/lon, and a peak time or a
    start+end pair. A peak given as a date (any lit or field observation)  is flagged `peak_is_date` and widened
    to the whole day in search_windows()."""
    out = pd.DataFrame({"site_id": df[site_col].astype(str)})
    out["gage_id"] = df[gage_col].map(normalize_gage_id) if gage_col in df else None
    out["lat"] = df[lat_col].astype(float) if lat_col in df else np.nan
    out["lon"] = df[lon_col].astype(float) if lon_col in df else np.nan
    for name, col in (("peak_time", peak_col), ("start", start_col), ("end", end_col)):
        out[name] = _to_utc(df[col], tz) if col in df else pd.NaT
    raw = df[peak_col].astype(str) if peak_col in df else pd.Series("", index=df.index)
    out["peak_is_date"] = raw.str.fullmatch(r"\d{4}-\d{2}-\d{2}")
    out["event_id"] = out["site_id"] + "_" + out.groupby("site_id").cumcount().astype(str)

    no_loc = out["gage_id"].isna() & out[["lat", "lon"]].isna().any(axis=1)
    no_time = out["peak_time"].isna() & (out["start"].isna() | out["end"].isna())
    if no_loc.any() or no_time.any():
        bad = out.loc[no_loc | no_time, "event_id"].tolist()
        raise ValueError(f"{len(bad)} events lack a location or a usable time: {bad[:10]}")
    return out


# Basins, functions to cache given basins(already have) or pull basin geometries in input table is just a 
# list of sites 

def geom_hash(geom_4326, precision: float = 1e-6) -> str: 
    """stable hash of the polygon. Key every downstream cache on this, so that a re-delineated bsin 
    keeps its id"""
    return hashlib.sha1(set_precision(geom_4326, precision).normalize().wkb).hexdigest()[:12]

def _nldi_basin(gid: str, timeout: int = 120): 
    CACHE.mkdir(parents=True, exist_ok=True)
    key = CACHE /f"nldi_{gid}.geojson"
    if not key.exists(): 
        req = urllib.request.Request(NLDI_BASIN.format(gid = gid), headers= {"Accept": "application/json"})
        key.write_text(urllib.request.urlopen(req, timeout=timeout).read().decode())
    feats = json.loads(key.read_text()).get("features") or []
    return shape(feats[0]["geometry"]) if feats else None


# function to ffetch the basins with caching so a re-run resumes and doesnt override 

# def _fetch(region: str, lat: float, lon: float, timeout: int=300, 
 #   CACHE.mkdir(parents=True, exist_ok= True)
 #   key = CACHE / f"{region}_{lat:.6f}_{lon:.6f}.json"
 #   if key.exists(): 
  #      return json.loads(key.read_text())
 #   url = f"{BASE}/{region}?lat={lat}&lon={lon}"
 #   last = None
 #   for attempt in range(retries): 
 #       try: 
  #          req = urllib.request.Request(url, headers=)


def _end(geom, last: bool) -> dict: 
    """Get the end and start point of a NHD line, helper function for getting time of conc"""
    if geom.geom_type == "MultiLineString": 
        geom = geom.geoms[-1] if last else geom.goems[0]
    return Point(geom.coords[-1] if last else geom.coords[0])

def flowlines(basin_5070) -> gpd.GeoDataFrame: 
    """Network + non-network NHD flowlines clipped to the basin, with attributes"""
    bb = tuple(gpd.GeoSeries([basin_5070], crs = EQUAL_AREA).to_crs(4326).total_bounds)
    parts = []
    for layer, is_net in LAYERS.items(): 
        try: 
            g = pynhd.WaterData(layer).bybox(bb, box_crs=4326)
        except ZeroMatchedError:                                     # this is no flowlines in the bounding box
            continue
        g.columns = g.columns.str.lower()
        parts.append(g.assign(network = is_net))

    if not parts:
        return gpd.GeoDataFrame(columns=["network", "len_km", "geometry"], geometry="geometry", crs = EQUAL_AREA)
    g = gpd.GeoDataFrame(pd.concat(parts, ignore_index= True), crs =parts[0].crs).to_crs(EQUAL_AREA)
    g["geometry"] = g.geometry.intersection(basin_5070)
    g = g[~geometry.is_empty].copy()
    g["len_km"] = g.length / 1000
    return g 


def sinuosity(main: gpd.GeoDataFrame) -> float: 
    """Mainstem channel length over the straight distance between head and outlet"""
    if main is None or main.empty: 
        return float("nan")
    m = main.sort_values(by="hydroseq") # hydroseq increases upstream 
    d = _end(m.geometry.iloc[0], last = True).distance(_end(m.geometry.iloc[-1], last = False)) 
    return float(m.length.sum() / d ) if d > 0 else float("nan")