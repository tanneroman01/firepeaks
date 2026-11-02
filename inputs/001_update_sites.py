"""Pull all gages with IV records in western US, get all fires in Western US (using NLDI where
mtbs is not up to date, and find all gages where fires have burned their catchement)"""

######## Warning!!! because we need to be polite to the USGS APIs with req to avoid 404s,
# this script takes multiple hours to run, do not run unless absolutely necessary 


import pandas as pd 
import numpy as np 
import geopandas as gpd 
import requests
from pathlib import Path
import yaml
from datetime import date, time, timedelta, timezone
import pynhd
import pygeohydro as gh
from pygeohydro import NWIS 

import os 

pd.set_option('display.max_columns', None)

# NIFC req constants 
NIFC = ("https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Interagency_Perimeters/FeatureServer/0/query")
SINCE = (date(2024, 12, 1)) 

# look at all the available NWIS stations within a bounding box (Western US, can we use the states file?)
# States path 

STATES = r"C:\Users\A02306692\Downloads\PFF\release\Release_repo_in_progress\data\geoBoundaries-USA-ADM1_simplified.geojson"

states = gpd.read_file(STATES)

# select western US states 
WESTERN = ["Washington", "Montana", "Wyoming", "Utah", "Oregon", "Idaho", "Colorado", "Arizona", "New Mexico", "California", "Nevada"]
states = states[states["shapeName"].isin(WESTERN)]
states["dissolve_col"] = 0
sd = states.dissolve(by = "dissolve_col")

# nwis = NWIS()

# Check how many gages in total exist in western US 
codes = states["shapeISO"].str.replace("US-", "").tolist()

queries = [{"stateCd": c, "parameterCd": "00060", "hasDataTypeCd": "iv"} for c in codes]
sites = nwis.get_info(queries)

sites = sites.drop_duplicates("site_no")
print(len(sites))

print(len(sites), sites["site_no"].is_unique)
print(sites["huc_cd"].str[:2].value_counts())

# pull all of the basin geometries 

import time, json, requests
from pathlib import Path

URL = "https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-{}/basin"
CACHE = Path("nldi_basins"); CACHE.mkdir(exist_ok=True)

def fetch_basin(site, pause=4.6, max_wait=900):
    f = CACHE / f"{site}.geojson"
    if f.exists():
        return "cached"
    wait = 60
    while True:
        r = requests.get(URL.format(site), params={"simplified": "true"}, timeout=120)

        if r.status_code == 200:
            f.write_text(r.text)
            remaining = int(r.headers.get("X-Ratelimit-Remaining", 100))
            if remaining < 20:
                print(f"{remaining} requests left in window, sleeping 10 min")
                time.sleep(600)
            else:
                time.sleep(pause)
            return "ok"

        if r.status_code == 429 or (r.status_code in (502, 503) and "RATE_LIMIT" in r.text):
            print(f"{site}: rate limited, sleeping {wait}s")
            time.sleep(wait)
            wait = min(wait * 2, max_wait)
            continue

        time.sleep(pause)
        return f"HTTP {r.status_code}"     # e.g. 404: site not in NLDI, retrying won't help


# take the cached files, merge them, and export to a gpkg
parts = []
for f in CACHE.glob("*.geojson"):
    g = gpd.GeoDataFrame.from_features(json.loads(f.read_text())["features"], crs=4269)
    g["site_no"] = f.stem
    parts.append(g[["site_no", "geometry"]])

basins = pd.concat(parts, ignore_index=True).to_crs(4326)
basins["area_sqkm"] = basins.to_crs(5070).area / 1e6

bad = {s: v for s, v in status.items() if v not in ("ok", "cached")}
print(len(basins), "basins;", len(bad), "failed:", pd.Series(bad).value_counts())

# export file 

basins.to_file(CACHE / "nldi_basins_merged.gpkg", layer = "basins", driver = "GPKG")


STATES = r"C:\Users\A02306692\Downloads\PFF\release\Release_repo_in_progress\data\geoBoundaries-USA-ADM1_simplified.geojson"
states = gpd.read_file(STATES)

WESTERN = ["Washington", "Montana", "Wyoming", "Utah", "Oregon", "Idaho", "Colorado", "Arizona", "New Mexico", "California", "Nevada"]
states = states[states["shapeName"].isin(WESTERN)]
states["dissolve_col"] = 0
sd = states.dissolve(by = "dissolve_col")
FIRES = r"C:\Users\A02306692\Downloads\PFF\altRREDI\inputs\mtbs_perims_DD.shp"
fires = gpd.read_file(FIRES)
sd = sd.to_crs(epsg = 4326)
fires =fires.to_crs(epsg= 4326)
# mask fires by western US states 
fc = fires.clip(mask = sd)


# load in basins with gages 

BASINS = r"C:\Users\A02306692\Downloads\nldi_basins\nldi_basins_merged.gpkg"
basins = gpd.read_file(BASINS)

print(len(basins))


def fire_perimeters_in(basins, since = SINCE): 
    """return a geodataframe of the perims since `since` inside the basins' bounding box"""
    xmin, ymin, xmax, ymax = basins.total_bounds
    params = {
        "where": f"attr_FireDiscoveryDateTime >= DATE '{since}'",
        "geometry": f"{xmin},{ymin},{xmax},{ymax}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": 4326,
        "spatialRel":"esriSpatialRelIntersects",
        "outFields": "OBJECTID,poly_IncidentName,poly_GISAcres,attr_IrwinID,attr_IncidentTypeCategory,attr_POOState,attr_FireDiscoveryDateTime,poly_PolygonDateTime",
        "orderByFields": "OBJECTID",
        "returnGeometry":"true",
        "outSR":4326,
        "f": "geojson",
    }

    # the service returns max 2000 records per request, so page through them
    pages, offset = [], 0
    while True:
        r = requests.post(NIFC, data = {**params, "resultOffset": offset}, timeout=120)
        r.raise_for_status()
        js = r.json()
        if "error" in js:   # arcgis can send a 200 with an error in the body
            raise RuntimeError(js["error"])
        if not js["features"]:
            break
        pages.append(gpd.GeoDataFrame.from_features(js["features"], crs = 4326))
        offset += len(js["features"])

    perims = pd.concat(pages, ignore_index = True)
    perims["attr_FireDiscoveryDateTime"] = pd.to_datetime(perims["attr_FireDiscoveryDateTime"], unit = "ms")
    return perims

perims = fire_perimeters_in(basins)
print(len(perims), "perimeters in bounding box")


# test each fire on its own agaisnt basin, keep every pair where that one fire covers
# at least 5% of the basin

b = basins[["site_no", "geometry"]].to_crs(EQUAL_AREA)
b["geometry"] = b.geometry.make_valid()
b["area_sqkm"] = b.area / 1e6

# every fire-basin pair that touches
b_idx, f_idx = all_fires.sindex.query(b.geometry, predicate = "intersects")
pieces = shapely.intersection(b.geometry.values[b_idx], all_fires.geometry.values[f_idx])

pairs = pd.concat([b.iloc[b_idx][["site_no", "area_sqkm"]].reset_index(drop = True),
                   all_fires.iloc[f_idx].drop(columns = "geometry").reset_index(drop = True)], axis = 1)
pairs["burned_sqkm"] = shapely.area(pieces) / 1e6
pairs["pct_burned"] = 100 * pairs["burned_sqkm"] / pairs["area_sqkm"]

# keep pairs where this single fire burned at least MIN_PCT of this basin
fire_basin_5pct = pairs[pairs["pct_burned"] >= MIN_PCT].sort_values(["site_no", "fire_date"]).reset_index(drop = True)


n_fires = fire_basin_5pct.drop_duplicates(["source", "fire_id"]).shape[0]
print(len(fire_basin_5pct), "fire-basin pairs;", n_fires, "fires;", fire_basin_5pct["site_no"].nunique(), "basins")