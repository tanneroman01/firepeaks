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



















def fire_perimeters_in(esri_polygon, since = SINCE):
    """Return a GeoJSON of fire perimeters that intersect the delineated basin"""
    params = {
        "geometry": esri_polygon,
        "geometryType": "esriGeometryPolygon",
        "inSR": 4326,
        "spatialRel":"esriSpatialRelIntersects",
        "outFields": "poly_IncidentName,poly_GISAcres,attr_FireDiscoveryDateTime,poly_PolygonDateTime",
        "returnGeometry":"true",
        "outSR":4326,
        "f": "geojson",
    }
    if since:
        params["where"] = f"attr_FireDiscoveryDateTime >= DATE '{since}'"

    r = requests.post(NIFC, data = params, timeout= 60)   # POST: a basin polygon is too long for a URL
    r.raise_for_status()
    return r.json()   # may hold an "error" key instead of features; the caller checks


def fire_table(fires, basin_geojson):
    """One row per fire: name, discovery date, acres, and the share of the basin inside the perimeter."""
    f = gpd.GeoDataFrame.from_features(fires["features"], crs = 4326).to_crs(5070)
    b = gpd.GeoDataFrame.from_features(basin_geojson["features"], crs = 4326).to_crs(5070).geometry.union_all()
    f["burned_pct"] = f.geometry.intersection(b).area / b.area * 100
    f["discovered"] = [datetime.fromtimestamp(ms / 1000, tz = timezone.utc).date() if ms else None
                       for ms in f["attr_FireDiscoveryDateTime"]]
    f["label"] = [f"{n} ({d}) - {p:.0f}% of basin burned"
                  for n, d, p in zip(f["poly_IncidentName"], f["discovered"], f["burned_pct"])]
    return f.sort_values("burned_pct", ascending = False).reset_index(drop = True)

