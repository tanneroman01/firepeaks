"""Pull all gages with IV records in western US, get all fires in Western US (using NIFC where
mtbs is not up to date, and find all gages where fires have burned their catchement), then keep the
GAGES-II reference basins

steps:
    sites   NWIS sites with IV discharge in the 11 western states    -> nwis_sites.csv
    basins  NLDI basin for each site, one cached geojson per site      -> {site}.geojson, nldi_status.csv
    merge   merge the cached basins                                    -> nldi_basins_merged.gpkg
    fires   NIFC WFIGS perimeters discovered since SINCE (MTBS is local) -> nifc_perims.gpkg
    pairs   each fire on its own vs each basin, keep wildfires burning >= MIN_PCT of the basin
                                        -> fire_basin_pairs_5pct.csv, basins_5pct.gpkg, wildfires_5pct.gpkg
    ref     keep GAGES-II reference (CLASS == Ref) basins, then those < MAX_KM2
                                        -> ref_*_5pct.csv/.gpkg, ref_*_5pct_lt1000km2.csv/.gpkg

"""

######## Warning!!! because we need to be polite to the USGS APIs with req to avoid 404s,
# the basins step takes multiple hours to run, do not run unless absolutely necessary
# (sites already in OUT are skipped, so a rerun only fetches what is missing)
# NIFC perimeters are live data, rerunning fires + pairs will not give the NIFC rows exactly

from __future__ import annotations
import io, json, sys, time, zipfile
from datetime import date
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import requests
import shapely
from firepeaks.inputs.common import EQUAL_AREA, exists 

# ---------------------------------- CONFIG ----------------------------------
NLDI_URL = "https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-{}/basin"
NIFC = ("https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Interagency_Perimeters/FeatureServer/0")
# ---------------------------------- STEPS ----------------------------------

def step_sites(out: Path, states_file: Path, states: list[str], overwrite: bool = False) -> pd.DataFrame | None: 
    """every NWIS site with instantaneous discharge (00060) in the names states"""
    if exists(out / "nwis_sites.csv", overwrite = overwrite): 
        return None
    from pygeohydro import NWIS
    st = gpd.read_file(states_file)
    st = st[st["shapeName"].isin(states)]
    codes = st["shapeISO"].str.replace("US-", "").tolist()   # NWIS wants the bare 2-letter code
    queries = [{"stateCd": c, "parameterCd": "00060", "hasDataTypeCd": "iv"} for c in codes]
    sites = NWIS().get_info(queries).drop_duplicates("site_no")
    print(len(sites), "sites;", sites["huc_cd"].str[:2].value_counts().to_dict())   # 4,359 on 2026-10-29
    sites.drop(columns="geometry", errors="ignore").to_csv(out / "nwis_sites.csv", index=False)
    return sites


def fetch_basin(site: str, out: Path, pause: float = 4.6, max_wait: int = 900) -> str:
    """NLDI basin for one site, cached as {site}.geojson. NLDI rate limits hard, so pause between
    requests and back off when told to"""
    f = out / f"{site}.geojson"
    if f.exists():
        return "cached"
    wait = 60
    while True:
        r = requests.get(NLDI_URL.format(site), params={"simplified": "true"}, timeout=120)

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


def step_basins(out: Path, pause: float = 4.6) -> pd.Series:
    """fetch (or reuse) every site's basin, already cached sites arent re-fetched"""
    ids = pd.read_csv(out / "nwis_sites.csv", dtype={"site_no": str})["site_no"].tolist()
    status = pd.Series({s: fetch_basin(s, out) for s in ids}, name="status")
    status.rename_axis("site_no").to_csv(out / "nldi_status.csv")
    print(status.value_counts().to_dict())     # 4,051 basins and 308 HTTP 404 on 2026-10-29
    return status


def step_merge(out: Path, overwrite: bool = False) -> gpd.GeoDataFrame | None:
    """one layer with every cached basin"""
    if exists(out / "nldi_basins_merged.gpkg", overwrite = overwrite): 
        return None
    parts = []
    for f in out.glob("*.geojson"):
        g = gpd.GeoDataFrame.from_features(json.loads(f.read_text())["features"], crs=4269)
        g["site_no"] = f.stem
        parts.append(g[["site_no", "geometry"]])
    basins = pd.concat(parts, ignore_index=True).to_crs(4326)
    basins["area_sqkm"] = basins.to_crs(EQUAL_AREA).area / 1e6
    basins.to_file(out / "nldi_basins_merged.gpkg", layer="basins", driver="GPKG")
    print(len(basins), "basins")
    return basins


def _nifc_date_fields() -> list[str]:
    """names of the date fields in the NIFC layer. arcgis returns them as epoch milliseconds"""
    meta = requests.get(NIFC, params={"f": "json"}, timeout=60).json()
    return [f["name"] for f in meta["fields"] if f["type"] == "esriFieldTypeDate"]


def step_fires(out: Path, since: date, overwrite: bool = False) -> gpd.GeoDataFrame:
    """every NIFC perimeter discovered since `since` inside the basins' bounding box, all fields"""
    if exists(out / "nifc_perims.gpkg", overwrite = overwrite): 
        return None
    basins = gpd.read_file(out / "nldi_basins_merged.gpkg")
    xmin, ymin, xmax, ymax = basins.total_bounds
    params = {
        "where": f"attr_FireDiscoveryDateTime >= DATE '{since}'",
        "geometry": f"{xmin},{ymin},{xmax},{ymax}",
        "geometryType": "esriGeometryEnvelope",
        "inSR": 4326,
        "spatialRel": "esriSpatialRelIntersects",
        "outFields": "*",
        "orderByFields": "OBJECTID",
        "returnGeometry": "true",
        "outSR": 4326,
        "f": "geojson",
    }
    # the service returns max 2000 records per request, so page through them
    pages, offset = [], 0
    while True:
        r = requests.post(NIFC + "/query", data={**params, "resultOffset": offset}, timeout=120)
        r.raise_for_status()
        js = r.json()
        if "error" in js:   # arcgis can send a 200 with an error in the body
            raise RuntimeError(js["error"])
        if not js["features"]:
            break
        pages.append(gpd.GeoDataFrame.from_features(js["features"], crs=4326))
        offset += len(js["features"])
    perims = pd.concat(pages, ignore_index=True)
    for c in _nifc_date_fields():
        if c in perims:
            perims[c] = pd.to_datetime(perims[c], unit="ms")
    perims.to_file(out / "nifc_perims.gpkg", layer="perims", driver="GPKG")
    print(len(perims), "perimeters in bounding box")   # 7,298 on 2026-10-30
    return perims


def step_pairs(out: Path, mtbs_file: Path, min_pct: float, wildfire_types: list[str], overwrite: bool = False) -> pd.DataFrame | None:
    """test each fire on its own against each basin, keep every wildfire that burned at least
    MIN_PCT of the basin. writes the pairs and the basins and fires that appear in them"""
    files = [out / "fire_basin_pairs_5pct.csv", out / "basins_5pct.gpkg", out / "wildfires_5pct.gpkg"]
    if exists(*files, overwrite = overwrite): 
        return None
    basins = gpd.read_file(out / "nldi_basins_merged.gpkg")
    mtbs = gpd.read_file(mtbs_file)
    nifc = gpd.read_file(out / "nifc_perims.gpkg")

    # make the two sources "mergable", every original column is kept for the fires file
    # (MTBS is not clipped to the western states: clipping would cut perimeters at state lines)
    mtbs = mtbs.assign(source="MTBS", fire_id=mtbs["Event_ID"], fire_name=mtbs["Incid_Name"],
                       fire_date=pd.to_datetime(mtbs["Ig_Date"]), fire_type=mtbs["Incid_Type"])
    nifc = nifc.assign(source="NIFC", fire_id=nifc["OBJECTID"].astype(int).astype(str), fire_name=nifc["poly_IncidentName"],
                       fire_date=nifc["attr_FireDiscoveryDateTime"], fire_type=nifc["attr_IncidentTypeCategory"])
    all_fires = pd.concat([mtbs.to_crs(EQUAL_AREA), nifc.to_crs(EQUAL_AREA)], ignore_index=True)
    all_fires = all_fires[all_fires["fire_type"].isin(wildfire_types)].reset_index(drop=True)
    all_fires["geometry"] = all_fires.geometry.make_valid()
    print(all_fires["source"].value_counts().to_dict(), "wildfires")

    # some basins have invalid geometry, which breaks the intersection
    b = basins[["site_no", "geometry"]].copy()
    b["geometry"] = b.geometry.make_valid()
    b = b.to_crs(EQUAL_AREA)
    b["geometry"] = b.geometry.make_valid()
    b["area_sqkm"] = b.area / 1e6

    # every fire-basin pair that touches
    b_idx, f_idx = all_fires.sindex.query(b.geometry, predicate="intersects")
    pieces = shapely.intersection(b.geometry.values[b_idx], all_fires.geometry.values[f_idx])
    cols = ["source", "fire_id", "fire_name", "fire_date", "fire_type"]
    pairs = pd.concat([b.iloc[b_idx][["site_no", "area_sqkm"]].reset_index(drop=True),
                       all_fires.iloc[f_idx][cols].reset_index(drop=True)], axis=1)
    pairs["burned_sqkm"] = shapely.area(pieces) / 1e6
    pairs["pct_burned"] = 100 * pairs["burned_sqkm"] / pairs["area_sqkm"]

    # keep pairs where this single fire burned at least MIN_PCT of this basin
    pairs = (pairs[pairs["pct_burned"] >= MIN_PCT]
             .sort_values(["site_no", "fire_date"], kind="stable").reset_index(drop=True))
    pairs.to_csv(out / "fire_basin_pairs_5pct.csv", index=False)

    per_basin = pairs.groupby("site_no").agg(n_fires=("fire_id", "size"), max_pct_burned=("pct_burned", "max"),
                                             first_fire=("fire_date", "min"), last_fire=("fire_date", "max"))
    b5 = (b.merge(per_basin, left_on="site_no", right_index=True)
          [["site_no", "area_sqkm", "n_fires", "max_pct_burned", "first_fire", "last_fire", "geometry"]]
          .sort_values("site_no").to_crs(4326))
    b5.to_file(out / "basins_5pct.gpkg", layer="basins", driver="GPKG")

    per_fire = pairs.groupby(["source", "fire_id"]).agg(n_basins=("site_no", "size"), max_pct_burned=("pct_burned", "max"))
    f5 = (all_fires.merge(per_fire, left_on=["source", "fire_id"], right_index=True)
          .sort_values("fire_date", kind="stable").to_crs(4326))
    f5.to_file(out / "wildfires_5pct.gpkg", layer="fires", driver="GPKG")

    print(f"{len(pairs)} fire-basin pairs; {len(f5)} fires; {len(b5)} basins")   
    return pairs


def run(cfg: dict, P: dict) -> None: 
    c = cfg["inputs"]["candidates"]
    out, ow = P["nldi_out"], c.get("overwrite", False)
    out.mkdir(parents = True, exist_ok = True)
    steps = {
        "sites": lambda: step_sites(out, P["states"], c["states"], ow),
        "basins": lambda: step_basins(out, c["nldi_pause_s"]),
        "merge": lambda: step_merge(out, ow),
        "fires": lambda: step_fires(out, date.fromisoformat(str(c["nifc_since"])), ow),
        "pairs": lambda: step_pairs(out, P["mtbs"], c["min_pct"], c["wildfire_types"], ow),

    for s in c["steps"]: 
        print()
