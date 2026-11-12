"""Pull all gages with IV records in western US, get all fires in Western US (using NIFC where
mtbs is not up to date, and find all gages where fires have burned their catchement), then keep the
GAGES-II reference basins

Reproduces the files in OUT (Downloads/nldi_basins) made 2026-10-29/30 in explore.ipynb and
explore_2.ipynb. Each step reads the previous step's files from OUT, so any step can be run alone

steps, in order:
    sites   NWIS sites with IV discharge in the 11 western states    -> nwis_sites.csv
    basins  NLDI basin for each site, one cached geojson per site      -> {site}.geojson, nldi_status.csv
    merge   merge the cached basins                                    -> nldi_basins_merged.gpkg
    fires   NIFC WFIGS perimeters discovered since SINCE (MTBS is local) -> nifc_perims.gpkg
    pairs   each fire on its own vs each basin, keep wildfires burning >= MIN_PCT of the basin
                                        -> fire_basin_pairs_5pct.csv, basins_5pct.gpkg, wildfires_5pct.gpkg
    ref     keep GAGES-II reference (CLASS == Ref) basins, then those < MAX_KM2
                                        -> ref_*_5pct.csv/.gpkg, ref_*_5pct_lt1000km2.csv/.gpkg

usage:
    python inputs/001_update_sites.py ref
    python inputs/001_update_sites.py sites basins merge fires pairs ref
"""

######## Warning!!! because we need to be polite to the USGS APIs with req to avoid 404s,
# the basins step takes multiple hours to run, do not run unless absolutely necessary
# (sites already in OUT are skipped, so a rerun only fetches what is missing)
# NIFC perimeters are live data, rerunning fires + pairs will not give the 2026-10-30 NIFC rows exactly

from __future__ import annotations
import io, json, sys, time, zipfile
from datetime import date
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import requests
import shapely

pd.set_option('display.max_columns', None)

REPO = Path(__file__).resolve().parents[1]

# ---------------------------------- CONFIG ----------------------------------
OUT = Path(r"C:\Users\A02306692\Downloads\nldi_basins")    # where the 2026-10-30 files live
STATES = r"C:\Users\A02306692\Downloads\PFF\release\Release_repo_in_progress\data\geoBoundaries-USA-ADM1_simplified.geojson"
MTBS = REPO / "inputs" / "mtbs_perims_DD.shp"
WESTERN = ["Washington", "Montana", "Wyoming", "Utah", "Oregon", "Idaho", "Colorado", "Arizona", "New Mexico", "California", "Nevada"]

NLDI_URL = "https://api.water.usgs.gov/nldi/linked-data/nwissite/USGS-{}/basin"
NIFC = ("https://services3.arcgis.com/T4QMspbfLg3qTGWY/arcgis/rest/services/"
        "WFIGS_Interagency_Perimeters/FeatureServer/0")
SINCE = date(2024, 12, 1)       # NIFC perimeters from here on, MTBS before
                                # NB: MTBS 2024 is incomplete, see report 2026-11-03 (fires Jun-Nov 2024 are missing from both)

EQUAL_AREA = 5070
MIN_PCT = 5                     # one fire must burn at least this % of the basin
WILDFIRE_TYPES = ["Wildfire", "WF"]   # MTBS Incid_Type, NIFC attr_IncidentTypeCategory

GAGESII_URL = ("https://www.sciencebase.gov/catalog/file/get/631405bbd34e36012efa304a"
               "?f=__disk__01%2F0a%2Fd2%2F010ad2fb628d42c764b81596fda93f0e86be478c")   # basinchar_and_report_sept_2011.zip
GAGESII_DIR = REPO / "inputs" / "gagesii_sept2011"
MAX_KM2 = 1000


# ---------------------------------- STEPS ----------------------------------

def step_sites(out: Path = OUT) -> pd.DataFrame:
    """every NWIS site with instantaneous discharge (00060) in the western states"""
    from pygeohydro import NWIS
    states = gpd.read_file(STATES)
    states = states[states["shapeName"].isin(WESTERN)]
    codes = states["shapeISO"].str.replace("US-", "").tolist()   # NWIS wants the bare 2-letter code
    queries = [{"stateCd": c, "parameterCd": "00060", "hasDataTypeCd": "iv"} for c in codes]
    sites = NWIS().get_info(queries).drop_duplicates("site_no")
    print(len(sites), "sites;", sites["huc_cd"].str[:2].value_counts().to_dict())   # 4,359 on 2026-10-29
    sites.drop(columns="geometry", errors="ignore").to_csv(out / "nwis_sites.csv", index=False)
    return sites


def fetch_basin(site: str, out: Path = OUT, pause: float = 4.6, max_wait: int = 900) -> str:
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


def step_basins(out: Path = OUT) -> pd.Series:
    """fetch (or reuse) every site's basin and record what happened to each"""
    ids = pd.read_csv(out / "nwis_sites.csv", dtype={"site_no": str})["site_no"].tolist()
    status = pd.Series({s: fetch_basin(s, out) for s in ids}, name="status")
    status.rename_axis("site_no").to_csv(out / "nldi_status.csv")
    print(status.value_counts().to_dict())     # 4,051 basins and 308 HTTP 404 on 2026-10-29
    return status


def step_merge(out: Path = OUT) -> gpd.GeoDataFrame:
    """one layer with every cached basin"""
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


def step_fires(out: Path = OUT, since: date = SINCE) -> gpd.GeoDataFrame:
    """every NIFC perimeter discovered since `since` inside the basins' bounding box, all fields"""
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


def step_pairs(out: Path = OUT) -> pd.DataFrame:
    """test each fire on its own against each basin, keep every wildfire that burned at least
    MIN_PCT of the basin. writes the pairs and the basins and fires that appear in them"""
    basins = gpd.read_file(out / "nldi_basins_merged.gpkg")
    mtbs = gpd.read_file(MTBS)
    nifc = gpd.read_file(out / "nifc_perims.gpkg")

    # make the two sources "mergable", every original column is kept for the fires file
    # (MTBS is not clipped to the western states: clipping would cut perimeters at state lines)
    mtbs = mtbs.assign(source="MTBS", fire_id=mtbs["Event_ID"], fire_name=mtbs["Incid_Name"],
                       fire_date=pd.to_datetime(mtbs["Ig_Date"]), fire_type=mtbs["Incid_Type"])
    nifc = nifc.assign(source="NIFC", fire_id=nifc["OBJECTID"].astype(int).astype(str), fire_name=nifc["poly_IncidentName"],
                       fire_date=nifc["attr_FireDiscoveryDateTime"], fire_type=nifc["attr_IncidentTypeCategory"])
    all_fires = pd.concat([mtbs.to_crs(EQUAL_AREA), nifc.to_crs(EQUAL_AREA)], ignore_index=True)
    all_fires = all_fires[all_fires["fire_type"].isin(WILDFIRE_TYPES)].reset_index(drop=True)
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

    print(f"{len(pairs)} fire-basin pairs; {len(f5)} fires; {len(b5)} basins")   # 2,763 / 1,134 / 1,471 on 2026-10-30
    return pairs


def gagesii_tables(folder: Path = GAGESII_DIR) -> pd.DataFrame:
    """GAGES-II (Falcone, Sept 2011) BasinID and Bas_Classif for the lower 48, one row per gage.
    downloads the release from ScienceBase once (55 MB) and keeps the two tables"""
    need = ["conterm_basinid.txt", "conterm_bas_classif.txt"]
    if not all((folder / n).exists() for n in need):
        folder.mkdir(parents=True, exist_ok=True)
        r = requests.get(GAGESII_URL, timeout=600, headers={"User-Agent": "firepeaks"})
        r.raise_for_status()
        with zipfile.ZipFile(io.BytesIO(r.content)) as outer:
            with zipfile.ZipFile(io.BytesIO(outer.read("spreadsheets-in-csv-format.zip"))) as inner:
                for n in need:
                    (folder / n).write_bytes(inner.read(n))
    read = lambda n: pd.read_csv(folder / n, dtype={"STAID": str}, encoding="latin-1")
    bid = read("conterm_basinid.txt")[["STAID", "STANAME", "DRAIN_SQKM", "HCDN-2009", "STATE"]]
    cls = read("conterm_bas_classif.txt")[["STAID", "CLASS", "AGGECOREGION", "HYDRO_DISTURB_INDX"]]
    return bid.merge(cls, on="STAID", validate="one_to_one")


def step_ref(out: Path = OUT, max_km2: float = MAX_KM2) -> pd.DataFrame:
    """keep the basins GAGES-II classes as reference (CLASS == Ref, not only HCDN-2009), then those
    smaller than max_km2 by their NLDI polygon area (the area burned % is computed against).
    a site not in GAGES-II has no class and is dropped"""
    pairs = pd.read_csv(out / "fire_basin_pairs_5pct.csv", dtype={"site_no": str, "fire_id": str})
    pairs["fire_date"] = pd.to_datetime(pairs["fire_date"], format="mixed")
    basins = gpd.read_file(out / "basins_5pct.gpkg")
    gii = gagesii_tables()

    b = basins.merge(gii, left_on="site_no", right_on="STAID", how="left").drop(columns="STAID")
    b["area_pct_diff_gii"] = 100 * (b["area_sqkm"] - b["DRAIN_SQKM"]) / b["DRAIN_SQKM"]
    print(f"{len(b)} basins with a >= {MIN_PCT}% wildfire: {int(b['CLASS'].notna().sum())} in GAGES-II "
          f"({b['CLASS'].value_counts().to_dict()}), {int(b['CLASS'].isna().sum())} not in GAGES-II")

    ref = b[b["CLASS"] == "Ref"].reset_index(drop=True)
    ref_pairs = pairs[pairs["site_no"].isin(ref["site_no"])].merge(
        ref[["site_no", "STANAME", "DRAIN_SQKM", "area_pct_diff_gii", "HCDN-2009", "AGGECOREGION", "HYDRO_DISTURB_INDX"]],
        on="site_no")
    ref.to_file(out / "ref_basins_5pct.gpkg", layer="basins", driver="GPKG")
    ref_pairs.to_csv(out / "ref_pairs_5pct.csv", index=False)

    small = ref[ref["area_sqkm"] < max_km2].reset_index(drop=True)
    small_pairs = ref_pairs[ref_pairs["site_no"].isin(small["site_no"])].reset_index(drop=True)
    tag = f"lt{int(max_km2)}km2"
    small.to_file(out / f"ref_basins_5pct_{tag}.gpkg", layer="basins", driver="GPKG")
    small_pairs.to_csv(out / f"ref_pairs_5pct_{tag}.csv", index=False)

    for name, bb, pp in [("reference", ref, ref_pairs), (f"reference < {max_km2:g} km2", small, small_pairs)]:
        print(f"{name}: {len(bb)} basins, {len(pp)} pairs, {pp.drop_duplicates(['source', 'fire_id']).shape[0]} fires "
              f"({pp['source'].value_counts().to_dict()})")
    return small_pairs


STEPS = {"sites": step_sites, "basins": step_basins, "merge": step_merge,
         "fires": step_fires, "pairs": step_pairs, "ref": step_ref}

if __name__ == "__main__":
    asked = sys.argv[1:]
    if not asked or any(s not in STEPS for s in asked):
        sys.exit(f"name the steps to run, in order, from: {' '.join(STEPS)}")
    for s in asked:
        print(f"== {s}")
        STEPS[s]()
