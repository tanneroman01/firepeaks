
##### FROZEN, ONLY RERUN IF YOU ARE ADDING MORE SITES 


from __future__ import annotations
from pathlib import Path
import numpy as np
import geopandas as gpd
REPO = Path(__file__).resolve().parents[1]
PFF = REPO.parent



def normalize_gage_id(g) -> str | None:
    """takes gage ids and converts to strings, removes decimals, etc. pandas will read in gage id's as float"""
    if g is None or (isinstance(g, float) and np.isnan(g)):  #dont pass nans
        return None
    s = str(g).strip().split(".")[0] # converts to string and splits, then indexes to the first item (removing decimals)
    return s.zfill(8) if len(s) < 8 else s

g = gpd.read_file(PFF / "GIS_geometries/Basins/gagesii_basins.gpkg")
sites = (g.assign(gage_id = g["GAGE_ID"].map(normalize_gage_id)).rename(columns= {"gii_km2": "drain_sqkm_gii"})
         [["gage_id", "drain_sqkm_gii"]].assign(site_id = lambda d: d.gage_id, source = "gagesii"))
(REPO / "inputs").mkdir(exist_ok=True)
sites.to_csv(REPO / "inputs/sites.csv", index = False)
print(len(sites), "sites!")
