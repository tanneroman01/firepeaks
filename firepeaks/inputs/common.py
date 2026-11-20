"""boiler plate functions used in inputs/ """
from __future__ import annotations 
from pathlib import Path
import numpy as np 
import geopandas as gpd 
import pandas as pd 

EQUAL_AREA = 5070 

def normalize_gage_id(g) -> str | None:
    """takes gage ids and converts to strings, removes decimals, etc. pandas will read in gage id's as float"""
    if g is None or (isinstance(g, float) and np.isnan(g)):  #dont pass nans
        return None
    s = str(g).strip().split(".")[0] # converts to string and splits, then indexes to the first item (removing decimals)
    return s.zfill(8) if len(s) < 8 else s

def study_basins(P: dict[str, Path]) -> gpd.GeoDataFrame: 
    """polygons of every site in paths.sites from paths.basins, EPSG 5070, one row per gage_id"""
    sites = pd.read_csv(P["sites"], dtype= str)
    b = gpd.read_file(P["basins"])
    b["gage_id"] = b["GAGE_ID"].map(normalize_gage_id)
    b = b[b.gage_id.isin(sites.gage_id)].to_crs(EQUAL_AREA)
    b["geometry"] = b.geometry.make_valid()
    b = b.dissolve("gage_id").reset_index()["gage_id", "geometry"]
    missing = set(sites.gage_id) - set(b.gage_id)
    if missing: 
        print(f"{len(missing)} gaged sites dont have a polygon, uh oh!")
    return b 


def exists(*files, overwrite: bool = False) -> bool: 
    """True when every file exists and overwrite is off, the caller skips its step if True"""
    if not overwrite and all(Path(f).exists() for f in files): 
        return True
    return False 

