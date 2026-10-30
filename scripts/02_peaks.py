"""dfsadfd"""

import pandas as pd 
import numpy as np 
import geopandas as gpd
import firepeaks.pre_fire_peaks as pre
import firepeaks.post_fire_peaks as post
from pathlib import Path 
import yaml 
import firepeaks.flow as flow 

#urlS 
IV_URL   = "https://nwis.waterservices.usgs.gov/nwis/iv/"
DV_URL   = "https://waterservices.usgs.gov/nwis/dv/"
SITE_URL = "https://waterservices.usgs.gov/nwis/site/"
UA = {"User-Agent": "firepeaks"}

# paths and config
REPO = Path(__file__).resolve().parents[1]  # path object to which we add child items
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) #loads in path dict from yaml
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s02"].mkdir(parents = True, exist_ok = True)

basins = gpd.read_file(P["s01"] / "basins.gpkg")
fires = pd.read_parquet(P["s01"] / "fires.parquet")
sites = pd.read_parquet(P["s01"] / "sites.parquet")


