"""dfsadfd"""

import pandas as pd 
import numpy as np 
import geopandas as gpd
from pathlib import Path 
import yaml 
import firepeaks.flow as flow 
import firepeaks.pre_fire_peaks as pre 
import firepeaks.post_fire_peaks as post 






# paths and config
REPO = Path(__file__).resolve().parents[1]  # path object to which we add child items
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) #loads in defaults dict from yaml
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s02"].mkdir(parents = True, exist_ok = True)
C, F, FL = cfg["peaks"], cfg["fires"], cfg["filters"]


basins = gpd.read_file(P["s01"] / "basins.gpkg")
fires = pd.read_parquet(P["s01"] / "fires.parquet")
sites = pd.read_parquet(P["s01"] / "sites.parquet")
area = sites.set_index("gage_id").area_km2
snow = pd.read_csv(P["snow_table"], dtype= str).rename(columns=str.lower)
snow = pd.Series(pd.to_numeric(snow.snow_pct_precip, errors= "coerce").values, index = snow.gage_id.str.zfill(8)) 

match = pd.Timedelta
