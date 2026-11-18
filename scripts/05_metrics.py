"""take the windows generated in step 4, aggregate precip over the windows 
and calculate all precip related metrics for modelling."""

import pandas as pd 
import numpy as np 
import geopandas as gpd 
from pathlib import Path
import yaml 
import firepeaks.aggregate as agg 
# paths and config 

REPO = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) 
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s05"].mkdir(parents = True, exist_ok = True)
C, F, FL = 

windows = pd.read_parquet(P["s04"] / "windows.parquet")
precip = 