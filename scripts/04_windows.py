""" Calculate a time of concentration at each gage in the dataset, locate the 
rise start of each event, define search windows as the rise start - Tc """


import pandas as pd 
import numpy as np 
import geopandas as gpd 
from pathlib import Path 
import firepeaks.window as window
import yaml 


# PATHS AND CONFIG READ IN FROM YAML 
REPO = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) #loads in defaults dict from yaml
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s04"].mkdir(parents = True, exist_ok = True)


peaks = pd.read_parquet(P["s02"])