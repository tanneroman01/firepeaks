"""Taking the 162 current gages and preparing them for peak identification"""

import pandas as pd 
import geopandas as gpd 
from pathlib import Path
from firepeaks import sites as S 
import yaml


# paths and constants config
REPO = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text())
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()}
P["s01"].mkdir(parents = True, exist_ok = True)

site_list = pd.read_csv(P["sites"], dtype={"gage_id": str, "site_id": str}) 
basins = S.load_basins(P["basins"], site_list, **cfg["sites"])

fires = S.fire_pairs(basins, P["mtbs"], wildfire_only= cfg["fires"]["wildfire_only"])
fires = S.search_periods(fires, cfg["fires"]["min_burned_pct"], cfg["fires"]["post_fire_yrs"])

chk, extra = S.check_against_events(fires, P["event_table"], basins.gage_id, cfg["fires"]["min_burned_pct"])
print(f"event-table pairs found: {chk.found.sum()}/{len(chk)}, "
      f"median burned % diff: {chk.pct_diff.abs().median(): .2f}, "
      f"{len(extra)} searched fires not in event table")

chk[~chk.found | (chk.pct_diff.abs() > 5)].to_csv(P["s01"] / "check_event_pairs.csv", index = False)
extra.to_csv(P["s01"] / "check_extra_fires.csv", index = False)

basins.to_file(P["s01"] / "basins.gpkg", driver="GPKG")
basins.drop(columns="geometry").to_parquet(P["s01"] / "sites.parquet")
fires.to_parquet(P["s01"] / "fires.parquet")
print (f"{len(basins)} basinss ({basins.area_flag.sum()} area -flagged ), "
       f"{len(fires)} gage-fire pairs, {fires.searched.sum()} searched") 
