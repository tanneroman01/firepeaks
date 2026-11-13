""" functions for storm features, following the RREDI definitions

StormTotalPrecipmm   sum of hourly depths over the storm 
StormDuration_h      span of the storm's tip series 
StormPeakInt_mmh     MI60
AntecedentPrecip_mm 

These features are aggregated by the per-hour spatial maximum and the per-gour spatial minimum 


"""


from __future__ import annotations 
import hashlib, warnings 
from pathlib import Path 
import numpy as np 
import pandas as pd 


BUCKET = "noaa-nws-aorc-v1-1-1km" 
VAR = "APCP_surface"
LAT0, LON0, DEG = 20.0, -130.0, 1.0 / 120.0 
TIP_MM = 0.1
TIP_ROUNDING = "round" # or "floor". This governs how an hourly depth becomes 
# a whole number of tips. StormGenerator breaks a storm 
# where the gap between consecutive tips exceeds 3h
# because the last tip of a wet hr lands on the hr boundary 
# and the first tip of the next wet hr lands 3600/N in, the gap is 
# (dry_hrs + 1/N) hrs. so: 
# dry_hrs = 3 -> gap = 3 + 1/N 

# this is all necessary because parts of RREDI were originally written for tipping bucket gages, so to keep 
# storm definitions consistent we treat each hrs total as tips spread evenly throughout that hr, e.g an hr with 0.4 mm becomes
# N = 4 tips 
# if we just wrote dry_hrs > 3m sorms separated by exactly 3 hrs would be merged 
# with "round" an hr with less than 0.05mm becomes zero tips 


STORM_GAP_H = 6   # this is the minimum interum time MIT, 3 or 6 hrs?  
MIN_DEPTH_MM = 1.0 
MIN_DURATION_H = 1.0 
ANTECEDENT_DAYS = 6 
CACHE = Path(__file__).resolve().parents[1] / "cache" / "aorc"
WET_AREA_FRAC = 0.05 # a fraction of cells that have precip for an hour to count as wet? this was defined as 
# 5% in the manuscript from Haley and Belize, should we keep? 


def _open(year: int): 
    """pull the AORC data for a basin, via NOAA's public AWS bucket
        one zarr store per year (this function needs network access if called)"""
    import s3fs, xarray as xr
    fs = s3fs.S3FileSystem(anon = True) # because the bucket is public, no key needed 
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore")
        return xr.open_zarr(s3fs.S3Map(f"{BUCKET}/{year}.zarr", s3 =fs), consolidated= True) 

    
def cells_in_basin(geom_4326, max_cells: int = 4000, mode: str = "centre"):
    """indices of the AORC cells to summarize over, takes all the cells whose
    centre falls inside the basin, not every cell whose footprint intersects the basin"""
    raise NotImplementedError("cells_in_basin: stub")


def watershed_series(geom_4326, start, end, tag: str) -> pd.DataFrame:
    """Hourly basin summary in both mean and max depth, wetfrac is
    the fraction of cells with precip"""
    raise NotImplementedError("watershed_series: stub")


def watershed_max_series(geom_4326, start, end, tag: str) -> pd.Series:
    """Hourly per-cell MAXIMUM precipitation (mm) over the basin, for [start, end]"""
    raise NotImplementedError("watershed_max_series: stub")


# ---- original WIP code (commented out so the module imports) ----
# def cells_in_basin(goem_4326, max_cells: int = 4000, mode: str = "centre"): 
#     """ indices of the AORC cells to summarize over, takes all the cells whose 
#     centerfalls inside the basin, not every cell whose footprint intersects the basin"""
#     from shapely.geometry import Point, box 
#     minx, miny, maxx, maxy = geom_4326.bounds
#     i0 = max(0, int(np.floor((miny - Lat0) / DEG))) 
#     i1 = 
#     j0 = 
#     j1 =
#     if (i1-i0 + 1) * (j1 - j0 +1 ) > max_cells: 
#         raise ValueError("basin bounding box to big for a single extraction")
#     h = DEG/ 2.0 
#     ii, jj = [], []
#     for i in range(j0, j1 + 1): 
#
#
#
#
#
# def watershed_series(geom_4326, start, end, tag: str) -> pd.DataFrame: 
#     """Hourly mazin summary in both mean and max depth, wetfrac is 
#     the fractiom of cells with precip"""
#     CACHE.mkdir(parents=True, exist_ok=True)
#     import warnings as _w
#     ii,jj, n = cells_in_basin(geom_4326, mode=mode)
#     # the cells set is part of the cache key
#     cells = hashlib.sha1(np.asarray([ii, jj]).tobytes()).hexdigest()[:10]
#     key = (CACHE / f)
#
#
#
#     def watershed_max_series(geom_4326, start, end, tag: str) -> pd.Series: 
#         """Hourly per-cell MAXIMUM precipitation (mm) over the basin, for [start, end]"""
#         CACHE.mkdir(parents=True, exist_ok= True)
#         key = CACHE / f"{tag}_{pd.Timestamp(start).date()}_{pd.Timestamp(end).date()}.csv"
#         if key.exists(): 
#             s = pd.read_csv(key, index_col=0, parse_dates= True).iloc[:, 0]
#             s.name = "mm"
#             return s
#         ii, jj, n = cells_in_basin(geom_4326, mode=mode)
#         parts = []
#         for year in range(pd.Timestamp(start).year, pd.Timestamp(end).year + 1): 
#             ds = _open(year)
