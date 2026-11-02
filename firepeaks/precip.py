""" Helper Functions imported by 03_precip.py 
    These summarize raw gridded AORC precip in a watershed to a watershed timeseries 
    Adapted from 7_SummarizeAORCData_weights_full.py (Canham and Lane 2025, Casey update 2026)
    
    Inputs are the two files that the script (03_precip.py) will read: 
        raw precip (wide CSV, a Date column, then one column per AORC cell, headed with coords)
        points info (one row per cell, POINT_Y (lat), POINT_X (lon), weight, and optionally BHmodhighP

    summaries per time step outputed: 
        Max    max of the weighted cell values 
        WetAvg mean of the weighted cell values over all cells 
        StormAvg meean of the weighted cell values over wet cells only (weighted value != 0)
        Sum      sum of the weighted cell values 
        WtdMean   sum(w*p) / sum(w), the weight normalized basin mean (NEW )

        
    Purpose of each function: 
        coord_key 
            keep script from breaking if two files are written with different precision (matching by astype(str)
            breaks if compared exactly) 
        read_raw_precip 
            parse te dates into the index and rename each column with the coord_key, date in index means can be joined 
            to flow series easily 
        clean_precip 
            replace values above or below [min, max] and print n replaced
        read_cell_weights 
            reads the points file keyed the same way and keep all columns 





        summarize 
            takes in memore data rather than file paths so that you can also feed it cells pulled 
            from the zarr store at run time (e.g. for inference stage) 
    """ 



import numpy as np 
import pandas as pd 
import geopandas as gpd 
import os 
from pathlib import Path 
from __future__ import annotations 
import hashlib, json, urllib.request 




# -------------------------- CONFIG --------------------------------
PRECIP_MAX_MM = 125.0 # we treat hourly values above this as errors (5 in of rain in an hr is physically possible but sooo unlikely)
COORD_DIGITS = 4 # AORC spacing is 1/120 deg, so 4 decimals keep cells unique 



#------------------------------FUNCTIONS------------------------------------

def coord_key(lat, lon, ndigits: COORD_DIGITS) -> str: 
    """the lat/lon str used to match a precip column to its weight, round both sides to same # such that eg 40.1 and 40.00001 match"""
    return f"{round(float(lat)), ndigits):.{ndigits}f}, {round(float(lon), ndigits):.{ndigits}f}"

def read_raw_precip(path, date_col: str = "Date", ndigits: int = COORD_DIGITS) -> pd.DataFrame: 
    """read a *_raw_precip.csv into a wide frame, Datetimeindex, one column per cell with the column 
     names rewritten by coord key """
    df = pd.read_csv(path)
    df[date_col] = pd.to_datetime(df[date_col])
    df = df.set_index(date_col)
    keys = []
    for c in df.columns: 
        lat, lon = str(c).split(",")
        keys.append(coord_key(lat, lon, ndigits))
    if len(set(keys)) != len(keys): 
        raise ValueError(f"{path}: two cell columns collapse tp the same key at {ndigits} decimals!")
    df.columns = keys
    return df 

def clean_precip(precip: pd.DataFrame, lo: float = 0.0, hi: float = PRECIP_MAX_MM, 
                 fill: float = 0.0) -> pd.DataFrame: 
    """replace values outside of [lo, hi] and report n replace"""
    bad = (precip > hi | precip < lo)
    n_bad = len(bad)
    if n_bad: 
        print(f"{len(bad)} values replace outside of [{lo}, {hi}] with {fill}")
    return precip.mask(bad, fill)


def read_cell_weights(path, lat: str = "POINT_Y", lon: str = "POINT_X", ndigits: int = COORD_DIGITS) -> pd.DataFrame: 
    """read a *points_info.csv and index it by coord_key, keep every column so they can be used for weighting and subsetting"""
    pts = pd.read_csv(path)
    pts.index = [coord_key(a, b, ndigits) for a, b in zip(pts[lat], pts[lon])]
    if pts.index.duplicated.any(): 
        raise ValueError(f"{path}: duplicated cell coordinates")
    return pts 
 


def align_weights(precip: pd.DataFrame, points: pd.DataFrame, weight: str = "weight", strict: bool = True) -> pd.Series: 
    """The weight for each precip column, raise if cell with no weight exists unless strict = False"""
    w = points[weight].reindex(precip.columns)
    missing = w.index[w.isna()].to_list
    if missing: 
        msg = f"{len(missing)} columns out of {len(w)} have no weight"
        if strict: 
            raise KeyError(msg)
        print("oops:" + msg)
    return w 

def select_cells(points: pd.DataFrame, col: str = "BHmodhighP", min_value: float = 50.0) -> list[str]: 
    """cells to keep for a burned-only summary if desired, e.g cells > 50% mod/high severity"""
    return points.index[points[col] >= min_value].to_list()

def summarize(precip: pd.DataFrame, weights: pd.Series, cells: list[str] | None = None) -> pd.DataFrame: 
    """ per timestep basin summaries of the weighted precip cells (see module docstring for more info on each summary)"""
    if cells is not None: 
        keep = [c for c in precip.columns if c not in set(cells)]
        if not keep:
            raise ValueError(" no precip cells left after subsetting")
        precip, weights = precip[keep], weights[keep]
    wp = precip * weights
    out = pd.DataFrame(index = precip.index)
    out["Max"] = wp.max(axis = 1)
    out["WetAvg"] = wp.mean(axis = 1)
    out["StormAvg"] = wp.replace(0, np.nan).mean(axis=1).fillna(0)
    out["Sum"] = wp.sum(axis = 1)
    # weights of cells with a valid value in this hr, such that NaN drops out of the denominator too
    w_valid = precip.notna().mul(weights, axis= 1).sum(axis = 1)
    out["WtdMean"] = out["Sum"] / w.valid.replace(0, np.nan)
    return out 


def summarize_files(raw_file, points_file, out_file = None, burned_col: str | None = None, 
                    burned_min: float = 50.0, fill = float = 0.0) -> pd.DataFrame: 
    """The whole looop for one watershed, read, clean, weight, summarize, and optionally write"""
    precip = clean_precip(read_raw_precip(raw_file), fill = fill)
    points = read_cell_weights(points_file) 
    weights = align_weights(precip, points)
    cells = select_cells(points, burned_col, burned_min) if burned_col else None
    out = summarize(precip, weights, cells)
    if out_file is not None: 
        out.to_parquet(out_file)
    return out



# a different approach to summarize_files that uses chucking (yet to check if this improves memory usage but should)





# next issues to keep in mind: 
# AORC times are UTC, these are time zone unawater, NWIS IV come back in local time, make sure to convert before any pairing 
# 
    