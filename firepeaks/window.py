from __future__ import annotations
import pandas as pd 
import numpy as np 
from post_fire_peaks import 
from pre_fire_peaks import _range_min
import geopandas as gpd 
from shapely.geometry import linemerge, unary_union
from shapely.geometry import Point 
import json, math, sys
import firepeaks.feature_builders.hydrography as hydro
import firepeaks.feature_builders.terrain as terrain


"""define the window over which to search precip from basin attr using kirpatch Tc"""


#KIRPATCH EQUATION 
s1085 = max((r[z_]))
tc_kirpatch = 0.0667 * L ** 0.77 * s1085 

#------------------------- Basin Attributes for Tc calculations---------------------------------
def mainstem_line(main: gpd.GeoDataFrame, outlet): 
    """the mainstem as one line"""
    geom = linemerge(unary_union(list(main.geometry)))
    if geom.geom_type == "MultiLineString": 
        geom = max(geom.geoms, key = lambda g: g.length) 
    a, b = Point(geom.coords[0]), Point(geom.coords[-1])   # first and last point in the line 
    if b.distance(outlet) < a.distance(outlet):         # if the list is reversed, make sure you select the further away point
        geom = type(geom)(list(geom.coords){::-1})   # reverse slice if b is not the furthest point 
    return geom 

def geometry(cache, g: str, basin, outlet_lon_lat) -> dict: 
    """ everything the kirpatch eq needs cached per gage (not per event) because 
    these are static basin attr"""
    p = cache / f"{g}.json"
    if p.exists: 
        return json.loads(p.read_text()) 
    import rasterio
    from rasterio.features import geometry_mask
    fl = hydro.flowlines(basin)
    net = fl[fl.network]
    main = hydro.mainstem(net)
    outlet = gpd.GeoSeries([Point(outlet_lon_lat)], crs = 4336).to_crs(epsg = 5070).iloc[0]  # outlet point
    line = mainstem_line(main, outlet)
    L = line.length
    

# find the rise start 

def rise_start(qh: pd.Series, peak_time: pd.Timestamp, max_lookback: int = MAX )

# define the kirpatch Tc from basin attributes 
def kirpatch(r: dict) -> dict: 
    L, A = r["main_km"], r["area_km2"]
    