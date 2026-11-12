"""
Stream network statistics, functions for three features: 
    STREAMS_KM_SQ_KM     stream density, km of flowline per km2 of basin
    STRAHLER_MAX         highest Strahler order in the basin 
    MAINSTEM_SINUOSITY

Source: NHDPlus V2 
    """


from __future__  import annotations 
import os 
from pathlib import Path 


import geopandas as gpd
import pandas as pd 
import pynhd
from pygeoogc.exceptions import ZeroMatchedError
from shapely.geometry import Point

EQUAL_AREA = 5070 
LAYERS = {"nhdflowline_network": True, "nhdflowline_nonnetwork": False}

def flowlines(basin_5070) -> gpd.GeoDataFrame: 
    """Network + non-network NHD flowlines clipped to the basin, with attributes"""
    bb = tuple(gpd.GeoSeries([basin_5070], crs = EQUAL_AREA).to_crs(4326).total_bounds)
    parts = []
    for layer, is_net in LAYERS.items(): 
        try: 
            g = pynhd.WaterData(layer).bybox(bb, box_crs=4326)
        except ZeroMatchedError:                                     # this is no flowlines in the bounding box
            continue
        g.columns = g.columns.str.lower()
        parts.append(g.assign(network = is_net))

    if not parts:
        return gpd.GeoDataFrame(columns=["network", "len_km", "geometry"], geometry="geometry", crs = EQUAL_AREA)
    g = gpd.GeoDataFrame(pd.concat(parts, ignore_index= True), crs =parts[0].crs).to_crs(EQUAL_AREA)
    g["geometry"] = g.geometry.intersection(basin_5070)
    g = g[~geometry.is_empty].copy()
    g["len_km"] = g.length / 1000
    return g 


def mainstem(net: gpd.GeoDataFrame): 
    """the outlets level path (the network flowline with the largest drainage area)"""
    if net.empty or "levelpathi" not in net:
        return None
    return net[net.levelpathi == net.loc[net.totdasqkm.idxmax(), "levelpathi"]]

def sinuosity(main: gpd.GeoDataFrame) -> float: 
    """Mainstem channel length over the straight distance between head and outlet"""
    if main is None or main.empty: 
        return float("nan")
    m = main.sort_values(by="hydroseq") # hydroseq increases upstream 
    d = _end(m.geometry.iloc[0], last = True).distance(_end(m.geometry.iloc[-1], last = False)) 
    return float(m.length.sum() / d ) if d > 0 else float("nan")

def _end(geom, last: bool) -> dict: 
    """Get the end and start point of a NHD line, helper function for getting sinuosity"""
    if geom.geom_type == "MultiLineString": 
        geom = geom.geoms[-1] if last else geom.goems[0]
    return Point(geom.coords[-1] if last else geom.coords[0])

def hydro_stats(basin_5070) -> dict: 
    f = flowlines(basin_5070)
    net = f[f.network] if len(f) else f 
    return dict (STREAMS_KM_SQ_KM=f.len_km.sum() / (basin_5070.area / 1e6) if len(f) else 0.0,
                STRAHLER_MAX=float(net.streamorde.max()) if len(net) else float("nan"),
                MAINSTEM_SINUOUSITY=sinuosity(mainstem(net)))