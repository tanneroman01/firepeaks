import pandas as pd 
import numpy as np 
from post_fire_peaks import 
from pre_fire_peaks import _range_min
import geopandas as gpd 
from shapely.geometry import linemerge, unary_union
from shapely.geometry import Point 



"""define the window over which to search precip from basin attr using kirpatch Tc"""


#KIRPATCH EQUATION 
s1085 = max((r[z_]))
tc_kirpatch = 0.0667 * L ** 0.77 * s1085 

#------------------------- Basin Attributes for Tc calculations---------------------------------
def mainstem_line(main: gpd.GeoDataFrame, outlet: Point): 
    """the mainstem as one line"""
    geom = linemerge(unary_union(lis))


# find the rise start, fetch basin attr for kirpatch Tc calc 

def rise_start(qh: pd.Series, peak_time: pd.Timestamp, max_lookback: int = MAX )


def kirpatch(r: dict) -> dict: 
    L, A = r["main_km"]