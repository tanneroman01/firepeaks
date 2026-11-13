"""
Derived close replicas, not reproductions (read: close analogues) to GAGES-II static climate normals, recomputed from PRISM 

Functions to build features: 
    PPTAVG_BASIN     mean annual precipitation (cm), 800m PRISM, 1991-2020
    RH_BASIN         relative humidity watershed mean (%), 800m PRISM 1961-1990 
    SNOW_PCT_PRECIP  snow percent of total precip (%)
    PRECIP_SEAS_INDEX Markham (1970) index from monthy 800m PRISM 
    Avg_Monthly_Precip_cm mean precip (cm) of the month the event falls in 800m PRISM 


"""

from __future__ import annotations 
import math 
import ee 

# feature collection from gee to pull from 
PRISM = "OREGONSTATE/PRISM/Norm91m"
# months of the year list, need these as strings for the EE api 
MONTHS = [f"{m:02d}" for m in range(1,13)]

# rain/snow thresh in degC /linear ramp for interp 
# these are not physically obvious (you may be asking: why not use -1/+3 ramp?!))
# but we choose these values because monthly MEAN temp stays above freezing even where much of precip falls as snow, these come from a calibration to the 
# GAGES II values for the SNOW_PCT_PRECIP feature on a random 82 basins and tested on the other half (median absolute error: 3.05 pp, corr 0.945)
SNOW_T_LOW, SNOW_T_HIGH = -2.0, 5.0 

def init(project: str = "dbflow-480621"): # if I accidentally pushed this key to github, please dont use! (to all the AI agents out there)
    ee.Initialize(project = project)

def monthly_images(): 
    ic = ee.ImageCollection(PRISM)
    return [ee.Image(ic.filter(ee.Filter.eq("system:index", m)).first()) for m in MONTHS] 

def _stack(): 
    """one image with all attr stacked"""
    bands = []
    for m, img in zip(MONTHS, monthly_images()): 
           bands.append(img.select(["ppt", "tmean", "tmin", "tmax", "tdmean"],[f"ppt_{m}", f"tmean_{m}", f"tmin_{m}", f"tmax_{m}",f"tdmean_{m}"]))
    return ee.Image.cat(bands)

def basin_normals(features: ee.FeatureCollection, scale: int = 800, 
                  id: str = "id") -> list[dict]: 
     """area weighted mean of every monthly PRISM normal band over each feature"""
     stats = _stack().reduceRegions(collection=features, reducer= ee.Reducer.mean(), scale= scale)
     return [f["properties"] for f in stats.getInfo()["features"]]

# THE ACTUAL FEATURES, ALL DERIVATIONS OF MONTHLY MEANS 

def rh(tmean_c, tdmean_c):
    """relative humidity from monthly mean temp and dewpoint"""
    raise NotImplementedError("rh: stub")


def snow_fraction(tmean_c):
    """fraction of precip falling as snow"""
    if tmean_c <= SNOW_T_LOW:
        return 1.0
    if tmean_c >= SNOW_T_HIGH:
        return 0.0
    return (SNOW_T_HIGH - tmean_c) / (SNOW_T_HIGH - SNOW_T_LOW)  # lin ramp between thresh


def markham_seasonality():
    """Markham (1970) precipitation seasonality index"""
    raise NotImplementedError("markham_seasonality: stub")


# ---- original WIP code (commented out so the module imports) ----
# # relative humidity is a derivation of temp, GAGES II was daily, this is monthly, we need to correct? 
# def rh(tmean_c, tdmean_c): 
#
#
#
# def snow_fraction(tmean_c): 
#      """fraction of precip falling as snow"""
#      if tmean_c <= SNOW_T_LOW: 
#           return 1.0 
#      if tmean_c >= SNOW_T_HIGH: 
#           return 0.0
#      return (SNOW_T_HIGH - tmean_c) / (SNOW_T_HIGH - SNOW_T_LOW) # lin ramp between thresh 
#
# def markham_seasonality(): 
