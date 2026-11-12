"""
Basin Shape and terrain statistics, proxy for GAGES-II BasinID / Topo statics 

This is not a standalone script, these are a set of functions to be imported to build_statics.py as part of the feature building pipeline. 

Every function takes a basin polygon in EPSG:5070 (equal area) and returns plain numbers (in meters usually)


"""

from __future__ import annotations # allows us to name a class before its defined, newer python has type hints native, but this helps if users run on pthron 3.6 or older 

import numpy as np 
from rasterio import features 
from rasterio.transform import from_origin 

# Contants 
EQUAL_AREA = 5070 
COMPACT_CELL_M = 30.0 
DEM_ASSET = "USGS/3DEP/10m_collection" # the DEM asset to use for terrain statistics, this is a 10m resolution DEM from USGS 3DEP collection
ALBERS_N = 0.6029 # the grid convergence for the Albers projection, used to correct aspect calculations
ALBERS_LON = -96.0 # the central meridian for the Albers projection, used to correct aspect calculations
TERRAIN_SCALE_M = 10 # the scale to use for terrain statistics, this is the resolution of the DEM in meters

def compactness(basin_5070) -> float:
    perim = basin_5070.length
    return float(basin_5070.area / perim ** 2 * 100) if perim else float("nan")

def dem(): 
    """Return the DEM asset to use for terrain statistics."""
    import ee
    return (ee.ImageCollection(DEM_ASSET).select("elevation").mosaic()
        .setDefaultProjection("EPSG:5070", None, 10)
        .reproject(crs="EPSG:5070", scale=TERRAIN_SCALE_M)) 

def terrain_stats(basin_5070, outlet_lonlat: tuple[float, float]) -> dict:
    """ seven GAGES II like terreain statistics for one basin, need ee.Initialize() first. 
        'outlet_lonlat' is the gage or pour point (WSG84 lon, lat): it gives gives the ELEV_SITE_M features and the longitude 
        for the aspect covergence correction"""
    import math
    import ee 
    import geopandas as gpd 
    geom = gpd.GeoSeries([basin_5070.simplify(30)], crs = EQUAL_AREA).to_crs(4326).iloc[0]
    region = ee.Geometry(geom.__geo_interface__)
    demm = dem()
    slope = ee.Terrain.slope(demm)
    asp = ee.Terrain.aspect(demm).multiply(math.pi / 180)
    img = ee.Image.cat([demm.rename("elev"), slope.multiply(math.pi /180).tan().multiply(100).rename("slope"), asp.sin().rename("s"), asp.cos().rename("c")])
    red = ee.Reducer.mean().combine(ee.Reducer.min(), sharedinputs = True).combine(ee.Reducer.max(), sharedInputs = True)
    r = img.reduceRegion(red, region, scale = TERRAIN_SCALE_M, crs = "EPSG:5070", maxPixels = 1e9, tileScale = 4).getInfo()
    lon, lat = outlet_lonlat
    site = demm.reduceRegion(ee.Reducer.first(), ee.Geometry.Point([lon, lat]), scale = TERRAIN_SCALE_M, 
                            crs = "EPSG:5070").getInfo().get("elevation")
    ang = math.atan2(r["s_mean"], r["c_mean"]) - math.radians(ALBERS_N * (lon - ALBERS_LON))
    lo, hi, mu = r["elev_min"], r["elev_max"], r["elev_mean"]
    return dict(ELEV_MEAN_M_BASIN = mu, ELEV_MAX_M_BASIN = hi, ELEV_SITE_M = site, 
                RRMEAN = (mu - lo) / (hi-lo) if hi > lo else float("nan"), 
                SLOPE_PCT = r["slope_mean"], ASPECT_NORTHNESS = math.cos(ang), ASPECT_EASTNESS = math.sin(ang))

# Defining the features in the above dict: 
# ELEV_MEAN_M_BASIN is the average elevation of the whole basin 
# ELEV_MAX_M_BASIN is the highest point in the basin 
# ELEV_SITE_M is the eleveation at the gages 
# RRMEAN is the relief ration (on 0-1)
# SLOPE_PCT is the avg slope, slope per pixel in degrees converted to percent with tan(slope) * 100
# ASPECT NORTHNESS is how much the basin, on average, faces north 
# ASPECT EASTNESS is how much the basin, on average, faces east 

