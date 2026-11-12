""" Grid-based static features via the USGS NHDPlus V2 catchement 
attribute releases 
    PET          mean annual potential ET, mm 
    WD_BASIN     mean annual days of measureable precip 
    BFI_BASIN    base-flow index, % 
    TOPWET       topographic wetness index 
    PERMAVE      STATSGO average permeability 
    CLAYAVE, SILTAVE, SANDAVE STATSGO % 


"""

from __future__ import annotations
import io 
import zipfile 
from pathlib import Path 
import geopandas as gpd 