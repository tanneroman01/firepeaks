"""fill out the rest of the event table, with all features, after runnign this, you can move on to modelling!"""

from pathlib import Path 
import numpy as np 
import geopandas as gpd 
import pandas as pd 
import yaml 
import firepeaks.feature_builders.terrain as terrain
import firepeaks.feature_builders.burn as burn 
import firepeaks.feature_builders.climate as climate 
import firepeaks.feature_builders.catchment as catchment 
import firepeaks.feature_builders.hydrography as hydro 
import firepeaks.feature_builders.veg as veg