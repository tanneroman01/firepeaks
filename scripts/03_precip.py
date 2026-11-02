"""Summarize precip over ____"""

import numpy as np 
import pandas as pd 
from __future__ import annotations 
from pathlib import Path 
from firepeaks.precip import summarize_files

# ______________ CONFIG_________________________________



for ws in ["DS", "US"]: 
    summarize_files()