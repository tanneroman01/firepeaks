"""Summarize precip over ____"""

from __future__ import annotations 
import numpy as np 
import pandas as pd  
from pathlib import Path 
from firepeaks.precip import summarize_files

# paths and config 
REPO = Path(__file__).resolve().parents[1]
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) #loads in defaults dict from yaml
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s03"].mkdir(parents = True, exist_ok = True)
