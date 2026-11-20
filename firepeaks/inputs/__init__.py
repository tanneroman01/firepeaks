"""input builders, run in order by scripts/001_build_inputs.py. Each module exposes run, modules are only imported when 
the step runs
"""
import importlib

STEPS = {
    "candidates": "firepeaks.inputs.candidates",
    "unaltered": "firepeaks.inputs.unaltered",
    "aorc": "firepeaks.inputs.aorc", 
    "snow": "firepeaks.inputs.snow",
}

def step(name:str): 
    return importlib.import_module(STEPS[name]).run
