"""call builders in inputs/ in sequence, getting all sites and the needed data for next steps. 
    Warning!!! because we need to be polite to the USGS APIs with req to avoid 404s,
# the basins step takes multiple hours to run, do not run unless absolutely necessary
# (sites already in OUT are skipped, so a rerun only fetches what is missing)
# NIFC perimeters are live data, rerunning fires + pairs will not give the NIFC rows exactly"""

##### FROZEN, ONLY RERUN IF YOU ARE ADDING MORE SITES 

from __future__ import annotations 
import argparse
import sys
import time 
from firepeaks.config import load
from firepeaks.inputs import STEPS, step 

def main(argv = None) -> None: 
    ap = argparse.ArgumentParser(description=__doc__, formatter_class = argparse.RawDescriptionHelpFormatter)
    ap.add_argument("steps", nargs="*", help= f"any of {', '.join(STEPS)}") 
    ap.add_argument("--config", default= None)
    ap.add_argument("--list", action = "store_true", help= "print the plan and exit")
    a = ap.parse_args(argv)
    cfg, P = load(a.config)
    plan = a.steps or [s for s, on in cfg["inputs"]["run"].items() if on]
    bad = [s for s in plan if s not in STEPS]
    if bad: 
        sys.exit(f"unkown step(s) {bad}; choose from {list(STEPS)}")
    print("plan:", "->".join(plan) if plan else "(nothin stwich on in inputs.run)")
    if a.list: 
        return
    t_all = time.time()
    for s in plan: 
        t = time.time()
        print(f"== {s} done ({time.time() - t:.0f}s)", flush = True)
    print(f"all finished ({time.time() - t_all:.0f}s)")

if __name__ == "__main__": 
    main()
    

