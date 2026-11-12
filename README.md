# Post Fire Peak Flow Storm Pairing and Attribution

A pipeline for post-fire peak flow storm attribution. 


The core modules live in firepeaks/, these modules are importable and used extensively in the scripts/, which is the actual pipeline that can be run-end to end. 

### Steps in scripts/ 

01_sites.py 
- input: a site list (gage IDs plus polygons), a MTBS perimeter 
- per-site values: area, the AORC cells, and a hash of the geometry 
- time of concentration: what method? will need flow path and slope 
- output: a parquet with one row per site and polygons in a gpkg linked by site ID 

02_peaks.py 
- peaks table, two different configs could live here 
1. Read a CSV of peaks you already have (such as the ~2,000 Haley already identified), date only peaks will need a different time precision (such as for those derived from literature or field observation)
2. Pull NWIS 15 min flow for the period, detect peaks over a threshold, place peaks atleast ___ apart, find each events start by walking back down the rising limb (same idea as RREDI's event start)

03_precip.py
- for each site, pull the AORC per-cell hourly stack over the union of all the event windows plus the antecedent lookack 
- cache it keyed on the geometry hash for calc's in next steps 

04_windows.py 
- search each peak backwards, and output a parquet of event, storm startm storm end, and other info 

05_metrics.py 
- within each events window, compute every temporal metric under every aggregation, then the spatial metrics 
- output the metrics w/ suffixes defining the aggregation method 

06_percentiles.py 
- calculate the percentile features with the long-term aggregation 


```mermaid
Flowchart TD 


```

