""" pull per-cell hourly AORC precip and store it as 
int16 parquet so that its much more memory efficient than full csvs

The summaries (max, mean, weighted mean, wet fraction, burned only) are computed 
from these files later """




from __future__ import annotations
import hashlib, os, threading, time, warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path
import numpy as np 
import pandas as pd 
import geopandas as gpd 

# CONFIG 
BUCKET = "noaa-nws-aorc-v1-1-1km"
VAR = "APCP_surface"
FILL = -32767         # missing value into the raw int16
SCALE = 0.1           # raw int16 * scale = mm
CHUNCK_LAT, CHUNK_LON = 128, 256   # zarr chunk shape (the time chuck is 144 h)
EA = 5070            # CONUS Albers of course 
STORE = Path("cache") / "aorc_cells"
_LOCK = threading.Lock()

@lru_cache(maxsize = 1)
def _s3(): 
    import s3fs
    return s3fs.S3FileSystem(anon = True) # one filesystem for the whole run 

@lru_cache(maxsize = 64)
def open_year(year: int): 
    """one year of AORC as raw int16 (mask_and_scale = False skips the float64 conversion"""
    import s3fs, xarray as xr 
    with warnings.catch_warnings(): 
        warnings.simplefilter("ignore")
        return xr.open_zarr(s3fs.S3Map(f"{BUCKET} /{year}.zarr", s3 =_s3()), 
                            consolidated=True, mask_and_scale=False) 


@lru_cache(maxsize = 1)
def grid(): 
    """cell-centre corrdinates read from the store, 
    the spacing is not exact, so dont compute from an origin"""
    ds = open_year(2000)
    return ds.latitude.values, ds.longitude.values

def basin_cells(geom, crs) -> pd.DataFrame: 
    """Every AORC cell whose outline touches the basin, 
    with frac share of the cell's area inside the basin (albers CONUS)
    cehter True if the cell centre is inside 
    col  grid index, column name used in yr files """
    from shapely import box 
    lat, lon = grid()
    dlat, dlon = np.diff(lat).mean(), np.diff(lon).mean()
    basin = gpd.GeoSeries([geom], crs = crs)
    minx, miny, maxx, maxy = basin.to_crs(4326).total_bounds
    i = np.where((lat + lon / 2 >= miny) & (lat - dlat / 2 <= maxy))[0]
    j = np.where((lon + dlon / 2 >= minx) & (lon - dlon / 2 <= maxx))[0]
    I, J = (a.ravel() for a in np.meshgrid(i, j, indexing="ij"))
    cy, cx = lat[I], lon[J]
    boxes = gpd.GeoSeries(box(cx - dlon / 2, cy - dlat / 2, cx + dlon / 2, cy + dlat / 2),
                          crs = 4326).to_crs(EA)
    b = basin.to_crs(EA).iloc[0]
    inter = boxes.intersection(b).area.to_numpy()
    keep = inter > 0 
    cells = pd.DataFrame({"i": I[keep], "j": J[keep], "lat": cy["keep"], "lon": cx[keep], 
                          "frac": inter[keep] / boxes.area.to_numpy()[keep]})
    pts = gpd.GeoSeries(gpd.points_from_xy(cells.lon, cells.lat), crs = 4326).to_crs(EA)
    cells ["centre"] = pts.within(b).to_numpy()
    cells["col"] = cells.i.astype(str) + "_" + cells.j.astype(str)
    return cells 

def register(name: str, geom, crs, store: Path = STORE) -> Path: 
    """Write cells.parquet for a basin, record name -> folder in basins.csv, 
    the folder is named by a hash of the cell set, sp a re-delinated basin gets a 
    new folder instead of the old basin's data"""
    cells = basin_cells(geom, crs) 
    cid = hashlib.sha1(cells[["i", "j"]].to_numpy("int64").tobytes()).hexdigest()[:12]
    folder = Path(store) / cid
    folder.mkdir(parents = True, exist_ok= True) 
    if not (folder / "cells.parquet").exists(): 
        cells.to_parquet(folder / "cells.parquet")
    index = Path(store) / "basins.csv"
    with _LOCK: 
        idx = pd.read_csv(index, dtype=str) if index.exists() else pd.DataFrame(columns = ["name", "cellset"])
        idx = idx[idx.name != str(name)]
        idx = pd.concat([idx, pd.DataFrame({"name": [str(name)], "cellset": [cid]})])
        idx.to_csv(index, index= False)
    return folder 

def folder_for(name: str, store: Path = STORE) -> Path: 
    """the folder for a basin name """
    idx = pd.read_csv(Path(store) / "basins.csv", dtype= str).set_index("name")
    return Path(store) / idx.loc[str(name), "cellset"]


def tile_groups(folders: list[Path]) -> list[list[Path]]:
    """group basins by the zarr tile their middle cell sits in, so neighbours share one fetch.
    a basin spanning two tiles still works, it just pulls the second tile itself"""
    groups: dict[tuple, list[Path]] = {}
    for f in folders:
        c = pd.read_parquet(f / "cells.parquet", columns=["i", "j"])
        key = (int(c.i.median()) // CHUNK_LAT, int(c.j.median()) // CHUNK_LON)
        groups.setdefault(key, []).append(f)
    return list(groups.values())


def pull_group_year(folders: list[Path], year: int) -> int:
    """fetch one year for a group of basins in a single read of the union of their cells,
    then split it into each basin's file. skips basins whose file already exists, writes
    through a .tmp so an interrupted run never leaves a half file. returns files written"""
    todo = [f for f in folders if not (f / f"{year}.parquet").exists()]
    if not todo:
        return 0
    cells = {f: pd.read_parquet(f / "cells.parquet", columns=["i", "j", "col"]) for f in todo}
    allc = pd.concat(cells.values()).drop_duplicates("col").sort_values(["i", "j"])
    da = open_year(year)[VAR]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        arr = da.isel(latitude=("cell", allc.i.to_numpy()),
                      longitude=("cell", allc.j.to_numpy())).values      # (hours, cells) int16
    pos = pd.Series(np.arange(len(allc)), index=allc.col.to_numpy())
    t = pd.DatetimeIndex(da.time.values, name="time_utc")
    for f, c in cells.items():
        df = pd.DataFrame(arr[:, pos[c.col].to_numpy()], index=t, columns=c.col.to_numpy())
        tmp = f / f"{year}.parquet.tmp"
        df.to_parquet(tmp, compression="zstd")
        os.replace(tmp, f / f"{year}.parquet")
    return len(todo)


def pull(folders: list[Path], years, workers: int = 4, retries: int = 3) -> list[tuple]:
    """pull every (tile group, year) in parallel. the work is network-bound, so threads help
    even on 4 cores. failures are retried with backoff, then returned rather than raised so
    one bad year doesn't stop a long run. rerun the same call to fill whatever failed"""
    jobs = [(g, y) for g in tile_groups(list(folders)) for y in years]

    def run(g, y):
        for k in range(retries):
            try:
                return pull_group_year(g, y)
            except Exception:
                if k == retries - 1:
                    raise
                time.sleep(2 ** k)

    failed, done, t0 = [], 0, time.time()
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(run, g, y): (g, y) for g, y in jobs}
        for fut in as_completed(futs):
            g, y = futs[fut]
            done += 1
            try:
                n = fut.result()
                print(f"[{done}/{len(jobs)}] {y} {len(g)} basins, {n} written  {time.time() - t0:.0f}s", flush=True)
            except Exception as e:
                failed.append((g, y, repr(e)))
                print(f"[{done}/{len(jobs)}] {y} FAILED {e!r}", flush=True)
    return failed


def read_precip(folder: Path, start, end, cols=None) -> pd.DataFrame:
    """hourly per-cell depth in mm (float32, NaN where AORC is missing) for [start, end] UTC"""
    start, end = pd.Timestamp(start), pd.Timestamp(end)
    parts = []
    for y in range(start.year, end.year + 1):
        p = Path(folder) / f"{y}.parquet"
        if not p.exists():
            raise FileNotFoundError(f"{p} not pulled yet")
        parts.append(pd.read_parquet(p, columns=cols).loc[start:end])
    raw = pd.concat(parts)
    return raw.where(raw != FILL).astype("float32") * np.float32(SCALE)


def read_cells(folder: Path) -> pd.DataFrame:
    """the cell table, indexed by column name, for weights (frac) and subsets (centre)"""
    return pd.read_parquet(Path(folder) / "cells.parquet").set_index("col")
