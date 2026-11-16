"""dfsadfd"""

import pandas as pd 
import numpy as np 
import geopandas as gpd
from pathlib import Path 
import yaml 
import firepeaks.pre_fire_peaks as pre 
import firepeaks.post_fire_peaks as post 
import firepeaks.flow as flow 


# paths and config
REPO = Path(__file__).resolve().parents[1]  # path object to which we add child items
cfg = yaml.safe_load((REPO / "config/default.yaml").read_text()) #loads in defaults dict from yaml
P = {k: (REPO / v).resolve() for k, v in cfg["paths"].items()} # takes key k and appends item v then adds to PATH object
P["s02"].mkdir(parents = True, exist_ok = True)
C, F, FL = cfg["peaks"], cfg["fires"], cfg["filters"]


fires = pd.read_parquet(P["s01"] / "fires.parquet")
sites = pd.read_parquet(P["s01"] / "sites.parquet")
area = sites.set_index("gage_id").area_km2


snow = pd.read_csv(P["snow_table"], dtype= str).rename(columns=str.lower)
snow = pd.Series(pd.to_numeric(snow.snow_pct_precip, errors= "coerce").values, index = snow.gage_id.str.zfill(8)) 

match = pd.Timedelta(days=C["nwis_match_days"]) # constant from config YAML 
now = pd.Timestamp.now(tz = "UTC") # for run logs 
pair_rows, peak_rows = [], []

for g, gf, in fires[fires.searched].groupby("gage_id"): 
    q = flow.load_iv(g, P["iv_cache"], C["record_start"], min(gf.search_end.max(), now), fetch = C["fetch_iv"]) 
    if q.dropna().empty: 
        print(f"{g}: no IV discharge on disk, skipped", flush= True)
        continue
    qh_all = post.hourly_max(q)
    ann_all = flow.annual_peaks(g, P["peak_cache"], fetch= C["fetch_iv"])
    s_pct = float(snow.get(g, np.nan))
    if not np.isfinite(s_pct): 
        print(f"{g}: now snow %, snowmelt filter not applied", flush= True)
    masked = fires[(fires.gage_id == g) & (fires.burned_pct >= F["mask_min_percent"])]     # earlier fires to mask

    for f in gf.itertuples():
        prior = masked.loc[masked.ig_date < f.ig_date, "ig_date"].tolist()
        qpre = pre.pre_fire_record(q, f.ig_date, prior, C["recovery_years"])
        th = pre.pot_threshold(qpre, C["pot_rate"], C["min_sep_h"], C["trough_ratio"], C["min_years"])
        ef = post.EventFilter(qpre, s_pct, FL["snow_pct_min"], FL["ac24_min"], FL["prom_min"])
        kw = dict(ef=ef, min_sep_h=C["min_sep_h"], trough_ratio=C["trough_ratio"])

        # pre-fire
        pk_pre = post.detect(qpre, th["q_thr"], qpre.index.min(), qpre.index.max(), qh=post.hourly_max(qpre), **kw)
        # post-fire
        t0, t1 = post.post_fire_period(f.ig_date, F["post_fire_yrs"])
        post_yrs = pre.record_years(q.loc[t0:t1])
        ev = post.detect(q, th["q_thr"], t0, t1, qh=qh_all, **kw)
        ev = ev.assign(gage_id=g, fire_id=f.fire_id, q_thr=th["q_thr"],
                       event_id=[f"{g}_{f.fire_id}_{t:%Y%m%d%H%M}" for t in ev.peak_time],
                       days_since_fire=(ev.peak_time - t0).dt.total_seconds() / 86400)
        peak_rows.append(ev)

        kept = ev[ev.excluded == ""]
        ann = ann_all[(ann_all.peak_dt >= t0) & (ann_all.peak_dt <= t1)]
        caught = sum(len(kept) > 0 and (kept.peak_time - a).abs().min() <= match for a in ann.peak_dt)
        pair_rows.append(dict(
            gage_id=g, fire_id=f.fire_id, incid_name=f.incid_name, ig_date=f.ig_date, burned_pct=f.burned_pct,
            next_fire_overlaps=f.next_fire_overlaps, area_km2=area.get(g, np.nan), snow_pct=s_pct, snowy=ef.snowy,
            prior_fires_masked=len(prior), q_thr=th["q_thr"], pot_rate=th["rate"], pre_years=th["record_years"],
            short_record=th["short_record"],
            pre_rate=(pk_pre.excluded == "").sum() / th["record_years"] if np.isfinite(th["q_thr"]) else np.nan,
            pre_removed_snowmelt=int((pk_pre.excluded == "snowmelt").sum()),
            pre_removed_minor_crest=int((pk_pre.excluded == "minor_crest").sum()),
            post_years_with_iv=post_yrs, n_peaks=len(kept),
            post_rate=len(kept) / post_yrs if post_yrs else np.nan,
            post_removed_snowmelt=int((ev.excluded == "snowmelt").sum()),
            post_removed_minor_crest=int((ev.excluded == "minor_crest").sum()),
            nwis_annual=len(ann), nwis_caught=caught,
            nwis_regulated=int(ann.peak_cd.fillna("").str.contains("6").sum())))
        print(f"{g} {f.fire_id}: thr {th['q_thr']:.0f} cfs, pre {th['record_years']:.1f} yr ({len(prior)} masked), "
              f"post {post_yrs:.2f} yr, {len(kept)} peaks ({len(ev) - len(kept)} excluded), "
              f"NWIS {caught}/{len(ann)}", flush=True)

pairs = pd.DataFrame(pair_rows)
peaks = pd.concat(peak_rows, ignore_index=True)
pairs.to_parquet(P["s02"] / "pairs.parquet")
peaks.to_parquet(P["s02"] / "peaks.parquet")
print(f"\n{len(pairs)} pairs ({pairs.short_record.sum()} short pre-fire records, no threshold), "
      f"{(peaks.excluded == '').sum()} peaks counted, {(peaks.excluded == 'snowmelt').sum()} snowmelt and "
      f"{(peaks.excluded == 'minor_crest').sum()} minor crest excluded; "
      f"NWIS annual peaks caught {pairs.nwis_caught.sum()}/{pairs.nwis_annual.sum()}")