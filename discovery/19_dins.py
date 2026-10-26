"""Were structures in single-egress pockets destroyed more often in real fires? CAL FIRE DINS, northern-ca.

External data (Discovery prize only; not a challenge layer and never a model input): the CAL FIRE Damage Inspection
(DINS) database, FeatureServer POSTFIRE_MASTER_DATA_SHARE layer 0 on services1.arcgis.com (public, no key). One
record per inspected structure since 2013. From 2018 on the inspectors also record every undamaged structure in the
inspected area, so a destruction rate is defined only for incidents that started in 2018 or later.
Query: LATITUDE / LONGITUDE inside the northern-ca AOI bounding box and INCIDENTSTARTDATE >= 2018-01-01; fields
DAMAGE, INCIDENTNAME, INCIDENTSTARTDATE, STRUCTURECATEGORY, LATITUDE, LONGITUDE (no address, no parcel number).
Cached as discovery/external/dins/northern-ca_2018on.csv; --offline reads the cache only.

Design, fixed on 26 Oct 2026 before the data was downloaded:
  points    inside the northern-ca AOI polygon; DAMAGE 'Inaccessible' dropped
  exposure  the structure's nearest drivable Overture segment (<= 150 m, as for buildings in 01_islands.py) lies in
            a pocket of >= 50 buildings (02_egress.py, the headline threshold). Roads, not buildings: most destroyed
            structures are no longer drawn in Overture, so matching to buildings would drop them selectively
  outcome   destroyed = DAMAGE 'Destroyed (>50%)'
  estimate  destroyed share in pockets vs outside: Mantel-Haenszel risk ratio over strata
            (a) incident; (b) incident x tract (the same fire and the same tract). Strata need both groups
  CI        cluster bootstrap over incident x tract cells, B = 2000, numpy default_rng(20261026), percentile 95%
  subsets   residential structures (STRUCTURECATEGORY containing 'Residence'); without the Camp Fire (the largest)
Outputs: discovery/out/dins_pockets.csv (estimates), discovery/out/dins_by_incident.csv
Usage  : python discovery/19_dins.py [--offline]
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402

from common import (DISC, OUT, assign_points_to_tracts, con, load_tracts, log, region_outline,  # noqa: E402
                    xy_to5070)

REGION = "northern-ca"
URL = ("https://services1.arcgis.com/jUJYIo9tSA7EHvfZ/ArcGIS/rest/services/POSTFIRE_MASTER_DATA_SHARE/"
       "FeatureServer/0/query")
FIELDS = ["DAMAGE", "INCIDENTNAME", "INCIDENTSTARTDATE", "STRUCTURECATEGORY", "LATITUDE", "LONGITUDE"]
CACHE_CSV = DISC / "external" / "dins" / f"{REGION}_2018on.csv"
NEAR_MAX_M = 150.0
MIN_POCKET = 50
SEED, B = 20261026, 2000


def download(bbox) -> pd.DataFrame:
    x0, y0, x1, y1 = bbox
    where = (f"LATITUDE > {y0:.4f} AND LATITUDE < {y1:.4f} AND LONGITUDE > {x0:.4f} AND LONGITUDE < {x1:.4f} "
             "AND INCIDENTSTARTDATE >= timestamp '2018-01-01 00:00:00'")
    rows, off = [], 0
    while True:
        q = urllib.parse.urlencode({"where": where, "outFields": ",".join(FIELDS), "returnGeometry": "false",
                                    "orderByFields": "OBJECTID", "resultOffset": off, "resultRecordCount": 2000,
                                    "f": "json"})
        with urllib.request.urlopen(f"{URL}?{q}", timeout=120) as r:
            js = json.loads(r.read().decode("utf-8"))
        feats = js.get("features", [])
        rows += [f["attributes"] for f in feats]
        log(f"  DINS records {len(rows)}")
        if not feats or not js.get("exceededTransferLimit"):
            break
        off += len(feats)
        time.sleep(0.5)
    d = pd.DataFrame(rows, columns=FIELDS)
    CACHE_CSV.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(CACHE_CSV, index=False)
    return d


def mh_rr(a, n1, c, n0):
    """Mantel-Haenszel risk ratio over strata (arrays per stratum; strata without both groups contribute 0)."""
    N = n1 + n0
    ok = (n1 > 0) & (n0 > 0)
    num = np.sum(a[ok] * n0[ok] / N[ok])
    den = np.sum(c[ok] * n1[ok] / N[ok])
    return num / den if den > 0 else np.nan


def cell_table(d: pd.DataFrame) -> pd.DataFrame:
    """Per incident x tract cell: destroyed / all structures in pockets (a / n1) and outside (c / n0)."""
    g = d.groupby(["cell", "pocket"]).destroyed.agg(["sum", "count"]).unstack("pocket", fill_value=0)
    t = pd.DataFrame(index=g.index)
    for out, (stat, grp) in {"a": ("sum", True), "n1": ("count", True), "c": ("sum", False), "n0": ("count", False)}.items():
        t[out] = g[(stat, grp)].astype(float) if (stat, grp) in g.columns else 0.0
    t["inc"] = d.groupby("cell").incident.first().reindex(t.index)
    return t


def estimates(t: pd.DataFrame) -> dict:
    by_inc = t.groupby("inc")[["a", "n1", "c", "n0"]].sum()
    return dict(crude_rr=(t.a.sum() / t.n1.sum()) / (t.c.sum() / t.n0.sum()),
                mh_rr_incident=mh_rr(*(by_inc[k].values for k in ("a", "n1", "c", "n0"))),
                mh_rr_incident_tract=mh_rr(*(t[k].values for k in ("a", "n1", "c", "n0"))),
                destroyed_share_pocket=t.a.sum() / t.n1.sum(), destroyed_share_other=t.c.sum() / t.n0.sum(),
                n_pocket=t.n1.sum(), n_other=t.n0.sum())


def main(offline: bool):
    tr, tg = load_tracts(REGION)
    region_poly, pieces = region_outline(tg)
    if CACHE_CSV.exists():
        d = pd.read_csv(CACHE_CSV)
    elif offline:
        sys.exit(f"--offline but no cache at {CACHE_CSV}")
    else:
        aoi = json.loads((DISC.parent / "data" / "boundaries" / f"{REGION}-aoi.geojson").read_text())
        bbox = shapely.bounds(shapely.union_all([shapely.geometry.shape(f["geometry"]) for f in aoi["features"]]))
        d = download(bbox)
    log("DINS records (2018+, bounding box)", len(d))
    d = d[d.DAMAGE.notna() & (d.DAMAGE != "Inaccessible") & d.LATITUDE.notna() & d.LONGITUDE.notna()].copy()
    x, y = xy_to5070(d.LONGITUDE.values, d.LATITUDE.values)
    inside = shapely.contains_xy(region_poly, x, y)
    d, x, y = d[inside].copy(), x[inside], y[inside]
    log("inside the AOI, damage known", len(d))

    spec = importlib.util.spec_from_file_location("egress", Path(__file__).resolve().parent / "02_egress.py")
    eg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(eg)
    G = eg.build_graph(con(), REGION, "roads", None, "d", region_poly, pieces)
    (qq, tt), _ = shapely.STRtree(G["geoms"]).query_nearest(shapely.points(x, y), max_distance=NEAR_MAX_M,
                                                            return_distance=True, all_matches=False)
    seg = np.full(len(d), -1, np.int64)
    seg[qq] = tt
    att = eg.attach_buildings(G, x, y, seg)
    pu = np.where(att >= 0, G["pocket"][G["eu"][np.maximum(att, 0)]], -1)
    pv = np.where(att >= 0, G["pocket"][G["ev"][np.maximum(att, 0)]], -1)
    pk = np.where(pu >= 0, pu, pv)
    P = pd.read_parquet(OUT / f"{REGION}_pockets.parquet")
    psize = pd.Series(np.where(P.exit_outside_region, 0, P.n_bldg.values), index=P.pocket.values)
    size = np.where(pk >= 0, psize.reindex(np.maximum(pk, 0)).fillna(0).values, 0)
    d["attached"] = att >= 0
    d["pocket"] = size >= MIN_POCKET
    d["destroyed"] = d.DAMAGE.str.startswith("Destroyed")
    d["incident"] = d.INCIDENTNAME.str.strip() + " " + pd.to_datetime(d.INCIDENTSTARTDATE, unit="ms").dt.year.astype(str)
    d["tract"] = assign_points_to_tracts(x, y, tg)
    log("attached to a drivable segment", int(d.attached.sum()), "of", len(d), "| destroyed share attached / not",
        round(d.destroyed[d.attached].mean(), 3), round(d.destroyed[~d.attached].mean(), 3))
    d = d[d.attached & (d.tract >= 0)].copy()
    d["cell"] = d.incident + "|" + d.tract.astype(str)

    rng = np.random.default_rng(SEED)
    rows = []
    subsets = {"all structures": np.ones(len(d), bool),
               "residential": d.STRUCTURECATEGORY.fillna("").str.contains("Residence").values,
               "without the Camp Fire": ~d.incident.str.startswith("Camp ").values}
    for lab, m in subsets.items():
        s = d[m]
        t = cell_table(s)
        est = estimates(t)
        cells = t.index
        # bootstrap: cells drawn with replacement; a repeated cell counts as a separate stratum in (b)
        reps = pd.DataFrame([estimates(t.iloc[rng.integers(0, len(t), len(t))].reset_index(drop=True))
                             for _ in range(B)])
        for k, v in est.items():
            lo, hi = (np.nanpercentile(reps[k], [2.5, 97.5]) if k in reps and "rr" in k else (np.nan, np.nan))
            rows.append(dict(subset=lab, measure=k, value=v, ci_lo=lo, ci_hi=hi, n_structures=len(s),
                             n_incidents=s.incident.nunique(), n_cells=len(cells)))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "dins_pockets.csv", index=False)
    inc = d.groupby("incident").apply(lambda s: pd.Series(dict(
        n=len(s), n_pocket=int(s.pocket.sum()), destroyed_pocket=s.destroyed[s.pocket].mean(),
        destroyed_other=s.destroyed[~s.pocket].mean())), include_groups=False).sort_values("n", ascending=False)
    inc.to_csv(OUT / "dins_by_incident.csv")
    pd.set_option("display.width", 200)
    print(R.round(4).to_string(index=False))
    print(inc.head(15).round(3).to_string())


if __name__ == "__main__":
    main("--offline" in sys.argv[1:])
