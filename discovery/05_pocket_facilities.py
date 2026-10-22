"""Critical facilities that sit behind a single road segment (inside a single-egress pocket).

Rebuilds the drivable graph exactly as 02_egress.py does (same functions), attaches Overture POIs of three
groups to the graph edge they sit on (nearest drivable segment within 150 m) and flags those inside a pocket
of >= 50 buildings whose only link to the meshed network is one segment.

  schools_childcare : elementary/middle/high/public/private/charter school, school, preschool, day care
  care_homes        : retirement_home, skilled_nursing, assisted_living_facility, hospice
  emergency         : fire_department, ambulance_and_ems_services, police_department, hospital, emergency_room

POIs are filtered to confidence >= 0.5 and de-duplicated on (lower-cased name, 200 m grid cell).
Usage: python discovery/05_pocket_facilities.py          (all regions; needs discovery/cache/<region>_bldg.npz)
Outputs: discovery/out/pocket_facilities.csv            one row per facility inside a pocket >= 50 buildings
         discovery/out/pocket_facilities_summary.csv    counts / shares by group x wildfire-hazard quartile
"""
from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import shapely

from common import CACHE, OUT, REGIONS, assign_points_to_tracts, con, load_tracts, log, p, region_outline, xy_to5070

eg = importlib.import_module("02_egress")

GROUPS = {
    "schools_childcare": ["elementary_school", "middle_school", "high_school", "school", "public_school",
                          "private_school", "charter_school", "preschool", "day_care_preschool",
                          "child_care_and_day_care"],
    "care_homes": ["retirement_home", "skilled_nursing", "assisted_living_facility", "hospice"],
    "emergency": ["fire_department", "ambulance_and_ems_services", "police_department", "hospital", "emergency_room"],
}
CAT2GROUP = {c: g for g, cs in GROUPS.items() for c in cs}
MIN_POCKET_BLDG = 50


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    ph = k / n
    den = 1 + z * z / n
    centre = (ph + z * z / (2 * n)) / den
    half = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return centre - half, centre + half


def run_region(region: str, whp_cut) -> pd.DataFrame:
    c = con()
    tr, tg = load_tracts(region)
    region_poly, pieces = region_outline(tg)
    st = pd.read_parquet(p(region, "strata"), columns=["GEOID", "usfs_WHP_mean", "aiannh_name"])
    st["GEOID"] = st.GEOID.astype(str)
    tr = tr.merge(st, on="GEOID", how="left")
    G = eg.build_graph(c, region, "roads", None, "d", region_poly, pieces)
    B = np.load(CACHE / f"{region}_bldg.npz")
    bx, by, bseg = B["x"].astype(np.float64), B["y"].astype(np.float64), B["seg"]
    batt = eg.attach_buildings(G, bx, by, bseg)
    pu = np.where(batt >= 0, G["pocket"][G["eu"][np.maximum(batt, 0)]], -1)
    pv = np.where(batt >= 0, G["pocket"][G["ev"][np.maximum(batt, 0)]], -1)
    bpocket = np.where(pu >= 0, pu, pv)
    psize = np.bincount(bpocket[bpocket >= 0], minlength=int(G["pocket"].max()) + 1)
    # pockets whose only link is the virtual OUT edge (exit leaves the extract) are excluded, as in 02
    u, v = G["u"], G["v"]
    pc = np.where(G["pocket"][u] >= 0, G["pocket"][u], G["pocket"][v])
    att = (G["core"][u] ^ G["core"][v]) & (pc >= 0) & ~G["virtual"]
    has_real_exit = np.zeros(len(psize), bool)
    has_real_exit[np.unique(pc[att])] = True
    egress_edge = pd.Series(np.where(att)[0]).groupby(pc[att]).first()
    names = c.execute(f"SELECT s.si, r.names.primary AS name FROM d_s s JOIN '{p(region, 'roads')}' r USING(id) "
                      "ORDER BY s.si").fetchdf()

    cats = ",".join(f"'{k}'" for k in CAT2GROUP)
    poi = c.execute(f"""SELECT id, names.primary AS name, categories.primary AS cat, confidence,
                               addresses[1].locality AS locality, ST_X(geometry) x, ST_Y(geometry) y
                        FROM '{p(region, 'pois')}' WHERE categories.primary IN ({cats}) AND confidence >= 0.5
                        ORDER BY id""").fetchdf()
    px, py = xy_to5070(poi.x.values, poi.y.values)
    poi["gx"], poi["gy"] = np.floor(px / 200), np.floor(py / 200)
    poi["key"] = poi.name.fillna("").str.lower().str.strip()
    keep = ~poi.duplicated(["key", "cat", "gx", "gy"])
    poi, px, py = poi[keep].reset_index(drop=True), px[keep.values], py[keep.values]
    tree = shapely.STRtree(G["geoms"])
    (qq, tt), _ = tree.query_nearest(shapely.points(px, py), max_distance=150.0, return_distance=True,
                                     all_matches=False)
    pseg = np.full(len(poi), -1, np.int64)
    pseg[qq] = tt
    patt = eg.attach_buildings(G, px, py, pseg)
    qu = np.where(patt >= 0, G["pocket"][G["eu"][np.maximum(patt, 0)]], -1)
    qv = np.where(patt >= 0, G["pocket"][G["ev"][np.maximum(patt, 0)]], -1)
    ppocket = np.where(qu >= 0, qu, qv)
    poi["region"] = region
    poi["group"] = poi.cat.map(CAT2GROUP)
    poi["on_network"] = patt >= 0
    poi["pocket_bldg"] = np.where(ppocket >= 0, psize[np.maximum(ppocket, 0)], 0)
    poi["pocket_real_exit"] = np.where(ppocket >= 0, has_real_exit[np.maximum(ppocket, 0)], False)
    poi["in_pocket50"] = (poi.pocket_bldg >= MIN_POCKET_BLDG) & poi.pocket_real_exit
    ee = egress_edge.reindex(np.maximum(ppocket, 0)).values
    seg_of = np.where(np.isnan(ee), -1, G["eseg"][np.nan_to_num(ee, nan=0).astype(int)])
    poi["egress_road"] = np.where(poi.in_pocket50 & (seg_of >= 0), names.name.values[np.maximum(seg_of, 0)], None)
    t = assign_points_to_tracts(px, py, tg)
    poi["GEOID"] = np.where(t >= 0, tr.GEOID.values[np.maximum(t, 0)], None)
    poi = poi.merge(tr[["GEOID", "usfs_WHP_mean", "aiannh_name"]], on="GEOID", how="left")
    poi["whp_q"] = pd.cut(poi.usfs_WHP_mean, whp_cut, labels=["Q1", "Q2", "Q3", "Q4"], include_lowest=True)
    log(region, "facilities", len(poi), "in pockets>=50", int(poi.in_pocket50.sum()))
    return poi.drop(columns=["gx", "gy", "key"])


def main():
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    q = tm.usfs_WHP_mean.quantile([0.25, 0.5, 0.75]).values
    whp_cut = [-np.inf, q[0], q[1], q[2], np.inf]   # same pooled quartiles as 03_stratify.py (pd.qcut)
    log("pooled WHP quartile cut points", np.round(q, 2))
    allp = pd.concat([run_region(r, whp_cut) for r in REGIONS], ignore_index=True)
    allp[allp.in_pocket50].sort_values(["group", "usfs_WHP_mean"], ascending=[True, False]).to_csv(
        OUT / "pocket_facilities.csv", index=False)
    rows = []
    for (g, wq), sub in allp[allp.on_network].groupby(["group", "whp_q"], observed=True):
        k, n = int(sub.in_pocket50.sum()), len(sub)
        lo, hi = wilson(k, n)
        rows.append(dict(group=g, whp_q=str(wq), facilities_on_network=n, in_single_egress_pocket=k,
                         share=k / n if n else np.nan, ci_lo=lo, ci_hi=hi))
    for g, sub in allp[allp.on_network].groupby("group"):
        k, n = int(sub.in_pocket50.sum()), len(sub)
        lo, hi = wilson(k, n)
        rows.append(dict(group=g, whp_q="all", facilities_on_network=n, in_single_egress_pocket=k,
                         share=k / n, ci_lo=lo, ci_hi=hi))
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "pocket_facilities_summary.csv", index=False)
    pd.set_option("display.width", 200)
    print(S.round(4).to_string())


if __name__ == "__main__":
    main()
