"""What a length-based coverage metric sees in pocket tracts (permitted Overture data only).

The scorecard's transport component compares major-road length per tract (motorway / trunk / primary / secondary)
with a reference layer that was removed from the challenge and is NOT used here. As a proxy for what the
component *observes*, we compute Overture road length per tract by class, assigning each drivable segment to the
tract containing its midpoint (EPSG:5070 lengths, computed in DuckDB), and compare tracts by the share of
residents living in single-egress pockets (>= 50 buildings).
Inputs : data/reference/<region>/<region>-overture-roads.parquet, data/strata/<region>/<region>-census-tracts.parquet,
         discovery/out/tract_metrics_all.parquet
Outputs: discovery/out/coverage_proxy_tracts.csv (per tract), discovery/out/coverage_proxy_summary.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from common import OUT, REGIONS, con, log, p  # noqa: E402

MAJOR = ("motorway", "trunk", "primary", "secondary")


def region_km(region: str) -> pd.DataFrame:
    c = con()
    q = f"""
    WITH s AS (
      SELECT class,
             ST_Length(ST_Transform(geometry, 'EPSG:4326', 'EPSG:5070', always_xy := true)) AS len_m,
             CASE WHEN ST_GeometryType(geometry)::VARCHAR = 'LINESTRING' THEN ST_LineInterpolatePoint(geometry, 0.5)
                  ELSE ST_PointOnSurface(geometry) END AS mid
      FROM '{p(region, 'roads')}'),
    t AS (SELECT GEOID, geometry FROM '{p(region, 'tracts')}')
    SELECT t.GEOID, s.class, sum(s.len_m) / 1000.0 AS km
    FROM s JOIN t ON ST_Intersects(t.geometry, s.mid)
    GROUP BY ALL"""
    df = c.execute(q).fetchdf()
    w = df.pivot_table(index="GEOID", columns="class", values="km", aggfunc="sum", fill_value=0.0)
    out = pd.DataFrame(index=w.index)
    out["km_all"] = w.sum(axis=1)
    out["km_major"] = w[[k for k in MAJOR if k in w.columns]].sum(axis=1)
    out["km_below_secondary"] = out.km_all - out.km_major
    out = out.reset_index()
    out["GEOID"] = out.GEOID.astype(str)
    log(region, "tracts with roads", len(out), "km", round(out.km_all.sum()))
    return out


def main():
    km = pd.concat([region_km(r) for r in REGIONS], ignore_index=True)
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    st = pd.concat([pd.read_parquet(p(r, "strata"), columns=["GEOID", "ALAND"]) for r in REGIONS], ignore_index=True)
    st["GEOID"] = st.GEOID.astype(str)
    d = tm.merge(km, on="GEOID", how="left").merge(st, on="GEOID", how="left")
    for c in ["km_all", "km_major", "km_below_secondary"]:
        d[c] = d[c].fillna(0.0)
    d["land_km2"] = d.ALAND / 1e6
    d["pocket_share"] = d.est_pop_pocket50 / d.pop_total
    d["km_all_per_km2"] = d.km_all / d.land_km2.clip(lower=0.01)
    d["km_major_per_km2"] = d.km_major / d.land_km2.clip(lower=0.01)
    d["km_all_per_1000res"] = d.km_all / (d.pop_total / 1000)
    d["share_km_below_secondary"] = d.km_below_secondary / d.km_all.clip(lower=1e-9)
    d["pocket_bin"] = pd.cut(d.pocket_share, [-1e-9, 0.0, 0.05, 0.20, 10], labels=["0%", "0-5%", "5-20%", ">20%"])
    d[["GEOID", "region", "pop_total", "land_km2", "km_all", "km_major", "km_below_secondary", "pocket_share",
       "whp_q", "ur"]].to_csv(OUT / "coverage_proxy_tracts.csv", index=False)
    rows = []
    for scope, sub0 in [("all tracts", d), ("WHP Q4 tracts", d[d.whp_q.astype(str) == "Q4"]),
                        ("rural tracts", d[d.ur == "Rural"])]:
        for b, sub in sub0.groupby("pocket_bin", observed=True):
            rows.append(dict(scope=scope, pocket_bin=str(b), n_tracts=len(sub),
                             median_km_major=sub.km_major.median(), median_km_all=sub.km_all.median(),
                             median_km_all_per_km2=sub.km_all_per_km2.median(),
                             median_km_major_per_km2=sub.km_major_per_km2.median(),
                             median_km_all_per_1000res=sub.km_all_per_1000res.median(),
                             median_share_km_below_secondary=sub.share_km_below_secondary.median()))
        for col in ["km_major", "km_major_per_km2", "km_all_per_km2", "km_all_per_1000res"]:
            rho, pv = spearmanr(sub0.pocket_share, sub0[col])
            rows.append(dict(scope=scope, pocket_bin=f"spearman(pocket_share, {col})", n_tracts=len(sub0),
                             median_km_major=rho, median_km_all=pv))
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "coverage_proxy_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(S.round(4).to_string())
    for g in ["06103000300", "06103000201", "40145030201", "06061022014"]:
        r = d[d.GEOID == g]
        if len(r):
            r = r.iloc[0]
            print(g, f"km_all {r.km_all:.1f} km_major {r.km_major:.1f} pocket_share {r.pocket_share:.3f} "
                     f"land {r.land_km2:.1f} km2")


if __name__ == "__main__":
    main()
