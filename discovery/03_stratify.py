"""Stratify the topology metrics by region / rural-urban / tribal / SVI / wildfire hazard, with tract bootstrap CIs.

Inputs : discovery/out/<region>_tract_islands.parquet, <region>_tract_egress.parquet, <region>_pockets.parquet,
         <region>_islands.parquet, <region>_clusters.parquet  (from 01_islands.py / 02_egress.py)
         data/strata/<region>/<region>-strata-tract-table.parquet
Outputs: discovery/out/tract_metrics_all.parquet (+ .csv)   per-tract metrics + strata, all regions
         discovery/out/strata_rates.csv                     population-weighted rates by stratum (95% CI)
         discovery/out/strata_ratios.csv                    rate ratios between strata (95% CI)
         discovery/out/within_region_ratios.csv             same comparisons inside each region
         discovery/out/top_pockets.csv, top_pockets_wildfire.csv, top_single_exit_clusters.csv, top_islands.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from common import OUT, REGIONS, p

SEED = 20261012
B = 2000

METRICS = {
    # name: (numerator column, denominator column, description)
    "pocket50": ("est_pop_pocket50", "pop_total", "residents in single-egress pockets (>=50 bldgs)"),
    "pocket50_noext": ("est_pop_pocket50_noext", "pop_total", "same, still single-egress with service/track added"),
    "pocket100": ("est_pop_pocket100", "pop_total", "residents in single-egress pockets (>=100 bldgs)"),
    "cluster_1exit": ("est_pop_cluster_1exit", "pop_total", "residents of settlement clusters with 1 exit"),
    "cluster_le2exit": ("est_pop_cluster_le2exit", "pop_total", "residents of settlement clusters with <=2 exits"),
    "cluster_unmapped": ("est_pop_cluster_unmapped", "pop_total", "residents of dense clusters whose streets are "
                         "mostly missing from the drivable layer"),
    "island_all": ("est_pop_island", "pop_total", "residents attached to drivable routing islands"),
    "island_strict": ("est_pop_island_strict", "pop_total", "same, not reconnected by service/track"),
    "far150": ("est_pop_far150", "pop_total", "residents >150 m from any drivable segment"),
    "bldg_pocket50": ("n_bldg_pocket50", "n_bldg", "buildings in single-egress pockets (>=50 bldgs)"),
}


def load_all() -> pd.DataFrame:
    parts = []
    for r in REGIONS:
        a = pd.read_parquet(OUT / f"{r}_tract_islands.parquet")
        b = pd.read_parquet(OUT / f"{r}_tract_egress.parquet").drop(columns=["pop_total", "region", "n_bldg"])
        s = pd.read_parquet(p(r, "strata"), columns=["GEOID", "ur_class", "pct_urban", "tribal_any", "tribal_legal",
                                                     "aiannh_name", "aiannh_kind", "svi_overall", "svi_housing_transport",
                                                     "usfs_WHP_mean"])
        s["GEOID"] = s.GEOID.astype(str)
        parts.append(a.merge(b, on="GEOID").merge(s, on="GEOID", how="left"))
    d = pd.concat(parts, ignore_index=True)
    d = d[d.pop_total > 0].copy()
    d["svi_q"] = pd.cut(d.svi_overall, [-0.001, 0.2, 0.4, 0.6, 0.8, 1.0], labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    d["svi_ht_q"] = pd.cut(d.svi_housing_transport, [-0.001, 0.2, 0.4, 0.6, 0.8, 1.0], labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    d["whp_q"] = pd.qcut(d.usfs_WHP_mean, 4, labels=["Q1", "Q2", "Q3", "Q4"])
    d["whp_q_region"] = d.groupby("region").usfs_WHP_mean.transform(lambda s: pd.qcut(s, 4, labels=["Q1", "Q2", "Q3", "Q4"]))
    d["whp_top"] = np.where(d.whp_q == "Q4", "Q4", "Q1-Q3")
    d["tribal"] = np.where(d.tribal_any.fillna(False).astype(bool), "tribal", "non-tribal")
    d["ur"] = d.ur_class.fillna("NA")
    return d


def boot_rates(d, num, den, group, rng):
    """Population-weighted rate sum(num)/sum(den) per group value, with B tract-bootstrap replicates."""
    out = {}
    for g, sub in d.groupby(group, observed=True):
        x, y = sub[num].values, sub[den].values
        n = len(sub)
        idx = rng.integers(0, n, size=(B, n))
        reps = x[idx].sum(1) / np.maximum(y[idx].sum(1), 1e-9)
        out[g] = (x.sum() / max(y.sum(), 1e-9), reps, n, y.sum(), x.sum())
    return out


def ci(a):
    return np.percentile(a, 2.5), np.percentile(a, 97.5)


def compare(res, a, b, extra):
    """Rate ratio and rate difference (a - b) with percentile bootstrap CIs. When the reference rate is zero in
    some replicates the ratio is unbounded; ratio_ci_hi is then reported as inf."""
    ra, rb = res[a][1], res[b][1]
    with np.errstate(divide="ignore", invalid="ignore"):
        rr = np.where(rb > 0, ra / np.where(rb > 0, rb, 1), np.inf)
    lo, hi = np.percentile(rr, [2.5, 97.5])
    dlo, dhi = ci(ra - rb)
    out = dict(extra)
    out.update(rate_a=res[a][0], rate_b=res[b][0], ratio=res[a][0] / res[b][0] if res[b][0] > 0 else np.inf,
               ratio_ci_lo=lo, ratio_ci_hi=hi, diff=res[a][0] - res[b][0], diff_ci_lo=dlo, diff_ci_hi=dhi,
               n_a=res[a][2], n_b=res[b][2], pop_a=res[a][3], pop_b=res[b][3], affected_a=res[a][4])
    return out


COMPARISONS = [  # (group column, numerator level, denominator level, label)
    ("ur", "Rural", "Urban", "rural vs urban"),
    ("tribal", "tribal", "non-tribal", "tribal vs non-tribal"),
    ("svi_q", "Q5", "Q1", "SVI Q5 (most vulnerable) vs Q1"),
    ("svi_ht_q", "Q5", "Q1", "SVI housing-type/transport theme Q5 vs Q1"),
    ("whp_q", "Q4", "Q1", "wildfire hazard (WHP) Q4 vs Q1, pooled quartiles"),
    ("whp_q_region", "Q4", "Q1", "wildfire hazard Q4 vs Q1, within-region quartiles"),
    ("whp_top", "Q4", "Q1-Q3", "wildfire hazard Q4 vs Q1-Q3, pooled quartiles"),
]


def main():
    rng = np.random.default_rng(SEED)
    d = load_all()
    keep = ["GEOID", "region", "pop_total", "n_bldg", "ur", "pct_urban", "tribal", "tribal_legal", "aiannh_name",
            "svi_overall", "svi_q", "svi_housing_transport", "svi_ht_q", "usfs_WHP_mean", "whp_q", "whp_q_region", "tract_locality"]
    cols = keep + [c for c in d.columns if c.startswith(("n_bldg_", "est_pop_", "share_"))]
    d[cols].to_parquet(OUT / "tract_metrics_all.parquet", index=False)
    d[cols].to_csv(OUT / "tract_metrics_all.csv", index=False)

    rates, ratios = [], []
    for mname, (num, den, desc) in METRICS.items():
        for group in ["region", "ur", "tribal", "svi_q", "svi_ht_q", "whp_q", "whp_q_region", "whp_top"]:
            res = boot_rates(d, num, den, group, rng)
            for g, (rate, reps, n, popd, popn) in res.items():
                lo, hi = ci(reps)
                rates.append(dict(metric=mname, desc=desc, group=group, level=str(g), n_tracts=n, denom=popd,
                                  numer=popn, rate=rate, ci_lo=lo, ci_hi=hi))
            for gcol, a, b, lab in COMPARISONS:
                if gcol != group or a not in res or b not in res:
                    continue
                ratios.append(compare(res, a, b, dict(metric=mname, comparison=lab)))
        # overall
        x, y = d[num].values, d[den].values
        idx = rng.integers(0, len(d), size=(B, len(d)))
        lo, hi = ci(x[idx].sum(1) / y[idx].sum(1))
        rates.append(dict(metric=mname, desc=desc, group="all", level="all", n_tracts=len(d), denom=y.sum(),
                          numer=x.sum(), rate=x.sum() / y.sum(), ci_lo=lo, ci_hi=hi))
    pd.DataFrame(rates).to_csv(OUT / "strata_rates.csv", index=False)
    pd.DataFrame(ratios).to_csv(OUT / "strata_ratios.csv", index=False)

    # within-region comparisons (removes between-region confounding)
    wr = []
    for r, sub in d.groupby("region"):
        for mname in ["pocket50", "pocket50_noext", "cluster_1exit", "island_all"]:
            num, den, _ = METRICS[mname]
            for gcol, a, b, lab in COMPARISONS:
                if gcol == "whp_q":
                    continue
                res = boot_rates(sub, num, den, gcol, rng)
                if a not in res or b not in res or res[a][2] < 5 or res[b][2] < 5:
                    continue
                wr.append(compare(res, a, b, dict(region=r, metric=mname, comparison=lab)))
    for u_, sub in d.groupby("ur"):
        if u_ not in ("Rural", "Urban"):
            continue
        for mname in ["pocket50", "pocket50_noext", "cluster_1exit", "island_all"]:
            num, den, _ = METRICS[mname]
            for gcol, a, b, lab in [("whp_q", "Q4", "Q1", "wildfire hazard Q4 vs Q1 (pooled quartiles)"),
                                    ("whp_top", "Q4", "Q1-Q3", "wildfire hazard Q4 vs Q1-Q3 (pooled quartiles)")]:
                res = boot_rates(sub, num, den, gcol, rng)
                if a in res and b in res:
                    wr.append(compare(res, a, b, dict(region=f"ALL ({u_} tracts only)", metric=mname, comparison=lab)))
    pd.DataFrame(wr).to_csv(OUT / "within_region_ratios.csv", index=False)

    # ---------------- named examples ----------------
    whp = d[["GEOID", "whp_q", "whp_q_region", "svi_overall", "ur", "tract_locality"]]
    P = pd.concat([pd.read_parquet(OUT / f"{r}_pockets.parquet") for r in REGIONS], ignore_index=True)
    P = P[(P.n_bldg >= 50) & ~P.exit_outside_region].merge(whp, on="GEOID", how="left", suffixes=("", "_t"))
    pc = ["region", "GEOID", "locality", "tract_locality", "aiannh_name", "n_bldg", "est_pop", "road_km", "egress_road",
          "egress_class", "ext_second_exit", "near_miss", "dangle_gap_m", "min_dist_outline_km", "usfs_WHP_mean",
          "whp_q", "svi_overall", "ur", "lat", "lon"]
    P.sort_values("n_bldg", ascending=False)[pc].head(60).to_csv(OUT / "top_pockets.csv", index=False)
    mix = P.groupby("egress_class").agg(pockets=("n_bldg", "size"), residents=("est_pop", "sum"))
    mix["share_residents"] = mix.residents / mix.residents.sum()
    mix.sort_values("residents", ascending=False).to_csv(OUT / "pocket_egress_class_mix.csv")
    P.groupby(["region"]).agg(pockets=("n_bldg", "size"), bldg=("n_bldg", "sum"), residents=("est_pop", "sum"),
                              share_ext_second_exit=("ext_second_exit", "mean"), share_near_miss=("near_miss", "mean"),
                              share_within_1km_outline=("min_dist_outline_km", lambda s: (s < 1).mean())
                              ).to_csv(OUT / "pocket_summary_by_region.csv")
    W = P[(P.whp_q == "Q4") & ~P.ext_second_exit].sort_values("n_bldg", ascending=False)
    W[pc].head(60).to_csv(OUT / "top_pockets_wildfire.csv", index=False)
    C = pd.concat([pd.read_parquet(OUT / f"{r}_clusters.parquet") for r in REGIONS], ignore_index=True)
    C = C.merge(whp[["GEOID", "whp_q", "svi_overall"]].rename(columns={"svi_overall": "svi"}), on="GEOID", how="left")
    cc = ["region", "GEOID", "locality", "aiannh_name", "n_bldg", "est_pop", "exits", "exits_ext", "share_in_pocket",
          "usfs_WHP_mean", "whp_q", "svi", "lat", "lon"]
    cc = cc + ["share_attached"]
    C[(C.exits == 1) & C.exits_valid].sort_values("n_bldg", ascending=False)[cc].head(60).to_csv(
        OUT / "top_single_exit_clusters.csv", index=False)
    C[~C.exits_valid].sort_values("n_bldg", ascending=False)[cc].head(60).to_csv(
        OUT / "top_unmapped_clusters.csv", index=False)
    I = pd.concat([pd.read_parquet(OUT / f"{r}_islands.parquet") for r in REGIONS], ignore_index=True)
    I.sort_values("n_bldg", ascending=False).head(80).to_csv(OUT / "top_islands.csv", index=False)

    # ---------------- console summary ----------------
    R = pd.DataFrame(rates)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(R[R.metric.isin(["pocket50", "pocket50_noext", "cluster_1exit", "island_all", "island_strict"])]
          .round(5).to_string())
    print(pd.DataFrame(ratios).round(4).to_string())
    print(pd.DataFrame(wr).round(4).to_string())
    print(mix.sort_values("residents", ascending=False).round(3).to_string())


if __name__ == "__main__":
    main()
