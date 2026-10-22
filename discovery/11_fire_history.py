"""Recent wildfire history vs single-egress pockets (challenge strata only).

The challenge strata carry tract-level fire aggregates, not perimeters or fire names:
  NIFC  nifc_wildfire_burned_km2_2015_2025, nifc_wildfire_burned_pct_land_2015_2025, nifc_wildfire_last_year
  MTBS  mtbs_wildfire_burned_km2_2015_2023, mtbs_wildfire_burned_pct_land_2015_2023, mtbs_wildfire_last_year
  FOD   fod_fires_2015_2020 (ignition points), fod_largest_fire_acres
So exposure is defined per tract:
  burned15      any NIFC (2015-2025) or MTBS (2015-2023) wildfire burned area inside the tract
  burned15_1pct >= 1% of the tract's land burned by wildfire since 2015 (max of NIFC and MTBS shares)
Comparison: residents of pockets >= 50 buildings vs all other residents, inside each pooled WHP quartile (rate =
share of the group living in an exposed tract). Paired tract bootstrap (B = 2000, seed 20261012) and a
Mantel-Haenszel rate ratio pooled over the four WHP quartiles.
Inputs : discovery/out/tract_metrics_all.parquet, discovery/out/<region>_pockets.parquet, strata tract tables
Outputs: discovery/out/fire_history.csv, discovery/out/fire_history_named.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import OUT, REGIONS, p  # noqa: E402

SEED, B = 20261012, 2000
NIFC = ["GEOID", "nifc_wildfire_burned_km2_2015_2025", "nifc_wildfire_burned_pct_land_2015_2025",
        "nifc_wildfire_last_year"]
MTBS = ["GEOID", "mtbs_wildfire_burned_km2_2015_2023", "mtbs_wildfire_burned_pct_land_2015_2023",
        "mtbs_wildfire_last_year"]
FOD = ["GEOID", "fod_fires_2015_2020", "fod_largest_fire_acres", "fod_last_fire_year"]
NAMED = [("Rancho Tehama", "06103000300"), ("Lake California", "06103000201"), ("Olympic Valley", "06061022014"),
         ("Westwood Avenue (Redding)", "06089011002"), ("Icicle Road", "53007960201"),
         ("Stirrup Road (Leavenworth)", "53007960202"), ("Casassa Road (Cle Elum)", "53037975201"),
         ("Salmon la Sac Road", "53037975101"), ("Seven Bays (Davenport)", "53043960200"),
         ("Crescent Bar (Quincy)", "53025010500"), ("Whitehorn Cove (Wagoner County)", "40145030201"),
         ("Basin Road (Mannford)", "40037020706"), ("Rimrock / Lake Montezuma", "04025001605"),
         ("Bagdad", "04025001300"), ("Sargent (FM 457)", "48321730503"), ("Port O'Connor", "48057000501")]


def table(region, kind, cols):
    f = Path(p(region, "strata")).parent / f"{region}-{kind}-tract-table.parquet"
    t = pd.read_parquet(f, columns=cols)
    t["GEOID"] = t.GEOID.astype(str)
    return t


def main():
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    fire = []
    for r in REGIONS:
        t = table(r, "nifc", NIFC).merge(table(r, "mtbs", MTBS), on="GEOID", how="outer") \
            .merge(table(r, "fpa-fod", FOD), on="GEOID", how="outer")
        fire.append(t)
    fire = pd.concat(fire, ignore_index=True)
    d = tm.merge(fire, on="GEOID", how="left")
    km = d.nifc_wildfire_burned_km2_2015_2025.fillna(0) + d.mtbs_wildfire_burned_km2_2015_2023.fillna(0)
    pct = np.fmax(d.nifc_wildfire_burned_pct_land_2015_2025.fillna(0), d.mtbs_wildfire_burned_pct_land_2015_2023.fillna(0))
    d["burned15"] = km > 0
    d["burned15_1pct"] = pct >= 1.0
    d["pocket_pop"] = d.est_pop_pocket50.clip(lower=0)
    d["other_pop"] = (d.pop_total - d.pocket_pop).clip(lower=0)
    rng = np.random.default_rng(SEED)
    rows = []
    for expo in ["burned15", "burned15_1pct"]:
        reps_num, reps_den = [], []  # for the Mantel-Haenszel pooled ratio
        mh_num = mh_den = 0.0
        for q in ["Q1", "Q2", "Q3", "Q4"]:
            sub = d[d.whp_q.astype(str) == q]
            e = sub[expo].values.astype(float)
            a, o = sub.pocket_pop.values, sub.other_pop.values
            ra, ro = (a * e).sum() / a.sum(), (o * e).sum() / o.sum()
            idx = rng.integers(0, len(sub), size=(B, len(sub)))
            A, O, E = a[idx], o[idx], e[idx]
            rra = (A * E).sum(1) / A.sum(1)
            rro = (O * E).sum(1) / O.sum(1)
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = rra / rro
            # MH components: exposed pocket residents x total other / total, etc.
            n_tot = a.sum() + o.sum()
            mh_num += (a * e).sum() * o.sum() / n_tot
            mh_den += (o * e).sum() * a.sum() / n_tot
            ntot_b = A.sum(1) + O.sum(1)
            reps_num.append((A * E).sum(1) * O.sum(1) / ntot_b)
            reps_den.append((O * E).sum(1) * A.sum(1) / ntot_b)
            rows.append(dict(exposure=expo, whp_q=q, n_tracts=len(sub), pocket_residents=a.sum(),
                             pocket_residents_exposed=(a * e).sum(), other_residents=o.sum(),
                             other_residents_exposed=(o * e).sum(), rate_pocket=ra, rate_other=ro,
                             ratio=ra / ro if ro > 0 else np.inf, ratio_ci_lo=np.nanpercentile(rr, 2.5),
                             ratio_ci_hi=np.nanpercentile(rr, 97.5)))
        mh = mh_num / mh_den
        mh_reps = np.sum(reps_num, axis=0) / np.sum(reps_den, axis=0)
        rows.append(dict(exposure=expo, whp_q="MH pooled over WHP quartiles", n_tracts=len(d),
                         pocket_residents=d.pocket_pop.sum(), pocket_residents_exposed=(d.pocket_pop * d[expo]).sum(),
                         other_residents=d.other_pop.sum(), other_residents_exposed=(d.other_pop * d[expo]).sum(),
                         rate_pocket=(d.pocket_pop * d[expo]).sum() / d.pocket_pop.sum(),
                         rate_other=(d.other_pop * d[expo]).sum() / d.other_pop.sum(), ratio=mh,
                         ratio_ci_lo=np.percentile(mh_reps, 2.5), ratio_ci_hi=np.percentile(mh_reps, 97.5)))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "fire_history.csv", index=False)
    # pockets in exposed tracts
    P = pd.concat([pd.read_parquet(OUT / f"{r}_pockets.parquet") for r in REGIONS], ignore_index=True)
    P = P[(P.n_bldg >= 50) & ~P.exit_outside_region].merge(d[["GEOID", "burned15", "burned15_1pct", "whp_q"]],
                                                           on="GEOID", how="left")
    cnt = P.groupby(["whp_q", "burned15"], observed=True).agg(pockets=("n_bldg", "size"), bldg=("n_bldg", "sum"),
                                                              residents=("est_pop", "sum")).reset_index()
    cnt.to_csv(OUT / "fire_history_pocket_counts.csv", index=False)
    # named communities: tract fire record
    nm = pd.DataFrame(NAMED, columns=["label", "GEOID"]).merge(
        d[["GEOID", "whp_q", "nifc_wildfire_burned_km2_2015_2025", "nifc_wildfire_burned_pct_land_2015_2025",
           "nifc_wildfire_last_year", "mtbs_wildfire_burned_km2_2015_2023", "mtbs_wildfire_last_year",
           "fod_fires_2015_2020", "fod_largest_fire_acres", "burned15", "burned15_1pct"]], on="GEOID", how="left")
    nm.to_csv(OUT / "fire_history_named.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(R.round(4).to_string())
    print(cnt.round(0).to_string())
    print(nm.round(2).to_string())


if __name__ == "__main__":
    main()
