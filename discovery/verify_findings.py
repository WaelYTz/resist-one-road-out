"""Recompute every headline number of the write-up from discovery/out and compare it with the published value.

Each check reads an output table written by the numbered scripts (01-18), recomputes the number and compares it
at the precision it is printed with in FINDINGS.md (§ given per check). Prints PASS / FAIL per check and a total;
exits 1 if any check fails. Runs in a few seconds and needs only pandas.
Usage: python -I discovery/verify_findings.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
results = []


def check(name, value, expected, decimals=0, where=""):
    ok = round(float(value), decimals) == round(float(expected), decimals)
    results.append(ok)
    shown = f"{value:,.{decimals}f}"
    print(f"{'PASS' if ok else 'FAIL'}  {name:62s} {shown:>14s}  published {expected:,.{decimals}f}  ({where.replace('§', 'sec. ')})")


def main():
    # §1, §4.2: pockets, buildings and residents behind one segment
    reg = pd.read_csv(OUT / "pocket_summary_by_region.csv")
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    check("single-egress pockets", reg.pockets.sum(), 13775, 0, "§4.2")
    check("buildings in pockets", reg.bldg.sum(), 1562885, 0, "§4.2")
    check("residents in pockets", reg.residents.sum(), 2884039, 0, "§4.2")
    check("share of residents in pockets, %", 100 * reg.residents.sum() / tm.pop_total.sum(), 7.13, 2, "§4.2")

    # §1, §4.3: wildfire-hazard gradient (residents in pockets >= 50 buildings, pooled WHP quartiles)
    sr = pd.read_csv(OUT / "strata_rates.csv")
    q = sr[(sr.metric == "pocket50") & (sr.group == "whp_q")].set_index("level")
    check("WHP Q4: residents in pockets, %", 100 * q.loc["Q4", "rate"], 10.15, 2, "§1")
    check("WHP Q1: residents in pockets, %", 100 * q.loc["Q1", "rate"], 3.28, 2, "§1")
    check("WHP Q4 / Q1 ratio (crude)", q.loc["Q4", "rate"] / q.loc["Q1", "rate"], 3.10, 2, "§1, §4.4")

    # §4.4: adjusted Mantel-Haenszel ratios with county-block bootstrap CIs
    ar = pd.read_csv(OUT / "adjusted_ratios.csv")
    sel = (ar.metric == "pocket50") & (ar.method == "MH") & (ar.quartiles == "pooled quartiles") \
        & (ar.ci == "county-block bootstrap")
    for adj, (r, lo, hi) in {"region x rural": (3.31, 2.29, 4.87),
                             "region x rural x density": (1.71, 1.23, 2.52),
                             "county": (5.25, 3.82, 7.27)}.items():
        row = ar[sel & (ar.adjusted_for == adj)].iloc[0]
        check(f"MH ratio adjusted for {adj}", row.ratio, r, 2, "§4.4")
        check(f"  lower CI bound ({adj})", row.ci_lo, lo, 2, "§4.4")
        check(f"  upper CI bound ({adj})", row.ci_hi, hi, 2, "§4.4")

    # §4.10: external validation (random WHP-Q4 sample; resolved pockets only)
    ov = pd.read_csv(OUT / "official_validation.csv")
    rq = ov[ov["sample"] == "random_q4"]
    osm = rq[rq.osm_class.isin(["a", "b", "c"])]
    check("OSM: random high-hazard pockets resolved", len(osm), 33, 0, "§4.10")
    check("OSM: no public second road in OSM either", osm.osm_class.isin(["a", "c"]).sum(), 29, 0, "§4.10")
    check("OSM: Overture gap (second road only in OSM)", (osm.osm_class == "b").sum(), 4, 0, "§4.10")
    off = rq[rq.off_class.isin(["A_off", "B_off", "C_off"])]
    check("HPMS: pockets resolved", len(off), 33, 0, "§4.10")
    check("HPMS: no officially public second exit", off.off_class.isin(["A_off", "C_off"]).sum(), 31, 0, "§4.10")

    # §4.11: who is affected (ACS 2020-2024)
    ac = pd.read_csv(OUT / "acs_evacuation.csv")
    a65 = ac[(ac.scope == "all tracts") & (ac.variable == "age65plus") & (ac.measure == "ratio")
             & (ac.group == "pocket/all")].iloc[0]
    check("65+ share, pocket residents / all residents", a65.value, 1.12, 2, "§4.11")

    # §4.13: recent wildfire history (MH pooled over WHP quartiles)
    fh = pd.read_csv(OUT / "fire_history.csv")
    f = fh[(fh.exposure == "burned15") & (fh.whp_q == "MH pooled over WHP quartiles")].iloc[0]
    check("burned since 2015: pocket vs same-hazard residents", f.ratio, 1.15, 2, "§4.13")

    # §4.13b: CAL FIRE damage inspections (external), destroyed share pocket vs other, same incident and tract
    dn = pd.read_csv(OUT / "dins_pockets.csv")
    r = dn[(dn.subset == "all structures") & (dn.measure == "mh_rr_incident_tract")].iloc[0]
    check("DINS: destroyed, pocket vs other, same incident x tract (RR)", r.value, 0.89, 2, "§4.13b")

    # §5b: action list
    al = pd.read_csv(OUT / "action_list_full.csv", dtype={"county_fips": str})
    check("action list: high-hazard single-exit pockets", len(al), 5095, 0, "§5b")
    check("action list: residents (thousands)", al.est_pop.sum() / 1000, 825, 0, "§5b")
    check("action list: counties", al.county_fips.nunique(), 195, 0, "§5b")

    n, k = len(results), sum(results)
    print(f"\n{'PASS' if k == n else 'FAIL'} {k}/{n}")
    sys.exit(0 if k == n else 1)


if __name__ == "__main__":
    main()
