"""Share of single-egress-pocket residents whose second exit exists in Overture only as a service / track /
unknown / living_street segment (i.e. it is removed by the drivable `roads` layer), by stratum.
share = 1 - est_pop_pocket50_noext / est_pop_pocket50, tract bootstrap (B = 2000, seed 20261012).
Input : discovery/out/tract_metrics_all.parquet (03_stratify.py).  Output: discovery/out/class_filter_share.csv
"""
import numpy as np
import pandas as pd

from common import OUT

SEED, B = 20261012, 2000


def main():
    d = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    rng = np.random.default_rng(SEED)
    rows = []
    for group in ["all", "region", "ur", "tribal", "whp_q", "whp_top"]:
        if group == "all":
            d["all"] = "all"
        if group == "whp_top":
            d["whp_top"] = np.where(d.whp_q.astype(str) == "Q4", "Q4", "Q1-Q3")
        for lev, sub in d.groupby(group, observed=True):
            a, b = sub.est_pop_pocket50.values, sub.est_pop_pocket50_noext.values
            idx = rng.integers(0, len(sub), size=(B, len(sub)))
            reps = 1 - b[idx].sum(1) / np.maximum(a[idx].sum(1), 1e-9)
            rows.append(dict(group=group, level=str(lev), n_tracts=len(sub), pocket_residents=a.sum(),
                             second_exit_only_low_class=a.sum() - b.sum(), share=1 - b.sum() / a.sum(),
                             ci_lo=np.percentile(reps, 2.5), ci_hi=np.percentile(reps, 97.5)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "class_filter_share.csv", index=False)
    print(out.round(4).to_string())


if __name__ == "__main__":
    main()
