"""Who lives in single-egress pockets? Age / disability / vehicle / housing-type proxies from permitted strata.

The challenge strata do NOT carry the SVI component variables (EP_AGE65, EP_DISABL, EP_NOVEH, EP_MOBILE, ...);
they carry the SVI percentile ranks of the four themes. The closest permitted proxies are therefore:
  svi_household          SVI theme 2 'household characteristics' (aged 65+, aged 17-, civilian disability,
                         single-parent households, limited English)
  svi_housing_transport  SVI theme 4 'housing type & transportation' (multi-unit, MOBILE HOMES, crowding,
                         NO VEHICLE, group quarters)
  cvi_baseline_health    Climate Vulnerability Index, baseline health domain
Rates are population-weighted shares of residents in pockets >= 50 buildings (same numerator as 03_stratify),
with a tract bootstrap (B = 2000, seed 20261012). Also the interaction with wildfire hazard (WHP Q4 vs Q1-Q3).
Inputs : discovery/out/tract_metrics_all.parquet, data/strata/<region>/<region>-strata-tract-table.parquet
Output : discovery/out/who_affected.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import OUT, REGIONS, p  # noqa: E402

SEED, B = 20261012, 2000
COLS = ["svi_household", "svi_housing_transport", "cvi_baseline_health"]


def rate_boot(sub, num, den, rng):
    x, y = sub[num].values, sub[den].values
    idx = rng.integers(0, len(sub), size=(B, len(sub)))
    reps = x[idx].sum(1) / np.maximum(y[idx].sum(1), 1e-9)
    return x.sum() / y.sum(), reps


def main():
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    st = pd.concat([pd.read_parquet(p(r, "strata"), columns=["GEOID"] + COLS) for r in REGIONS], ignore_index=True)
    st["GEOID"] = st.GEOID.astype(str)
    d = tm.merge(st[["GEOID"] + [k for k in COLS if k not in tm.columns]], on="GEOID", how="left")
    d["whp_top"] = np.where(d.whp_q.astype(str) == "Q4", "WHP Q4", "WHP Q1-Q3")
    for c in COLS:
        if c.startswith("svi"):
            d[c + "_q"] = pd.cut(d[c], [-0.001, 0.2, 0.4, 0.6, 0.8, 1.0], labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
        else:
            d[c + "_q"] = pd.qcut(d[c], 5, labels=["Q1", "Q2", "Q3", "Q4", "Q5"])
    rng = np.random.default_rng(SEED)
    rows = []
    for metric, num in [("pocket50", "est_pop_pocket50"), ("cluster_1exit", "est_pop_cluster_1exit")]:
        for c in COLS:
            q = c + "_q"
            for scope, sub0 in [("all tracts", d), ("WHP Q4 tracts", d[d.whp_top == "WHP Q4"]),
                                ("WHP Q1-Q3 tracts", d[d.whp_top == "WHP Q1-Q3"])]:
                res = {}
                for lev, sub in sub0.groupby(q, observed=True):
                    res[str(lev)] = (*rate_boot(sub, num, "pop_total", rng), len(sub))
                for lev, (rate, reps, n) in res.items():
                    rows.append(dict(metric=metric, variable=c, scope=scope, level=lev, n_tracts=n, rate=rate,
                                     ci_lo=np.percentile(reps, 2.5), ci_hi=np.percentile(reps, 97.5)))
                if "Q5" in res and "Q1" in res:
                    ra, rb = res["Q5"][1], res["Q1"][1]
                    with np.errstate(divide="ignore", invalid="ignore"):
                        rr = np.where(rb > 0, ra / np.where(rb > 0, rb, 1), np.inf)
                    rows.append(dict(metric=metric, variable=c, scope=scope, level="ratio Q5/Q1",
                                     n_tracts=res["Q5"][2] + res["Q1"][2], rate=res["Q5"][0] / res["Q1"][0]
                                     if res["Q1"][0] > 0 else np.inf, ci_lo=np.percentile(rr, 2.5),
                                     ci_hi=np.percentile(rr, 97.5)))
    # interaction: high household-vulnerability (theme 2 rank >= 0.6) x wildfire hazard
    d["hh_high"] = np.where(d.svi_household >= 0.6, "household theme >= 0.6", "household theme < 0.6")
    for metric, num in [("pocket50", "est_pop_pocket50"), ("cluster_1exit", "est_pop_cluster_1exit")]:
        for (h, w_), sub in d.groupby(["hh_high", "whp_top"]):
            rate, reps = rate_boot(sub, num, "pop_total", rng)
            rows.append(dict(metric=metric, variable="svi_household x WHP", scope=w_, level=h, n_tracts=len(sub),
                             rate=rate, ci_lo=np.percentile(reps, 2.5), ci_hi=np.percentile(reps, 97.5)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "who_affected.csv", index=False)
    show = out.copy()
    for c in ["rate", "ci_lo", "ci_hi"]:
        show[c] = np.where(show.level.str.startswith("ratio"), show[c], show[c] * 100)
    pd.set_option("display.width", 220)
    pd.set_option("display.max_rows", 400)
    print(show.round(3).to_string())


if __name__ == "__main__":
    main()
