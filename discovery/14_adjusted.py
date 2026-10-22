"""Headline ratio (WHP Q4 vs Q1, residents in single-egress pockets) adjusted for region and rural/urban, with a
county-block bootstrap.

The headline 3.10x (FINDINGS §1, §4.3) is crude and its CI comes from a tract bootstrap that ignores spatial
clustering. This script reports, for the share of residents in pockets (>= 50 buildings, and the stricter
"still single-egress with service/track added"):
  crude       sum(pocket residents) / sum(pop) in Q4, divided by the same in Q1
  MH          Mantel-Haenszel rate ratio Q4 vs Q1 over strata (pocket residents as events, pop as denominator)
  Poisson     log-linear model  log E[pocket residents] = log(pop) + stratum + WHP quartile  (all four quartiles,
              fitted by IRLS on stratum x quartile cells, which is exact for categorical covariates); exp(b_Q4)
Adjustment sets: none / region / rural-urban / region x rural-urban (the main one: 10 strata).
  stratum     for region x rural-urban, also the Q4/Q1 ratio inside each stratum and its share of the MH weight
CIs: county-block bootstrap. Counties (state+county FIPS) are resampled with replacement within each region,
keeping all tracts of a drawn county; B = 2000, numpy default_rng(20261014). Percentile 95% intervals.
Sensitivity (crude ratio, pocket50, pooled quartiles): the same replicates with only one region's counties
resampled and the others held fixed, to show which region drives the interval width.

Input : discovery/out/tract_metrics_all.parquet (from 03_stratify.py)
Output: discovery/out/adjusted_ratios.csv
Usage : python -I discovery/14_adjusted.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
SEED = 20261014
B = 2000
METRICS = {"pocket50": "est_pop_pocket50", "pocket50_noext": "est_pop_pocket50_noext"}
QUARTILES = {"whp_q": "pooled quartiles", "whp_q_region": "within-region quartiles"}
ADJ = {"none": [], "region": ["region"], "rural": ["ur"], "region x rural": ["region", "ur"],
       # urban-form check: residents per km2 of land (strata ALAND), terciles within region
       "region x rural x density": ["region", "ur", "dens3"]}
Q = ["Q1", "Q2", "Q3", "Q4"]


def poisson_rr(y, t, s, q, n_s):
    """exp(beta_Q4) from a Poisson GLM with offset log(t), stratum dummies s (0..n_s-1) and quartile dummies q
    (0..3, Q1 reference). Cells with t == 0 are dropped. IRLS, a few iterations suffice."""
    m = t > 0
    y, t, s, q = y[m], t[m], s[m], q[m]
    X = np.zeros((len(y), n_s + 3))
    X[np.arange(len(y)), s] = 1.0
    for k in (1, 2, 3):
        X[q == k, n_s + k - 1] = 1.0
    X = X[:, X.any(0)] if X[:, n_s + 2].any() else None
    if X is None:
        return np.nan
    beta = np.zeros(X.shape[1])
    rate0 = max(y.sum(), 1e-9) / t.sum()
    beta[: X.shape[1] - 3] = np.log(rate0)  # stratum columns first
    for _ in range(50):
        eta = X @ beta + np.log(t)
        mu = np.exp(eta)
        z = eta - np.log(t) + (y - mu) / np.maximum(mu, 1e-12)
        w = mu
        XtW = X.T * w
        new = np.linalg.solve(XtW @ X + 1e-10 * np.eye(X.shape[1]), XtW @ z)
        if np.max(np.abs(new - beta)) < 1e-10:
            beta = new
            break
        beta = new
    return float(np.exp(beta[-1]))


def mh_rr(y, t, s, q, n_s):
    """Mantel-Haenszel rate ratio Q4 (q==3) vs Q1 (q==0) over strata; strata missing either level drop out."""
    a = np.bincount(s[q == 3], y[q == 3], n_s)
    n1 = np.bincount(s[q == 3], t[q == 3], n_s)
    b = np.bincount(s[q == 0], y[q == 0], n_s)
    n0 = np.bincount(s[q == 0], t[q == 0], n_s)
    T = n1 + n0
    ok = (n1 > 0) & (n0 > 0)
    num = (a * n0 / np.where(ok, T, 1))[ok].sum()
    den = (b * n1 / np.where(ok, T, 1))[ok].sum()
    return num / den if den > 0 else np.nan


def crude_rr(y, t, q):
    r4 = y[q == 3].sum() / t[q == 3].sum()
    r1 = y[q == 0].sum() / t[q == 0].sum()
    return r4 / r1, r4, r1


def main():
    d = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    d["county"] = d.GEOID.str[:5]
    strata = Path(__file__).resolve().parents[1] / "data" / "strata"
    land = pd.concat([pd.read_parquet(p, columns=["GEOID", "ALAND"]) for p in sorted(strata.glob("*/*-strata-tract-table.parquet"))])
    d = d.merge(land.drop_duplicates("GEOID"), on="GEOID", how="left")
    dens = d.pop_total / (d.ALAND.fillna(0) / 1e6).clip(lower=0.01)
    d["dens3"] = dens.groupby(d.region).transform(lambda s: pd.qcut(s.rank(method="first"), 3, labels=["low", "mid", "high"])).astype(str)
    rng = np.random.default_rng(SEED)

    # bootstrap weights: number of times each county is drawn, resampling counties within region
    counties = d[["region", "county"]].drop_duplicates().sort_values(["region", "county"]).reset_index(drop=True)
    cidx = {c: i for i, c in enumerate(counties.county)}
    W = np.zeros((B, len(counties)))
    for r, sub in counties.groupby("region"):
        ids = sub.index.values
        draws = rng.integers(0, len(ids), size=(B, len(ids)))
        for b in range(B):
            W[b, ids] = np.bincount(draws[b], minlength=len(ids))
    tc = d.county.map(cidx).values
    print(f"{len(d)} tracts, {len(counties)} counties; county-block bootstrap B={B}, seed {SEED}")

    rows = []
    for mname, col in METRICS.items():
        for qcol, qlab in QUARTILES.items():
            qq = d[qcol].astype(str).map({k: i for i, k in enumerate(Q)}).values
            for alab, acols in ADJ.items():
                skey = d[acols].astype(str).agg("|".join, axis=1) if acols else pd.Series("all", index=d.index)
                scode, slev = pd.factorize(skey, sort=True)
                n_s = len(slev)
                # cells = county x stratum x quartile, so each bootstrap replicate is a weighted sum over counties
                cell = (tc * n_s + scode) * 4 + qq
                ncell = len(counties) * n_s * 4
                Yc = np.bincount(cell, d[col].values, ncell).reshape(len(counties), n_s * 4)
                Tc = np.bincount(cell, d.pop_total.values.astype(float), ncell).reshape(len(counties), n_s * 4)
                s_of = np.repeat(np.arange(n_s), 4)
                q_of = np.tile(np.arange(4), n_s)
                y0, t0 = Yc.sum(0), Tc.sum(0)
                cr, r4, r1 = crude_rr(y0, t0, q_of)
                est = dict(crude=cr, MH=mh_rr(y0, t0, s_of, q_of, n_s), Poisson=poisson_rr(y0, t0, s_of, q_of, n_s))
                reps = {k: np.empty(B) for k in est}
                Yb, Tb = W @ Yc, W @ Tc
                for b in range(B):
                    reps["crude"][b] = crude_rr(Yb[b], Tb[b], q_of)[0]
                    if alab != "none":
                        reps["MH"][b] = mh_rr(Yb[b], Tb[b], s_of, q_of, n_s)
                        reps["Poisson"][b] = poisson_rr(Yb[b], Tb[b], s_of, q_of, n_s)
                n_ok = int(((t0[q_of == 3] > 0) & (t0[q_of == 0] > 0)).sum())
                base = dict(metric=mname, quartiles=qlab, adjusted_for=alab, ci="county-block bootstrap",
                            n_strata=n_s, n_strata_with_q1_and_q4=n_ok, rate_q4=r4, rate_q1=r1,
                            n_tracts=len(d), n_counties=len(counties), B=B, seed=SEED)
                for k, v in est.items():
                    if (alab == "none") != (k == "crude"):
                        continue
                    lo, hi = np.nanpercentile(reps[k], [2.5, 97.5])
                    rows.append(dict(base, method=k, ratio=v, ci_lo=lo, ci_hi=hi))
                if alab == "region x rural":
                    # stratum-specific Q4/Q1 ratios and their share of the MH weight b * n1 / (n1 + n0)
                    i4, i1 = np.arange(n_s) * 4 + 3, np.arange(n_s) * 4
                    w = np.where((t0[i4] > 0) & (t0[i1] > 0), y0[i1] * t0[i4] / np.maximum(t0[i4] + t0[i1], 1), 0)
                    for s_ in range(n_s):
                        a4, a1 = i4[s_], i1[s_]
                        if t0[a4] == 0 or t0[a1] == 0:
                            continue
                        # replicates that draw no Q1 or no Q4 tract of the stratum are dropped; a zero Q1 rate
                        # makes the ratio unbounded, and ci_hi is then inf (as in 03_stratify.py)
                        with np.errstate(divide="ignore", invalid="ignore"):
                            rb = (Yb[:, a4] / Tb[:, a4]) / (Yb[:, a1] / Tb[:, a1])
                        rb = rb[(Tb[:, a4] > 0) & (Tb[:, a1] > 0)]
                        lo, hi = np.percentile(np.where(np.isinf(rb), 1e300, rb), [2.5, 97.5])
                        lo, hi = (np.inf if x >= 1e299 else x for x in (lo, hi))
                        rows.append(dict(base, method="stratum", stratum=slev[s_],
                                         ratio=(y0[a4] / t0[a4]) / (y0[a1] / t0[a1]), ci_lo=lo, ci_hi=hi,
                                         mh_weight=w[s_] / w.sum(), rate_q4=y0[a4] / t0[a4],
                                         rate_q1=y0[a1] / t0[a1], pop_q4=t0[a4], pop_q1=t0[a1]))
                if mname == "pocket50" and qcol == "whp_q" and alab == "none":
                    # which region drives the width: resample one region's counties, hold the others fixed
                    for r, sub in counties.groupby("region"):
                        W1 = np.ones_like(W)
                        W1[:, sub.index.values] = W[:, sub.index.values]
                        Y1, T1 = W1 @ Yc, W1 @ Tc
                        rb = np.array([crude_rr(Y1[b], T1[b], q_of)[0] for b in range(B)])
                        lo, hi = np.nanpercentile(rb, [2.5, 97.5])
                        rows.append(dict(base, method="crude", ratio=cr, ci_lo=lo, ci_hi=hi,
                                         ci=f"county-block bootstrap, {r} counties only"))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "adjusted_ratios.csv", index=False)
    pd.set_option("display.width", 200)
    main_rows = (R.method != "stratum") & (R.ci == "county-block bootstrap")
    print(R.loc[main_rows, ["metric", "quartiles", "adjusted_for", "method", "ratio", "ci_lo", "ci_hi", "n_strata",
                            "n_strata_with_q1_and_q4"]].round(3).to_string(index=False))
    print(R.loc[(R.method == "stratum") & (R.metric == "pocket50"),
                ["quartiles", "stratum", "ratio", "ci_lo", "ci_hi", "mh_weight", "pop_q1"]].round(3).to_string(index=False))
    print(R.loc[R.ci != "county-block bootstrap", ["ci", "ratio", "ci_lo", "ci_hi"]].round(3).to_string(index=False))
    q1 = d[d.whp_q == "Q1"].groupby("county").pop_total.sum()
    print(f"largest county share of pooled-Q1 residents: {q1.idxmax()} {q1.max() / q1.sum():.3f}")
    # population by stratum x quartile, for the write-up
    print(d.pivot_table(index=["region", "ur"], columns="whp_q", values="pop_total", aggfunc="sum",
                        observed=True).round(-3).to_string())


if __name__ == "__main__":
    main()
