"""Evacuation vulnerability of pocket residents: households without a vehicle, limited-English households and
residents aged 65+, from the American Community Survey (ACS 2024 5-year, tract level). External data.

The challenge strata carry SVI theme ranks only (08_who_affected.py), not the component shares. This script reads
the component tables themselves and asks: do residents of single-egress pockets live in tracts where more (or
fewer) households lack a car, speak limited English, or where more residents are 65+?

Source: U.S. Census Bureau, American Community Survey 2020-2024 5-year estimates (release id acs2024_5yr), tables
  B08201  household size by vehicles available            no_vehicle      = B08201_002 / B08201_001   (households)
  C16002  household language by limited-English status    limited_english = (C16002_004 + _007 + _010 + _013)
                                                                            / C16002_001                (households)
  B01001  sex by age                                      age65plus       = (B01001_020..025 + B01001_044..049)
                                                                            / B01001_001                (persons)
  Served by the Census Reporter API (https://api.censusreporter.org/1.0/data/show/acs2024_5yr), which returns the
  Bureau's estimates and 90% margins of error. The Bureau's own API (api.census.gov) now redirects requests
  without a key to missing_key.html, so it is not used. Only the 9,767 populated tracts of tract_metrics_all are
  requested (explicit 14000US<GEOID> ids, sorted, BATCH per request), not whole states. Each gzip response is
  cached in discovery/external/acs/ with a .meta.json (URL, UTC retrieval time, HTTP status, bytes on the wire,
  release). ACS tables are U.S. Government works (public domain, 17 U.S.C. 105).

Ecological estimate (tract level): the share for pocket residents is sum_t w_t r_t / sum_t w_t, where r_t is the
  tract's ACS rate and w_t its estimated pocket residents (est_pop_pocket50 from 03_stratify.py). Groups:
    pocket       w = est_pop_pocket50
    action list  w = est_pop of the tract's pockets in out/action_list_full.csv (15_full_action_list.py)
    all          w = pop_total
    other        w = pop_total - est_pop_pocket50
    expected     the all-resident rate of each region x rural stratum, weighted by the group's (pocket or action
                 list) residents per stratum: what the group's share would be if it matched its stratum
                 neighbours, so the ratio to it removes where pockets are (region, rural class)
  Household rates are weighted by residents, so they read "the share of households without a car in the tract of
  the average pocket resident". A tract's rate is applied to all of its residents: differences between the pocket
  and the rest of its tract are invisible. Tracts whose ACS denominator is 0 drop out of that variable.
Scopes: all tracts; WHP Q4 / Q1-Q3 (pooled quartiles); rural / urban (ur, as in 03_stratify.py); each region.
CIs: county-block bootstrap (counties resampled with replacement within region, all tracts of a drawn county kept,
  as in 14_adjusted.py), B = 2000, numpy default_rng(20261018); percentile 95% intervals (ci_lo, ci_hi).
  Sensitivity (ci_lo_acs, ci_hi_acs): the same replicates with each tract's rate also redrawn from its ACS sampling
  error, normal with SE = MOE/1.645 (Census proportion formula; numerator MOE of a sum = root of summed squared
  MOEs, which overstates it), clipped to [0, 1].
  Robustness (scope "all tracts", groups "... (pooled ACS counts)"): ACS counts pooled instead of rates averaged
  (pocket = sum_t s_t num_t / sum_t s_t den_t with s_t = est_pop_pocket50 / pop_total), so household rates are
  weighted by households. Same county-block CIs, no ACS-noise variant.
  est_persons (age65plus only): the group's share x its estimated residents; the CI is the share's CI x the same
  fixed count (resampling counties would also resample the population total).

Inputs : discovery/out/tract_metrics_all.parquet (03_stratify.py), discovery/out/action_list_full.csv
         (15_full_action_list.py), discovery/external/acs/*.json.gz (downloaded here)
Output : discovery/out/acs_evacuation.csv  (long format: scope, variable, measure, group, value, CIs, n_tracts)
Usage  : python -I discovery/17_acs_evacuation.py            downloads missing batches (~10 MB), then analyses
         python -I discovery/17_acs_evacuation.py --offline  cache only; stops if a batch is missing
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

DISC = Path(__file__).resolve().parent
OUT = DISC / "out"
RAW = DISC / "external" / "acs"  # downloads, never imported
RELEASE = "acs2024_5yr"
API = f"https://api.censusreporter.org/1.0/data/show/{RELEASE}"
TABLES = ["B08201", "C16002", "B01001"]
USER_AGENT = "bias-bounty-mapping-equity-discovery/1.0 (research script; ACS tract tables, small batches)"
BATCH = 200
SEED = 20261018
B = 2000
CHUNK = 250  # bootstrap replicates per block in the ACS-noise sensitivity
AGE65 = [f"B01001{i:03d}" for i in [*range(20, 26), *range(44, 50)]]
VARS = {  # name: (numerator columns, denominator column, unit)
    "no_vehicle": (["B08201002"], "B08201001", "households"),
    "limited_english": (["C16002004", "C16002007", "C16002010", "C16002013"], "C16002001", "households"),
    "age65plus": (AGE65, "B01001001", "persons"),
}
NEEDED = sorted({c for num, den, _ in VARS.values() for c in [*num, den]})
GROUP_LABEL = {"pocket": "pocket residents", "action": "action-list residents", "all": "all residents",
               "other": "residents outside pockets"}


# ----------------------------------------------------------------------------------------------- download
def batches(geoids):
    g = sorted(geoids)
    return [g[i:i + BATCH] for i in range(0, len(g), BATCH)]


def cache_path(i, chunk):
    return RAW / f"{RELEASE}_b{i:03d}_{chunk[0]}_{chunk[-1]}.json.gz"


def fetch(i, chunk) -> int:
    """Download one batch; store the gzip body and a .meta.json beside it. Returns bytes on the wire."""
    path = cache_path(i, chunk)
    q = urllib.parse.urlencode({"table_ids": ",".join(TABLES), "geo_ids": ",".join("14000US" + g for g in chunk)},
                               safe=",")
    url = f"{API}?{q}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                body, enc, status = r.read(), r.headers.get("Content-Encoding", ""), r.status
            break
        except urllib.error.HTTPError as e:
            msg = e.read()[:500]
            if e.code < 500 or attempt == 3:
                raise SystemExit(f"batch {i}: HTTP {e.code}: {msg!r}")
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == 3:
                raise SystemExit(f"batch {i}: {e}")
        time.sleep(10 * (attempt + 1))
    raw = body if enc == "gzip" else gzip.compress(body)
    js = json.loads(gzip.decompress(raw))
    assert js["release"]["id"] == RELEASE, js.get("release")
    path.write_bytes(raw)
    meta = dict(url=url, retrieved_utc=dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                http_status=status, content_encoding=enc, bytes_wire=len(body), bytes_json=len(gzip.decompress(raw)),
                release=js["release"], tables=TABLES, n_geo_requested=len(chunk), n_geo_returned=len(js["data"]))
    path.with_suffix("").with_suffix(".meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print(f"  batch {i:03d}: {len(chunk)} tracts, {len(body) / 1e3:.0f} kB on the wire", flush=True)
    return len(body)


def load_acs(geoids, offline: bool) -> pd.DataFrame:
    RAW.mkdir(parents=True, exist_ok=True)
    chunks = batches(geoids)
    missing = [(i, c) for i, c in enumerate(chunks) if not cache_path(i, c).exists()]
    if missing and offline:
        raise SystemExit(f"--offline: {len(missing)} of {len(chunks)} batches missing in {RAW}, "
                         f"e.g. {cache_path(*missing[0]).name}")
    wire = 0
    for k, (i, c) in enumerate(missing):
        wire += fetch(i, c)
        if k + 1 < len(missing):
            time.sleep(1.0)
    rows, metas = [], []
    for i, c in enumerate(chunks):
        p = cache_path(i, c)
        js = json.loads(gzip.decompress(p.read_bytes()))
        assert js["release"]["id"] == RELEASE
        metas.append(json.loads(p.with_suffix("").with_suffix(".meta.json").read_text(encoding="utf-8")))
        for gid, tabs in js["data"].items():
            rec = {"GEOID": gid.removeprefix("14000US")}
            for col in NEEDED:
                t = tabs.get(col[:6], {})
                rec[col] = t.get("estimate", {}).get(col)
                rec[col + "_moe"] = t.get("error", {}).get(col)
            rows.append(rec)
    m = pd.DataFrame(metas)
    print(f"ACS cache: {len(chunks)} batches, {m.bytes_wire.sum() / 1e6:.2f} MB on the wire "
          f"({m.bytes_json.sum() / 1e6:.1f} MB JSON), retrieved {m.retrieved_utc.min()} to {m.retrieved_utc.max()}; "
          f"downloaded now: {wire / 1e6:.2f} MB; release {metas[0]['release']}")
    a = pd.DataFrame(rows).drop_duplicates("GEOID")
    return a.apply(lambda s: pd.to_numeric(s, errors="coerce") if s.name != "GEOID" else s)


def tract_rates(a: pd.DataFrame) -> pd.DataFrame:
    """Rate and sampling SE per tract and variable (Census proportion MOE formula; ratio formula if negative)."""
    out = a[["GEOID"]].copy()
    for v, (num, den, _) in VARS.items():
        n = a[num].sum(axis=1, min_count=len(num))
        n_moe = np.sqrt((a[[c + "_moe" for c in num]].fillna(0) ** 2).sum(axis=1))
        d, d_moe = a[den], a[den + "_moe"].fillna(0)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(d > 0, n / d, np.nan)
            rad = n_moe ** 2 - p ** 2 * d_moe ** 2
            rad = np.where(rad < 0, n_moe ** 2 + p ** 2 * d_moe ** 2, rad)
            se = np.where(d > 0, np.sqrt(rad) / d / 1.645, np.nan)
        out[v], out[v + "_se"], out[v + "_den"] = np.clip(p, 0, 1), se, d
    return out


# ----------------------------------------------------------------------------------------------- analysis
def county_weights(d, rng):
    counties = d[["region", "county"]].drop_duplicates().sort_values(["region", "county"]).reset_index(drop=True)
    W = np.zeros((B, len(counties)))
    for _, sub in counties.groupby("region"):
        ids = sub.index.values
        draws = rng.integers(0, len(ids), size=(B, len(ids)))
        for b in range(B):
            W[b, ids] = np.bincount(draws[b], minlength=len(ids))
    tc = d.county.map({c: i for i, c in enumerate(counties.county)}).values
    return W, tc, len(counties)


def pct(x):
    return np.nanpercentile(x, [2.5, 97.5])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true", help="read the cached responses only")
    args = ap.parse_args()

    d = pd.read_parquet(OUT / "tract_metrics_all.parquet",
                        columns=["GEOID", "region", "pop_total", "est_pop_pocket50", "whp_q", "ur"])
    al = pd.read_csv(OUT / "action_list_full.csv", dtype={"GEOID": str})
    acs = load_acs(d.GEOID.tolist(), args.offline)
    r = tract_rates(acs)
    d = d.merge(r, on="GEOID", how="left", validate="1:1")
    d = d.merge(al.groupby("GEOID").est_pop.sum().rename("est_pop_action"), left_on="GEOID", right_index=True,
                how="left")
    d["est_pop_action"] = d.est_pop_action.fillna(0.0)
    assert np.isclose(d.est_pop_action.sum(), al.est_pop.sum()), "action-list tract missing from tract_metrics"
    d["county"] = d.GEOID.str[:5]
    d["stratum"] = d.region + "|" + d.ur
    s_code, s_lev = pd.factorize(d.stratum, sort=True)
    n_s = len(s_lev)
    n_miss = int(d.no_vehicle_den.isna().sum())
    print(f"{len(d)} tracts; matched in ACS: {len(d) - n_miss}; zero households (rate undefined): "
          f"{int((d.no_vehicle_den == 0).sum())} tracts, {d.loc[d.no_vehicle_den == 0, 'est_pop_pocket50'].sum():.0f}"
          f" pocket residents; zero population: {int((d.age65plus_den == 0).sum())} tracts; ACS population / "
          f"strata pop_total = {d.age65plus_den.sum() / d.pop_total.sum():.4f}")

    wts = {"pocket": d.est_pop_pocket50.values, "action": d.est_pop_action.values, "all": d.pop_total.values * 1.0,
           "other": np.clip(d.pop_total.values - d.est_pop_pocket50.values, 0, None)}
    scopes = {"all tracts": np.ones(len(d), bool), "WHP Q4": (d.whp_q.astype(str) == "Q4").values,
              "WHP Q1-Q3": (d.whp_q.astype(str) != "Q4").values, "rural": (d.ur == "Rural").values,
              "urban": (d.ur == "Urban").values}
    scopes.update({f"region {g}": (d.region == g).values for g in sorted(d.region.unique())})

    rng = np.random.default_rng(SEED)
    W, tc, n_c = county_weights(d, rng)
    W1 = np.vstack([np.ones((1, n_c)), W])  # row 0 = point estimate
    cell = tc * n_s + s_code
    S1 = np.zeros((len(d), n_s))
    S1[np.arange(len(d)), s_code] = 1.0
    # ACS sampling noise: one standard-normal draw per replicate and tract per variable, shared by all groups
    Z = {v: rng.standard_normal((B, len(d))).astype(np.float32) for v in VARS}
    print(f"county-block bootstrap: {n_c} counties, B={B}, seed {SEED}")

    rows = []
    for v, (_, _, unit) in VARS.items():
        rate, se = d[v].values, d[v + "_se"].values
        ok = ~np.isnan(rate)
        r0, se0 = np.where(ok, rate, 0.0), np.where(ok, se, 0.0)
        for scope, msk in scopes.items():
            use = msk & ok
            N, D, Nn = {}, {}, {}
            for g, w in wts.items():
                ww = np.where(use, w, 0.0)
                Dc = np.bincount(cell, ww, n_c * n_s).reshape(n_c, n_s)
                Nc = np.bincount(cell, ww * r0, n_c * n_s).reshape(n_c, n_s)
                D[g], N[g] = W1 @ Dc, W1 @ Nc  # (B+1) x strata
                # with ACS noise: rates redrawn per replicate, so sum tract by tract in blocks of replicates
                nn = np.empty((B, n_s))
                for b0 in range(0, B, CHUNK):
                    b1 = min(b0 + CHUNK, B)
                    rr = np.clip(r0 + Z[v][b0:b1] * se0, 0, 1)
                    nn[b0:b1] = (W[b0:b1][:, tc] * ww * rr) @ S1
                Nn[g] = np.vstack([N[g][:1], nn])
            n_tr = {g: int((use & (w > 0)).sum()) for g, w in wts.items()}
            if D["pocket"][0].sum() == 0:
                continue

            def share(g, NN):
                return NN[g].sum(1) / D[g].sum(1)

            def expected(g, NN):  # all-resident stratum rates, weighted by group g's residents per stratum
                with np.errstate(divide="ignore", invalid="ignore"):
                    ra = np.where(D["all"] > 0, NN["all"] / D["all"], 0.0)
                return (D[g] * ra).sum(1) / D[g].sum(1)

            groups = ["pocket", "all", "other"] + (["action"] if D["action"][0].sum() > 0 else [])
            est = {}
            for NN, tag in [(N, ""), (Nn, "_acs")]:
                with np.errstate(divide="ignore", invalid="ignore"):
                    sh = {g: share(g, NN) for g in groups}
                    for g in [g for g in ("pocket", "action") if g in groups]:
                        e = expected(g, NN)
                        sh[f"expected|{g}"] = e
                        sh[f"ratio|{g}/all"] = sh[g] / sh["all"]
                        sh[f"diff_pp|{g}-all"] = 100 * (sh[g] - sh["all"])
                        sh[f"ratio|{g}/expected"] = sh[g] / e
                        if g == "pocket":
                            sh["ratio|pocket/other"] = sh["pocket"] / sh["other"]
                        if v == "age65plus":  # share x the group's (fixed) residents; CI from the share's CI
                            sh[f"est_persons|{g}"] = sh[g] * D[g][0].sum()
                est[tag] = sh
            for k, x in est[""].items():
                lo, hi = pct(x[1:])
                lo2, hi2 = pct(est["_acs"][k][1:])
                if "|" in k:
                    measure, group = k.split("|")
                    if measure == "expected":
                        measure, group = "share", f"expected for {GROUP_LABEL[group]} (region x rural mix)"
                    elif measure == "est_persons":
                        group = GROUP_LABEL[group]
                else:
                    measure, group = "share", GROUP_LABEL[k]
                base_g = k.split("|")[-1].split("/")[0].split("-")[0]
                rows.append(dict(scope=scope, variable=v, unit=unit, measure=measure, group=group, value=x[0],
                                 ci_lo=lo, ci_hi=hi, ci_lo_acs=lo2, ci_hi_acs=hi2,
                                 n_tracts=n_tr.get(base_g, np.nan), weight=D[base_g][0].sum()
                                 if base_g in D else np.nan, release=RELEASE, B=B, seed=SEED))
    # robustness: pool ACS counts instead of averaging rates (each tract's households / persons scaled by its
    # pocket share), so household rates are weighted by households rather than residents
    sp = np.clip(d.est_pop_pocket50.values / d.pop_total.values, 0, 1)
    for v, (num, den, unit) in VARS.items():
        ok = (d[v + "_den"] > 0).values
        nn = np.where(ok, d[v].values * d[v + "_den"].values, 0.0)
        dd = np.where(ok, d[v + "_den"].values, 0.0)
        cs = {k: W1 @ np.bincount(tc, x, n_c) for k, x in
              [("pn", sp * nn), ("pd", sp * dd), ("an", nn), ("ad", dd)]}
        sh_p, sh_a = cs["pn"] / cs["pd"], cs["an"] / cs["ad"]
        for group, x in [("pocket residents (pooled ACS counts)", sh_p), ("all residents (pooled ACS counts)", sh_a),
                         ("pocket/all (pooled ACS counts)", sh_p / sh_a)]:
            lo, hi = pct(x[1:])
            rows.append(dict(scope="all tracts", variable=v, unit=unit, measure="ratio" if "/" in group else "share",
                             group=group, value=x[0], ci_lo=lo, ci_hi=hi, n_tracts=int(ok.sum()), release=RELEASE,
                             B=B, seed=SEED))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "acs_evacuation.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    show = R.copy()
    pctm = show.measure == "share"
    for c in ["value", "ci_lo", "ci_hi", "ci_lo_acs", "ci_hi_acs"]:
        show.loc[pctm, c] = show.loc[pctm, c] * 100
    main_sc = show.scope.isin(["all tracts", "WHP Q4", "WHP Q1-Q3", "rural", "urban"])
    print(show.loc[main_sc & (show.measure != "est_persons"),
                   ["scope", "variable", "measure", "group", "value", "ci_lo", "ci_hi", "ci_lo_acs", "ci_hi_acs",
                    "n_tracts"]].round(3).to_string(index=False))
    print(show.loc[show.measure == "est_persons", ["scope", "group", "value", "ci_lo", "ci_hi"]].round(0)
          .to_string(index=False))
    reg = show[show.scope.str.startswith("region") & show.measure.isin(["share", "ratio"])
               & show.group.isin(["pocket residents", "all residents", "pocket/all", "pocket/expected"])]
    print(reg[["scope", "variable", "measure", "group", "value", "ci_lo", "ci_hi"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
