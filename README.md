# One road out: wildfire-country communities behind a single road segment

Best Bias Discovery entry for the Zindi **Bias Bounty Mapping Equity Challenge** (team **resist**). This repository
holds the full code, cached external checks and outputs for the bias discovery. The write-up itself was submitted
through the challenge's special-prize form; `discovery/FINDINGS.md` is the complete technical report.

## The finding

Overture's own road network leaves **2.88 million residents (7.1%) in 13,775 pockets** whose only connection to the
rest of the network is one road segment. The challenge's coverage scorecard measures road length and feature
counts, so it cannot see this.

- **It is concentrated where the wildfire hazard is.** 10.2% of residents in the highest wildfire-hazard quartile
  live in such pockets, against 3.3% in the lowest: **3.10x**. Adjusted for region and rural class it is 3.31x
  (county-block CI 2.29–4.87). Comparing only tracts of the same county it is 5.25x (3.82–7.27), so it is not a
  between-county effect. Within the same population-density tercile it is still 1.71x (1.23–2.52).
- **Checked against two outside sources.** In current OpenStreetMap, 29 of 33 random high-hazard pockets (88%, CI
  72–95%) have no public second road either. FHWA's public-road inventory (HPMS 2024) agrees in 31 of 33.
- **Fires that happened, not only a hazard map.** Within the same hazard quartile, pocket residents live 1.15x
  (1.09–1.22) more often in tracts that actually burned since 2015 (NIFC / MTBS perimeters).
- **Who is exposed.** Pocket residents are older: 65+ share 1.12x, about 486,000 people.
- **What to do.** `discovery/out/action_list_full.csv` lists 5,095 high-hazard single-exit pockets county by county
  (about 825,000 residents in 195 counties). For California, this is the list SB 99 asks counties to produce.

## Check every headline number in seconds

```bash
pip install pandas pyarrow
python discovery/verify_findings.py      # recomputes the 26 headline numbers from discovery/out: PASS 26/26
```

## Reproduce

Python 3.11 and the pinned packages in `requirements.txt`. The challenge data is public; place it under `data/` at
the repository root:

```bash
pip install -r requirements.txt awscli
aws s3 sync s3://us-west-2.opendata.source.coop/humane-intelligence/bias-bounty-mapping-equity-challenge/reference/ data/reference/ --no-sign-request
aws s3 sync s3://us-west-2.opendata.source.coop/humane-intelligence/bias-bounty-mapping-equity-challenge/strata/ data/strata/ --no-sign-request
```

Then run the scripts in the order given in `discovery/FINDINGS.md` §6, from the repository root. The OpenStreetMap,
ACS and HPMS checks rerun offline from the cached responses in `discovery/external/` (`--offline`). Runtimes,
seeds and parameters are in the same section, and §7 describes every output file.

## Data and licences

- **Challenge data** (Overture Maps extracts and strata tables): not included; downloaded as above. Outputs derived
  from Overture follow Overture's licences (ODbL for the transportation theme, CDLA-Permissive-2.0 for the others).
- **OpenStreetMap** (validation only): © OpenStreetMap contributors, available under the Open Database Licence
  (https://www.openstreetmap.org/copyright). Cached Overpass responses are in `discovery/external/osm_raw/`.
- **U.S. Census Bureau ACS 2020–2024** via Census Reporter, and **FHWA HPMS 2024**: U.S. Government works, public
  domain. Cached in `discovery/external/acs/` and `discovery/external/hpms_raw/`.
- No removed reference layer (TIGER roads, Microsoft buildings, HIFLD, CBP) is used anywhere.

The team's leaderboard model (Part I of the write-up) is not part of this repository.
