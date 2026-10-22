"""Face-validity checks on named places (looked up by Overture POI locality strings only, no external data).

Prints, for a few well-known places, what the Overture graph says: exits of the largest settlement cluster,
single-egress pockets, and roadless islands. Output: discovery/out/validation_named_places.csv
"""
from __future__ import annotations

import pandas as pd

from common import OUT, REGIONS

CHECKS = [  # (region, locality string as it appears in Overture POI addresses, what is publicly known)
    ("northern-ca", "Paradise", "Camp Fire 2018 town; several evacuation routes (Skyway, Clark, Pentz, Neal)"),
    ("northern-ca", "Magalia", "ridge community north of Paradise; few routes out"),
    ("northern-ca", "Olympic Valley", "Palisades Tahoe valley; one access road (Olympic Valley Rd)"),
    ("northern-ca", "Corning", "Rancho Tehama Reserve; one access road"),
    ("northern-ca", "Cottonwood", "Lake California subdivision; one main access road"),
    ("eastern-wa", "Stehekin", "no road connection to the outside (boat / plane / trail only)"),
    ("eastern-wa", "Leavenworth", "Icicle Creek canyon; one road in"),
    ("eastern-wa", "Quincy", "Crescent Bar; one access road"),
]


def main():
    C = pd.concat([pd.read_parquet(OUT / f"{r}_clusters.parquet") for r in REGIONS], ignore_index=True)
    P = pd.concat([pd.read_parquet(OUT / f"{r}_pockets.parquet") for r in REGIONS], ignore_index=True)
    I = pd.concat([pd.read_parquet(OUT / f"{r}_islands.parquet") for r in REGIONS], ignore_index=True)
    rows = []
    for region, loc, known in CHECKS:
        c = C[(C.region == region) & (C.locality == loc)].sort_values("n_bldg", ascending=False).head(1)
        p = P[(P.region == region) & (P.locality == loc) & ~P.exit_outside_region].sort_values("n_bldg", ascending=False).head(2)
        i = I[(I.region == region) & (I.locality == loc)].sort_values("n_bldg", ascending=False).head(1)
        rows.append(dict(
            region=region, locality=loc, known=known,
            largest_cluster_bldg=int(c.n_bldg.iloc[0]) if len(c) else None,
            largest_cluster_exits=int(c.exits.iloc[0]) if len(c) else None,
            largest_cluster_exits_with_service_track=int(c.exits_ext.iloc[0]) if len(c) else None,
            top_pockets="; ".join(f"{int(r.n_bldg)} bldg via {r.egress_road} ({r.egress_class}), GEOID {r.GEOID}"
                                  for r in p.itertuples()),
            island="; ".join(f"{int(r.n_bldg)} bldg, {r.len_km:.1f} km road, category={r.category}, gap {r.gap_m / 1000:.1f} km"
                             for r in i.itertuples()),
        ))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "validation_named_places.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    print(out.to_string())


if __name__ == "__main__":
    main()
