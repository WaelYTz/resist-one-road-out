"""Action list for emergency planners: high-hazard single-egress pockets per county, with what sits inside.

Input: out/top_pockets_wildfire.csv (WHP-Q4 pockets with no second link in any Overture vehicle class, from
02_egress / 03_stratify) and out/pocket_facilities.csv (Overture POIs inside pockets, from 05_pocket_facilities).
Output: out/action_list.csv and out/action_list.md (the table pasted into FINDINGS.md §5b).
Usage: python -I discovery/12_action_list.py
"""
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
COUNTY = {  # county FIPS -> name (US Census), for the counties that appear in the list
    "04013": "Maricopa AZ", "04019": "Pima AZ", "04021": "Pinal AZ", "04025": "Yavapai AZ",
    "06061": "Placer CA", "06089": "Shasta CA", "06103": "Tehama CA",
    "40037": "Creek OK", "40041": "Delaware OK", "40099": "Murray OK", "40145": "Wagoner OK",
    "48027": "Bell TX", "48049": "Brown TX", "48053": "Burnet TX", "48061": "Cameron TX", "48121": "Denton TX",
    "48321": "Matagorda TX", "48337": "Montague TX", "48339": "Montgomery TX", "48363": "Palo Pinto TX",
    "48367": "Parker TX", "48373": "Polk TX", "48423": "Smith TX", "48427": "Starr TX", "48439": "Tarrant TX",
    "48453": "Travis TX", "48471": "Walker TX",
    "53007": "Chelan WA", "53025": "Grant WA", "53041": "Lewis WA", "53053": "Pierce WA", "53057": "Skagit WA",
}
KEEP = {"fire_department", "ambulance_and_ems_services", "elementary_school", "middle_school", "high_school",
        "school", "private_school", "public_school", "child_care_and_day_care", "day_care_preschool", "preschool",
        "retirement_home", "assisted_living_facility", "nursing_home", "senior_citizen_services", "hospital"}


def main():
    t = pd.read_csv(OUT / "top_pockets_wildfire.csv", dtype={"GEOID": str})
    t["county"] = t["GEOID"].str[:5].map(COUNTY).fillna(t["GEOID"].str[:5])
    f = pd.read_csv(OUT / "pocket_facilities.csv", dtype={"GEOID": str})
    f = f[(f["in_pocket50"] == True) & f["cat"].isin(KEEP)]  # noqa: E712
    fac = (f.groupby(["GEOID", "egress_road"])["name"]
             .apply(lambda s: "; ".join(sorted(set(map(str, s)))[:4])).rename("facilities_inside"))
    t = t.merge(fac, left_on=["GEOID", "egress_road"], right_index=True, how="left")
    t = t.sort_values(["county", "est_pop"], ascending=[True, False])
    cols = ["county", "locality", "GEOID", "n_bldg", "est_pop", "egress_road", "egress_class", "lat", "lon",
            "usfs_WHP_mean", "aiannh_name", "facilities_inside"]
    t[cols].to_csv(OUT / "action_list.csv", index=False)
    top = t.groupby("county", sort=False).head(2).sort_values("est_pop", ascending=False).head(25)
    lines = ["| County | Community (locality) | Tract | Buildings | Est. residents | Single link | Lat, lon | Inside |",
             "|---|---|---|---|---|---|---|---|"]
    for r in top.itertuples():
        inside = r.facilities_inside if isinstance(r.facilities_inside, str) else "–"
        lines.append(f"| {r.county} | {r.locality} | {r.GEOID} | {r.n_bldg:,} | ~{round(r.est_pop, -1):,.0f} | "
                     f"{r.egress_road} ({r.egress_class}) | {r.lat:.4f}, {r.lon:.4f} | {inside} |")
    (OUT / "action_list.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(t)} pockets, {t['county'].nunique()} counties -> action_list.csv; top {len(top)} -> action_list.md")


if __name__ == "__main__":
    main()
