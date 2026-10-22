"""Full action list: EVERY highest-hazard (WHP Q4, pooled quartiles) single-egress pocket that stays single-egress
when service/track/living_street/unknown roads are added, with county, location, size and what sits inside.

12_action_list.py keeps only the 60 largest of these (top_pockets_wildfire.csv). This script exports all of them,
so a county planner can take the whole list.

Selection (same rule as 03_stratify.py's top_pockets_wildfire.csv, without the head(60)):
  pocket >= 50 buildings, exit not via the region outline (~exit_outside_region), tract in WHP Q4 (pooled
  quartiles of tract_metrics_all), and ~ext_second_exit (at most 50% of pocket nodes join the extended core).
Facilities: Overture POIs of out/pocket_facilities.csv (fire/EMS/police/hospital, schools/child care, care homes;
  confidence >= 0.5, de-duplicated) matched to the pocket on (region, single-link road name, pocket building count);
  where several pockets share that key, the one with the nearest centroid.
OSM class: from out/external_validation.csv where the pocket was checked (a/b/c, §4.10); blank otherwise.
County names: US Census county FIPS codes (state lists below, in FIPS order); no external data is read.

Inputs : discovery/out/<region>_pockets.parquet (02_egress.py), tract_metrics_all.parquet (03_stratify.py),
         pocket_facilities.csv (05_pocket_facilities.py), external_validation.csv (07_external_osm.py)
Outputs: discovery/out/action_list_full.csv       one row per pocket
         discovery/out/action_list_full.geojson   same, as points (pocket centroid, EPSG:4326)
         discovery/out/action_list_by_county.csv  pockets, buildings, residents, facilities per county
Usage  : python -I discovery/15_full_action_list.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "out"
REGIONS = ["eastern-ok", "maricopa-az", "northern-ca", "south-central-tx", "eastern-wa"]

# County names by FIPS. For CA, OK, TX and WA the county codes are 001, 003, ... in this order; AZ is explicit.
_SEQ = {
    "06": ("CA", "Alameda,Alpine,Amador,Butte,Calaveras,Colusa,Contra Costa,Del Norte,El Dorado,Fresno,Glenn,"
                 "Humboldt,Imperial,Inyo,Kern,Kings,Lake,Lassen,Los Angeles,Madera,Marin,Mariposa,Mendocino,Merced,"
                 "Modoc,Mono,Monterey,Napa,Nevada,Orange,Placer,Plumas,Riverside,Sacramento,San Benito,"
                 "San Bernardino,San Diego,San Francisco,San Joaquin,San Luis Obispo,San Mateo,Santa Barbara,"
                 "Santa Clara,Santa Cruz,Shasta,Sierra,Siskiyou,Solano,Sonoma,Stanislaus,Sutter,Tehama,Trinity,"
                 "Tulare,Tuolumne,Ventura,Yolo,Yuba", 58),
    "40": ("OK", "Adair,Alfalfa,Atoka,Beaver,Beckham,Blaine,Bryan,Caddo,Canadian,Carter,Cherokee,Choctaw,Cimarron,"
                 "Cleveland,Coal,Comanche,Cotton,Craig,Creek,Custer,Delaware,Dewey,Ellis,Garfield,Garvin,Grady,Grant,"
                 "Greer,Harmon,Harper,Haskell,Hughes,Jackson,Jefferson,Johnston,Kay,Kingfisher,Kiowa,Latimer,"
                 "Le Flore,Lincoln,Logan,Love,McClain,McCurtain,McIntosh,Major,Marshall,Mayes,Murray,Muskogee,Noble,"
                 "Nowata,Okfuskee,Oklahoma,Okmulgee,Osage,Ottawa,Pawnee,Payne,Pittsburg,Pontotoc,Pottawatomie,"
                 "Pushmataha,Roger Mills,Rogers,Seminole,Sequoyah,Stephens,Texas,Tillman,Tulsa,Wagoner,Washington,"
                 "Washita,Woods,Woodward", 77),
    "48": ("TX", "Anderson,Andrews,Angelina,Aransas,Archer,Armstrong,Atascosa,Austin,Bailey,Bandera,Bastrop,Baylor,"
                 "Bee,Bell,Bexar,Blanco,Borden,Bosque,Bowie,Brazoria,Brazos,Brewster,Briscoe,Brooks,Brown,Burleson,"
                 "Burnet,Caldwell,Calhoun,Callahan,Cameron,Camp,Carson,Cass,Castro,Chambers,Cherokee,Childress,Clay,"
                 "Cochran,Coke,Coleman,Collin,Collingsworth,Colorado,Comal,Comanche,Concho,Cooke,Coryell,Cottle,"
                 "Crane,Crockett,Crosby,Culberson,Dallam,Dallas,Dawson,Deaf Smith,Delta,Denton,DeWitt,Dickens,"
                 "Dimmit,Donley,Duval,Eastland,Ector,Edwards,Ellis,El Paso,Erath,Falls,Fannin,Fayette,Fisher,Floyd,"
                 "Foard,Fort Bend,Franklin,Freestone,Frio,Gaines,Galveston,Garza,Gillespie,Glasscock,Goliad,"
                 "Gonzales,Gray,Grayson,Gregg,Grimes,Guadalupe,Hale,Hall,Hamilton,Hansford,Hardeman,Hardin,Harris,"
                 "Harrison,Hartley,Haskell,Hays,Hemphill,Henderson,Hidalgo,Hill,Hockley,Hood,Hopkins,Houston,Howard,"
                 "Hudspeth,Hunt,Hutchinson,Irion,Jack,Jackson,Jasper,Jeff Davis,Jefferson,Jim Hogg,Jim Wells,"
                 "Johnson,Jones,Karnes,Kaufman,Kendall,Kenedy,Kent,Kerr,Kimble,King,Kinney,Kleberg,Knox,Lamar,Lamb,"
                 "Lampasas,La Salle,Lavaca,Lee,Leon,Liberty,Limestone,Lipscomb,Live Oak,Llano,Loving,Lubbock,Lynn,"
                 "McCulloch,McLennan,McMullen,Madison,Marion,Martin,Mason,Matagorda,Maverick,Medina,Menard,Midland,"
                 "Milam,Mills,Mitchell,Montague,Montgomery,Moore,Morris,Motley,Nacogdoches,Navarro,Newton,Nolan,"
                 "Nueces,Ochiltree,Oldham,Orange,Palo Pinto,Panola,Parker,Parmer,Pecos,Polk,Potter,Presidio,Rains,"
                 "Randall,Reagan,Real,Red River,Reeves,Refugio,Roberts,Robertson,Rockwall,Runnels,Rusk,Sabine,"
                 "San Augustine,San Jacinto,San Patricio,San Saba,Schleicher,Scurry,Shackelford,Shelby,Sherman,"
                 "Smith,Somervell,Starr,Stephens,Sterling,Stonewall,Sutton,Swisher,Tarrant,Taylor,Terrell,Terry,"
                 "Throckmorton,Titus,Tom Green,Travis,Trinity,Tyler,Upshur,Upton,Uvalde,Val Verde,Van Zandt,"
                 "Victoria,Walker,Waller,Ward,Washington,Webb,Wharton,Wheeler,Wichita,Wilbarger,Willacy,Williamson,"
                 "Wilson,Winkler,Wise,Wood,Yoakum,Young,Zapata,Zavala", 254),
    "53": ("WA", "Adams,Asotin,Benton,Chelan,Clallam,Clark,Columbia,Cowlitz,Douglas,Ferry,Franklin,Garfield,Grant,"
                 "Grays Harbor,Island,Jefferson,King,Kitsap,Kittitas,Klickitat,Lewis,Lincoln,Mason,Okanogan,Pacific,"
                 "Pend Oreille,Pierce,San Juan,Skagit,Skamania,Snohomish,Spokane,Stevens,Thurston,Wahkiakum,"
                 "Walla Walla,Whatcom,Whitman,Yakima", 39),
}
COUNTY = {"04001": "Apache AZ", "04003": "Cochise AZ", "04005": "Coconino AZ", "04007": "Gila AZ",
          "04009": "Graham AZ", "04011": "Greenlee AZ", "04012": "La Paz AZ", "04013": "Maricopa AZ",
          "04015": "Mohave AZ", "04017": "Navajo AZ", "04019": "Pima AZ", "04021": "Pinal AZ",
          "04023": "Santa Cruz AZ", "04025": "Yavapai AZ", "04027": "Yuma AZ"}
for _st, (_ab, _names, _n) in _SEQ.items():
    _l = _names.split(",")
    assert len(_l) == _n, (_ab, len(_l))
    COUNTY.update({f"{_st}{2 * i + 1:03d}": f"{nm} {_ab}" for i, nm in enumerate(_l)})
# spot checks against the hand-written table of 12_action_list.py
assert COUNTY["48471"] == "Walker TX" and COUNTY["40145"] == "Wagoner OK" and COUNTY["06103"] == "Tehama CA"
assert COUNTY["53057"] == "Skagit WA" and COUNTY["48363"] == "Palo Pinto TX" and COUNTY["40099"] == "Murray OK"

FAC_GROUPS = ["emergency", "schools_childcare", "care_homes"]


def main():
    t = pd.read_parquet(OUT / "tract_metrics_all.parquet", columns=["GEOID", "whp_q", "whp_q_region", "ur"])
    P = pd.concat([pd.read_parquet(OUT / f"{r}_pockets.parquet") for r in REGIONS], ignore_index=True)
    P = P[(P.n_bldg >= 50) & ~P.exit_outside_region].merge(t, on="GEOID", how="left")
    n_q4 = int((P.whp_q == "Q4").sum())
    A = P[(P.whp_q == "Q4") & ~P.ext_second_exit].copy()
    A["county_fips"] = A.GEOID.str[:5]
    A["county"] = A.county_fips.map(COUNTY)
    assert A.county.notna().all(), A.loc[A.county.isna(), "county_fips"].unique()

    # facilities inside: candidates share (region, single-link road name, pocket building count); when several
    # pockets share that key, the facility goes to the pocket with the nearest centroid
    f = pd.read_csv(OUT / "pocket_facilities.csv", dtype={"GEOID": str})
    f = f[f.in_pocket50 == True].copy()  # noqa: E712
    f["_road"] = f.egress_road.fillna("")
    f = f.rename(columns={"pocket_bldg": "n_bldg"})
    P["_road"] = P.egress_road.fillna("")
    key = ["region", "_road", "n_bldg"]
    m = f[["id", "name", "group", "x", "y"] + key].merge(P[key + ["pocket", "lat", "lon"]], on=key, how="inner")
    m["d"] = np.hypot((m.x - m.lon) * np.cos(np.radians(m.lat)), m.y - m.lat)
    n_cand = m.groupby("id").pocket.transform("size")
    n_amb = m.loc[n_cand > 1, "id"].nunique()
    m = m.sort_values(["id", "d"]).drop_duplicates("id")
    print(f"facilities in pockets >= 50 bldg: {len(f)}; matched to a pocket: {len(m)} "
          f"({n_amb} with several candidate pockets on the shared key, resolved by nearest centroid)")
    fac = m.groupby(["region", "pocket"]).agg(
        facilities_inside=("name", lambda s: "; ".join(sorted(set(map(str, s))))), n_facilities=("name", "size"))
    for g in FAC_GROUPS:
        fac[f"n_{g}"] = m[m.group == g].groupby(["region", "pocket"]).size()
    A = A.merge(fac, left_on=["region", "pocket"], right_index=True, how="left")
    for c in ["n_facilities"] + [f"n_{g}" for g in FAC_GROUPS]:
        A[c] = A[c].fillna(0).astype(int)
    print(f"facilities attributed to listed pockets: {int(A.n_facilities.sum())}")

    # OSM validation
    e = pd.read_csv(OUT / "external_validation.csv", dtype={"GEOID": str})
    e = e[["region", "pocket", "sample", "label", "osm_class", "osm_subclass"]].rename(
        columns={"sample": "osm_sample", "label": "osm_label"})
    A = A.merge(e, on=["region", "pocket"], how="left")

    A["partial_second_link_share"] = A.ext_core_share.round(3)
    A = A.sort_values(["county", "est_pop"], ascending=[True, False]).reset_index(drop=True)
    A["est_pop"] = A.est_pop.round(0).astype(int)
    A["road_km"] = A.road_km.round(2)
    cols = ["county", "county_fips", "region", "GEOID", "pocket", "locality", "lat", "lon", "n_bldg", "est_pop",
            "road_km", "egress_road", "egress_class", "partial_second_link_share", "near_miss", "usfs_WHP_mean",
            "whp_q_region", "ur", "svi_overall", "aiannh_name", "n_facilities", "n_emergency",
            "n_schools_childcare", "n_care_homes", "facilities_inside", "osm_class", "osm_subclass", "osm_sample",
            "osm_label"]
    A["usfs_WHP_mean"] = A.usfs_WHP_mean.round(1)
    A["lat"], A["lon"] = A.lat.round(5), A.lon.round(5)
    A[cols].to_csv(OUT / "action_list_full.csv", index=False)

    feats = []
    for r in A[cols].itertuples(index=False):
        props = {k: (None if (isinstance(v, float) and np.isnan(v)) else
                     (v.item() if hasattr(v, "item") else v)) for k, v in zip(cols, r)}
        feats.append({"type": "Feature", "geometry": {"type": "Point", "coordinates": [props["lon"], props["lat"]]},
                      "properties": props})
    gj = {"type": "FeatureCollection", "name": "action_list_full",
          "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:OGC:1.3:CRS84"}}, "features": feats}
    (OUT / "action_list_full.geojson").write_text(json.dumps(gj, ensure_ascii=False), encoding="utf-8")

    C = (A.assign(fac=A.n_facilities, emer=A.n_emergency, sch=A.n_schools_childcare, care=A.n_care_homes,
                  with_fac=A.n_facilities > 0, ge1000=A.est_pop >= 1000)
         .groupby(["county", "county_fips", "region"])
         .agg(pockets=("pocket", "size"), buildings=("n_bldg", "sum"), est_residents=("est_pop", "sum"),
              pockets_ge_1000_residents=("ge1000", "sum"), largest_pocket_residents=("est_pop", "max"),
              pockets_with_facilities=("with_fac", "sum"), facilities=("fac", "sum"), emergency=("emer", "sum"),
              schools_childcare=("sch", "sum"), care_homes=("care", "sum"))
         .reset_index().sort_values("est_residents", ascending=False))
    C.to_csv(OUT / "action_list_by_county.csv", index=False)

    print(f"WHP-Q4 pockets: {n_q4}; still single-egress with service/track: {len(A)} pockets, "
          f"{A.n_bldg.sum():,} buildings, ~{A.est_pop.sum():,} residents, {A.county.nunique()} counties "
          f"in {A.region.nunique()} regions")
    print("per region:\n", A.groupby("region").agg(pockets=("pocket", "size"), residents=("est_pop", "sum"),
                                                   counties=("county", "nunique")).to_string())
    print(f"pockets >= 1000 residents: {int((A.est_pop >= 1000).sum())}; with >= 1 facility: "
          f"{int((A.n_facilities > 0).sum())} (emergency {int((A.n_emergency > 0).sum())}, schools/child care "
          f"{int((A.n_schools_childcare > 0).sum())}, care homes {int((A.n_care_homes > 0).sum())})")
    print("OSM-checked in list:", A.osm_class.value_counts().to_dict())
    k = C.est_residents.cumsum() / C.est_residents.sum()
    print(f"counties holding 50% of residents: {int((k < 0.5).sum()) + 1}; 80%: {int((k < 0.8).sum()) + 1}")
    print(C.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
