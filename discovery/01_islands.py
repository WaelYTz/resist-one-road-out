"""Routing islands in the Overture drivable road graph.

For one region:
  1. Build the drivable graph from <region>-overture-roads.parquet: segments sharing a connector_id are joined.
  2. Connected components; giant = component with most segments. "Edge" components have a segment within
     EDGE_BUF_M of the region outline (or outside it) and may be joined outside the extract -> never called islands.
     Island = non-giant, non-edge component.
  3. Re-test every island in the extended graph (drivable + service/track/unknown/living_street from
     -roads-unfiltered) and in the all-classes graph (adds footway/path/cycleway/...).
  4. Diagnose the cut: dangling endpoints of the island near a giant segment (undershoot / missing connector),
     geometric crossings with the giant without a shared connector, and the gap distance to the giant.
  5. Buildings (bbox centre of every Overture footprint) -> nearest drivable segment within 150 m (EPSG:5070)
     -> component -> per-tract counts and estimated residents (pop_total x share of tract buildings).

Usage: python discovery/01_islands.py <region>
Outputs: discovery/out/<region>_tract_islands.parquet, discovery/out/<region>_islands.parquet,
         discovery/cache/<region>_bldg.npz, discovery/cache/<region>_seg.npz
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import shapely

from common import (CACHE, DRIVABLE, EXTRA_VEHICLE, OUT, assign_points_to_tracts, components_from_connectors, con,
                    edge_flags, graph_codes, load_tracts, log, p, region_outline, to5070, wkb_to_geoms, xy_to4326,
                    xy_to5070)

NEAR_MAX_M = 150.0
NEARMISS_M = 10.0


def mode_str(s: pd.Series):
    s = s.dropna()
    s = s[s.astype(str).str.len() > 0]
    if s.empty:
        return None
    vc = s.value_counts()
    return vc.index[0]


def main(region: str):
    T0 = time.time()
    c = con()
    tr, tg = load_tracts(region)
    region_poly, pieces = region_outline(tg)
    log(region, "tracts", len(tr))

    # ---------- drivable graph ----------
    si, ci, at, n = graph_codes(c, region, "roads", None, "d")
    lab = components_from_connectors(si, ci, n)
    seg = c.execute("SELECT s.si, s.class, r.names.primary AS name, ST_AsWKB(r.geometry) w FROM d_s s "
                    f"JOIN '{p(region, 'roads')}' r USING(id) ORDER BY s.si").fetchdf()
    assert (seg.si.values == np.arange(n)).all()
    geoms = to5070(wkb_to_geoms(seg.w))
    seg = seg.drop(columns="w")
    length = shapely.length(geoms)
    edge, tree = edge_flags(geoms, region_poly, pieces)
    log("drivable segments", n, "components", lab.max() + 1, "edge segs", int(edge.sum()),
        "mid-outside", int((~shapely.contains_xy(region_poly, *shapely.get_coordinates(
            shapely.line_interpolate_point(geoms, 0.5, normalized=True)).T)).sum()))

    comp = pd.DataFrame({"comp": lab, "len": length, "edge": edge}).groupby("comp").agg(
        n_seg=("len", "size"), len_m=("len", "sum"), edge=("edge", "any"))
    giant = int(comp.n_seg.idxmax())
    comp["giant"] = comp.index == giant
    comp["island"] = (~comp.giant) & (~comp.edge)
    log("giant n_seg", comp.loc[giant, "n_seg"], f"({comp.loc[giant, 'n_seg'] / n:.4f})",
        "island comps", int(comp.island.sum()), "island segs", int(comp.n_seg[comp.island].sum()),
        "edge non-giant comps", int((comp.edge & ~comp.giant).sum()))

    # ---------- extended graphs (unfiltered layer) ----------
    res = {}
    for tag, classes in (("e", DRIVABLE + EXTRA_VEHICLE), ("a", None)):
        esi, eci, _, en = graph_codes(c, region, "roads-unfiltered", classes, tag)
        elab = components_from_connectors(esi, eci, en)
        m = c.execute(f"SELECT d.si dsi, e.si esi FROM d_s d JOIN {tag}_s e USING(id) ORDER BY d.si").fetchnumpy()
        mapping = np.full(n, -1, np.int64)
        mapping[m["dsi"]] = m["esi"]
        missing = int((mapping < 0).sum())
        dl = np.where(mapping >= 0, elab[np.maximum(mapping, 0)], -1)
        vals, cnt = np.unique(dl[lab == giant], return_counts=True)
        egiant = vals[np.argmax(cnt)]
        res[tag] = dict(dlab=dl, egiant=egiant, elab=elab, en=en, mapping=mapping)
        log(f"graph {tag}: segs {en} comps {elab.max() + 1} missing drivable ids {missing}")
        del esi, eci

    # edge test for extended components that contain islands: need geometry of their non-drivable segments
    island_comps = comp.index[comp.island].values
    isl_seg_mask = comp.island.values[lab]
    e = res["e"]
    ext_lab_of_island = pd.Series(e["dlab"][isl_seg_mask]).groupby(lab[isl_seg_mask]).first()
    cand_ext = np.setdiff1d(ext_lab_of_island.values, [e["egiant"]])
    ext_edge = set()
    # drivable edge segments inside those extended components
    for L in np.unique(e["dlab"][edge]):
        ext_edge.add(int(L))
    extra_idx = np.where(np.isin(e["elab"], cand_ext))[0]
    c.register("extra_idx", pd.DataFrame({"si": extra_idx.astype(np.int32)}))
    ex = c.execute(f"SELECT s.si, ST_AsWKB(r.geometry) w FROM e_s s JOIN extra_idx USING(si) "
                   f"JOIN '{p(region, 'roads-unfiltered')}' r USING(id) ORDER BY s.si").fetchdf()
    if len(ex):
        eg = to5070(wkb_to_geoms(ex.w))
        eedge, _ = edge_flags(eg, region_poly, pieces)
        for L in np.unique(e["elab"][ex.si.values[eedge]]):
            ext_edge.add(int(L))
    log("extended comps holding islands", len(cand_ext), "of which edge", len(set(cand_ext) & ext_edge))

    # ---------- diagnostics of each island ----------
    gmask = lab == giant
    gidx = np.where(gmask)[0]
    tree_g = shapely.STRtree(geoms[gidx])
    iidx = np.where(isl_seg_mask)[0]
    (q, _), dist = tree_g.query_nearest(geoms[iidx], return_distance=True, all_matches=False)
    gap = pd.Series(dist).groupby(lab[iidx[q]]).min()
    ii, _ = tree_g.query(geoms[iidx], predicate="intersects")
    cross = pd.Series(np.ones(len(ii))).groupby(lab[iidx[ii]]).sum()
    # dangling endpoints (connector used by a single drivable segment, at segment end)
    deg = np.bincount(ci)
    endrow = ((at <= 1e-9) | (at >= 1 - 1e-9)) & (deg[ci] == 1) & isl_seg_mask[si]
    if endrow.any():
        pts = shapely.line_interpolate_point(geoms[si[endrow]], at[endrow], normalized=True)
        (q2, _), d2 = tree_g.query_nearest(pts, return_distance=True, all_matches=False)
        dgap = pd.Series(d2).groupby(lab[si[endrow]][q2]).min()
    else:
        dgap = pd.Series(dtype=float)
    comp["gap_m"] = gap.reindex(comp.index)
    comp["cross_giant"] = cross.reindex(comp.index).fillna(0)
    comp["dangle_gap_m"] = dgap.reindex(comp.index)
    comp["ext_lab"] = pd.Series(e["dlab"]).groupby(lab).first().reindex(comp.index)
    comp["all_lab"] = pd.Series(res["a"]["dlab"]).groupby(lab).first().reindex(comp.index)
    comp["conn_ext"] = comp.ext_lab == e["egiant"]
    comp["ext_edge"] = comp.ext_lab.isin(list(ext_edge)) & ~comp.conn_ext
    comp["conn_all"] = comp.all_lab == res["a"]["egiant"]
    near_miss = (comp.dangle_gap_m <= NEARMISS_M) | (comp.cross_giant > 0)
    cat = np.select(
        [comp.conn_ext, comp.ext_edge, comp.conn_all, near_miss],
        ["class_filter", "ext_edge", "footpath_only", "near_miss"], "disconnected")
    comp["category"] = np.where(comp.island, cat, np.where(comp.giant, "giant", "edge"))
    log("island categories:", comp[comp.island].category.value_counts().to_dict())

    # ---------- buildings ----------
    b = c.execute(f"SELECT (bbox.xmin+bbox.xmax)/2 x, (bbox.ymin+bbox.ymax)/2 y FROM '{p(region, 'buildings')}'").fetchnumpy()
    bx, by = xy_to5070(b["x"], b["y"])
    del b
    nb = len(bx)
    btr = assign_points_to_tracts(bx, by, tg)
    log("buildings", nb, "unassigned to tract", int((btr < 0).sum()))
    bseg = np.full(nb, -1, np.int32)
    bd = np.full(nb, np.nan, np.float32)
    CH = 1_000_000
    for s in range(0, nb, CH):
        pts = shapely.points(bx[s:s + CH], by[s:s + CH])
        (qq, tt), dd = tree.query_nearest(pts, max_distance=NEAR_MAX_M, return_distance=True, all_matches=False)
        bseg[s + qq] = tt
        bd[s + qq] = dd
        log(f"  nearest {s + len(pts)}/{nb}")
    np.savez_compressed(CACHE / f"{region}_bldg.npz", x=bx.astype(np.float32), y=by.astype(np.float32),
                        tract=btr, seg=bseg, dist=bd)
    np.savez_compressed(CACHE / f"{region}_seg.npz", comp=lab, edge=edge, giant=giant, length=length)

    bcomp = np.where(bseg >= 0, lab[np.maximum(bseg, 0)], -1)
    bcat = np.where(bseg >= 0, comp.category.reindex(np.maximum(bcomp, 0)).values, "far")
    ntr = len(tr)
    ok = btr >= 0
    def cnt(mask):
        return np.bincount(btr[ok & mask], minlength=ntr)
    T = tr[["GEOID", "pop_total"]].copy()
    T["region"] = region
    T["n_bldg"] = cnt(np.ones(nb, bool))
    T["n_bldg_far150"] = cnt(bcat == "far")
    T["n_bldg_edgecomp"] = cnt(bcat == "edge")
    for k in ["class_filter", "ext_edge", "footpath_only", "near_miss", "disconnected"]:
        T[f"n_bldg_isl_{k}"] = cnt(bcat == k)
    T["n_bldg_island"] = T[[f"n_bldg_isl_{k}" for k in ["class_filter", "ext_edge", "footpath_only", "near_miss", "disconnected"]]].sum(axis=1)
    T["n_bldg_island_strict"] = T[[f"n_bldg_isl_{k}" for k in ["footpath_only", "near_miss", "disconnected"]]].sum(axis=1)
    per_b = np.where(T.n_bldg > 0, T.pop_total / T.n_bldg.clip(lower=1), 0.0)
    T["est_pop_island"] = T.n_bldg_island * per_b
    T["est_pop_island_strict"] = T.n_bldg_island_strict * per_b
    T["est_pop_far150"] = T.n_bldg_far150 * per_b
    T["share_bldg_island"] = T.n_bldg_island / T.n_bldg.clip(lower=1)
    T["share_bldg_island_strict"] = T.n_bldg_island_strict / T.n_bldg.clip(lower=1)
    T.to_parquet(OUT / f"{region}_tract_islands.parquet", index=False)

    # ---------- island table ----------
    isl = comp[comp.island].copy()
    bi = (bcomp >= 0) & comp.island.reindex(np.maximum(bcomp, 0)).values & ok
    pb = per_b[btr[bi]]
    isl["n_bldg"] = pd.Series(np.ones(bi.sum())).groupby(bcomp[bi]).sum().reindex(isl.index).fillna(0)
    isl["est_pop"] = pd.Series(pb).groupby(bcomp[bi]).sum().reindex(isl.index).fillna(0)
    maj = pd.Series(btr[bi]).groupby(bcomp[bi]).agg(lambda v: np.bincount(v).argmax())
    isl["tract_idx"] = maj.reindex(isl.index)
    mids = shapely.line_interpolate_point(geoms[iidx], 0.5, normalized=True)
    mx, my = shapely.get_x(mids), shapely.get_y(mids)
    cx = pd.Series(mx).groupby(lab[iidx]).mean().reindex(isl.index)
    cy = pd.Series(my).groupby(lab[iidx]).mean().reindex(isl.index)
    # tract of island centroid when no buildings
    noT = isl.tract_idx.isna()
    if noT.any():
        isl.loc[noT, "tract_idx"] = assign_points_to_tracts(cx[noT].values, cy[noT].values, tg)
    isl["GEOID"] = [tr.GEOID.iloc[int(k)] if k >= 0 else None for k in isl.tract_idx.fillna(-1)]
    lon, lat = xy_to4326(cx.values, cy.values)
    isl["lon"], isl["lat"] = lon, lat
    sdf = pd.DataFrame({"comp": lab[iidx], "name": seg.name.values[iidx], "class": seg["class"].values[iidx]})
    isl["road_names"] = sdf.groupby("comp").name.agg(lambda s: "; ".join(s.dropna().value_counts().index[:3])).reindex(isl.index)
    isl["classes"] = sdf.groupby("comp")["class"].agg(lambda s: ",".join(sorted(set(s)))).reindex(isl.index)
    # nearest POI locality
    pois = c.execute(f"SELECT ST_X(geometry) x, ST_Y(geometry) y, addresses[1].locality AS locality FROM '{p(region, 'pois')}' "
                     "WHERE addresses[1].locality IS NOT NULL ORDER BY id").fetchdf()
    px, py = xy_to5070(pois.x.values, pois.y.values)
    ptree = shapely.STRtree(shapely.points(px, py))
    (q3, t3), d3 = ptree.query_nearest(shapely.points(cx.values, cy.values), max_distance=8000, return_distance=True)
    loc = pd.Series(pois["locality"].values[t3]).groupby(q3).first()
    isl["locality"] = loc.reindex(np.arange(len(isl))).values
    isl["locality_dist_m"] = pd.Series(d3).groupby(q3).first().reindex(np.arange(len(isl))).values
    isl["region"] = region
    isl = isl.reset_index().rename(columns={"comp": "comp_id"})
    isl["len_km"] = isl.len_m / 1000
    isl.drop(columns=["len_m", "edge", "giant", "island"]).to_parquet(OUT / f"{region}_islands.parquet", index=False)
    # tract locality (mode of POI locality within tract) for naming
    ptr = assign_points_to_tracts(px, py, tg)
    tl = pd.DataFrame({"t": ptr, "loc": pois["locality"].values})
    tl = tl[tl.t >= 0].groupby("t")["loc"].agg(mode_str)
    T["tract_locality"] = tl.reindex(np.arange(ntr)).values
    T.to_parquet(OUT / f"{region}_tract_islands.parquet", index=False)

    s = T[["n_bldg", "n_bldg_island", "n_bldg_island_strict", "n_bldg_far150", "n_bldg_edgecomp", "est_pop_island",
           "est_pop_island_strict", "pop_total"]].sum()
    log(region, "summary:", s.round(0).to_dict(), f"elapsed {time.time() - T0:.0f}s")


if __name__ == "__main__":
    main(sys.argv[1])
