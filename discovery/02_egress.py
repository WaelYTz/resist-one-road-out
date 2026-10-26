"""Single-egress communities in the Overture drivable road graph.

Graph: nodes = Overture connectors, edges = the pieces of each drivable segment between consecutive connectors
(ordered by `at`). A virtual OUT node is joined to every connector within EDGE_BUF_M of the region outline (or
outside it), so roads that leave the extract are treated as possible exits (conservative: fewer false
single-egress results near the clip edge).

Two complementary measures:
  A. Pockets (exact, graph-native). core = largest 2-edge-connected component (bridges removed; OUT included).
     A pocket is a connected set of non-core nodes; by construction it hangs off the core by exactly ONE edge
     (a bridge): one road segment whose closure cuts every building in the pocket off the meshed network.
  B. Settlement clusters. Buildings (bbox centres) are binned on a 200 m grid; cells with >= 2 buildings are
     8-connected into clusters; clusters with >= 50 buildings are kept. exits = max number of edge-disjoint
     paths (unit-capacity max-flow) from the graph edges the cluster's buildings attach to, to core nodes
     farther than R_M from every cluster cell (or to a road leaving a 10 km local window / the region).
Both are recomputed on the extended graph (drivable + service/track/unknown/living_street from
-roads-unfiltered) to see whether a second exit exists only as a low-class road.

Usage: python discovery/02_egress.py <region>   (needs discovery/cache/<region>_bldg.npz from 01_islands.py)
"""
from __future__ import annotations

import sys
import time

import numpy as np
import pandas as pd
import shapely
from scipy import ndimage
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components, maximum_flow
from scipy.spatial import cKDTree

from common import (CACHE, DRIVABLE, EDGE_BUF_M, EXTRA_VEHICLE, OUT, con, find_bridges, graph_codes, load_tracts,
                    log, nearest_locality, p, region_outline, to5070, wkb_to_geoms, xy_to4326)

R_M = 3000.0
WIN_M = 10000.0
CELL_M = 200.0
MIN_CELL_BLDG = 2
MIN_CLUSTER_BLDG = 50
BIG = 1_000_000
NEARMISS_M = 10.0
MIN_SHARE_ATTACHED = 0.5


def build_graph(c, region, layer, classes, tag, region_poly, pieces):
    si, ci, at, n = graph_codes(c, region, layer, classes, tag)
    w = c.execute(f"SELECT s.si, ST_AsWKB(r.geometry) w FROM {tag}_s s JOIN '{p(region, layer)}' r USING(id) "
                  "ORDER BY s.si").fetchdf()
    assert (w.si.values == np.arange(n)).all()
    geoms = to5070(wkb_to_geoms(w.w))
    del w
    length = shapely.length(geoms)
    nn = int(ci.max()) + 1
    first = np.unique(ci, return_index=True)[1]
    pts = shapely.line_interpolate_point(geoms[si[first]], at[first], normalized=True)
    nx_ = np.full(nn, np.nan)
    ny_ = np.full(nn, np.nan)
    nx_[ci[first]] = shapely.get_x(pts)
    ny_[ci[first]] = shapely.get_y(pts)
    same = si[1:] == si[:-1]
    eu = ci[:-1][same]
    ev = ci[1:][same]
    eseg = si[:-1][same]
    ea0 = at[:-1][same]
    elen = (at[1:] - at[:-1])[same] * length[eseg]
    # OUT node
    ptree = shapely.STRtree(pieces)
    npts = shapely.points(nx_, ny_)
    ii, _ = ptree.query(npts, predicate="dwithin", distance=EDGE_BUF_M)  # (input = node, tree = outline piece)
    nearb = np.zeros(nn, bool)
    nearb[ii] = True
    shapely.prepare(region_poly)
    outside = ~shapely.contains_xy(region_poly, nx_, ny_)
    bnodes = np.where(nearb | outside)[0]
    OUTN = nn
    u = np.r_[eu, bnodes]
    v = np.r_[ev, np.full(len(bnodes), OUTN)]
    is_virtual = np.r_[np.zeros(len(eu), bool), np.ones(len(bnodes), bool)]
    N = nn + 1
    br = find_bridges(N, u, v)
    keep = ~br & (u != v)
    A = coo_matrix((np.ones(keep.sum()), (u[keep], v[keep])), shape=(N, N))
    _, tecc = connected_components(A, directed=False)
    sizes = np.bincount(tecc)
    core_lab = {int(np.argmax(sizes)), int(tecc[OUTN])}
    core = np.isin(tecc, list(core_lab))
    # connected components of the full graph (to know what is attached to the core at all)
    A2 = coo_matrix((np.ones(len(u)), (u, v)), shape=(N, N))
    _, cc = connected_components(A2, directed=False)
    core_cc = cc[np.where(core)[0][0]]
    # pockets: components of the graph induced on non-core nodes that are connected to the core
    nc = ~core & (cc == core_cc)
    m_in = nc[u] & nc[v]
    A3 = coo_matrix((np.ones(m_in.sum()), (u[m_in], v[m_in])), shape=(N, N))
    _, pk = connected_components(A3, directed=False)
    pocket = np.where(nc, pk, -1)
    log(f"graph {tag}: nodes {nn} edges {len(eu)} boundary-nodes {len(bnodes)} bridges {int(br[~is_virtual].sum())} "
        f"core nodes {int(core.sum())} ({core.mean():.3f}) pockets {len(np.unique(pocket[pocket >= 0]))}")
    cid = c.execute(f"SELECT ci, cid FROM {tag}_c ORDER BY ci").fetchdf()
    return dict(n=n, geoms=geoms, length=length, nx=nx_, ny=ny_, eu=eu, ev=ev, eseg=eseg, ea0=ea0, elen=elen,
                u=u, v=v, virtual=is_virtual, bridge=br, core=core, pocket=pocket, N=N, OUTN=OUTN, cc=cc,
                core_cc=core_cc, cid=cid.cid.values, bnodes=bnodes)


def attach_buildings(G, bx, by, bseg):
    """For each building with a nearest drivable segment, the graph edge (index into eu/ev) it sits on."""
    ok = np.where(bseg >= 0)[0]
    f = shapely.line_locate_point(G["geoms"][bseg[ok]], shapely.points(bx[ok], by[ok]), normalized=True)
    key_e = G["eseg"] + 0.5 * G["ea0"]
    key_b = bseg[ok] + 0.5 * f
    idx = np.searchsorted(key_e, key_b, side="right") - 1
    idx = np.clip(idx, 0, len(key_e) - 1)
    bad = G["eseg"][idx] != bseg[ok]
    # position before first connector of the segment -> first edge of the segment
    nxt = np.clip(idx + 1, 0, len(key_e) - 1)
    fix = bad & (G["eseg"][nxt] == bseg[ok])
    idx[fix] = nxt[fix]
    bad = G["eseg"][idx] != bseg[ok]
    out = np.full(len(bx), -1, np.int64)
    out[ok[~bad]] = idx[~bad]
    return out


def cluster_buildings(bx, by):
    x0, y0 = bx.min(), by.min()
    ix = ((bx - x0) // CELL_M).astype(np.int64)
    iy = ((by - y0) // CELL_M).astype(np.int64)
    W, H = ix.max() + 1, iy.max() + 1
    cnt = np.zeros((H, W), np.int32)
    np.add.at(cnt, (iy, ix), 1)
    lab, nlab = ndimage.label(cnt >= MIN_CELL_BLDG, structure=np.ones((3, 3), bool))
    blab = lab[iy, ix] - 1  # -1 = not in a cluster
    sizes = np.bincount(blab[blab >= 0], minlength=nlab)
    keep = sizes >= MIN_CLUSTER_BLDG
    blab = np.where((blab >= 0) & keep[np.maximum(blab, 0)], blab, -1)
    # cell centres per cluster
    cy, cx = np.nonzero(lab)
    cl = lab[cy, cx] - 1
    m = keep[cl]
    cells = pd.DataFrame({"cl": cl[m], "x": x0 + (cx[m] + 0.5) * CELL_M, "y": y0 + (cy[m] + 0.5) * CELL_M})
    return blab, cells


def cluster_exits(G, att_nodes_by_cl, cells_by_cl):
    """Unit-capacity max-flow from each cluster's attach nodes to the far core within a local window.
    Local window = cluster extent + WIN_M; nodes are found through an x-sorted index and edges through a
    node->edge incidence CSR, so each cluster costs O(local graph) rather than O(region graph)."""
    N = G["N"]
    u, v = G["u"], G["v"]
    real = ~G["virtual"] & (u != v)
    ur, vr = u[real], v[real]
    m = len(ur)
    nxv, nyv = G["nx"], G["ny"]
    core = G["core"]
    on_boundary = np.zeros(N, bool)
    on_boundary[G["bnodes"]] = True
    xorder = np.argsort(nxv, kind="stable")
    xs = nxv[xorder]
    inc_n = np.r_[ur, vr]
    inc_e = np.r_[np.arange(m), np.arange(m)]
    o = np.argsort(inc_n, kind="stable")
    inc_e = inc_e[o]
    ptr = np.searchsorted(inc_n[o], np.arange(N + 1))
    inwin = np.zeros(N, bool)
    lid = np.full(N, -1, np.int64)
    res = {}
    for cl in sorted(att_nodes_by_cl):
        src_nodes = att_nodes_by_cl[cl]
        cx, cy = cells_by_cl[cl]
        x0, x1 = cx.min() - WIN_M, cx.max() + WIN_M
        y0, y1 = cy.min() - WIN_M, cy.max() + WIN_M
        lo, hi = np.searchsorted(xs, [x0, x1 + 1e-6])
        cand = xorder[lo:hi]
        loc = np.sort(cand[(nyv[cand] >= y0) & (nyv[cand] <= y1)])
        if len(loc) == 0:
            res[cl] = 0
            continue
        inwin[loc] = True
        lid[loc] = np.arange(len(loc))
        st_, en_ = ptr[loc], ptr[loc + 1]
        cnt = en_ - st_
        idx = np.repeat(st_ - np.r_[0, np.cumsum(cnt)[:-1]], cnt) + np.arange(cnt.sum())
        le = np.unique(inc_e[idx])
        eu_, ev_ = ur[le], vr[le]
        both = inwin[eu_] & inwin[ev_]
        cross_u = inwin[eu_] & ~inwin[ev_]
        cross_v = ~inwin[eu_] & inwin[ev_]
        S, T = len(loc), len(loc) + 1
        nb_ = int(both.sum())
        rows = [lid[eu_[both]], lid[ev_[both]], lid[eu_[cross_u]], lid[ev_[cross_v]]]
        cols = [lid[ev_[both]], lid[eu_[both]], np.full(cross_u.sum(), T), np.full(cross_v.sum(), T)]
        caps = [np.ones(2 * nb_ + int(cross_u.sum()) + int(cross_v.sum()), np.int64)]
        d, _ = cKDTree(np.column_stack([cx, cy])).query(np.column_stack([nxv[loc], nyv[loc]]), k=1,
                                                         distance_upper_bound=R_M)
        sl = lid[src_nodes]
        sl = np.unique(sl[sl >= 0])
        is_src = np.zeros(len(loc), bool)
        is_src[sl] = True
        far = np.isinf(d) & core[loc] & ~is_src
        rows.append(np.where(far)[0]); cols.append(np.full(far.sum(), T)); caps.append(np.full(far.sum(), BIG))
        ob = on_boundary[loc]
        rows.append(np.where(ob)[0]); cols.append(np.full(ob.sum(), T)); caps.append(np.ones(ob.sum(), np.int64))
        rows.append(np.full(len(sl), S)); cols.append(sl); caps.append(np.full(len(sl), BIG))
        M = csr_matrix((np.concatenate(caps).astype(np.int32), (np.concatenate(rows), np.concatenate(cols))),
                       shape=(len(loc) + 2, len(loc) + 2))
        M.sum_duplicates()
        M.setdiag(0)
        M.eliminate_zeros()
        res[cl] = int(min(maximum_flow(M, S, T).flow_value, BIG))
        inwin[loc] = False
        lid[loc] = -1
    return res


def main(region: str):
    T0 = time.time()
    c = con()
    tr, tg = load_tracts(region)
    region_poly, pieces = region_outline(tg)
    st = pd.read_parquet(p(region, "strata"), columns=["GEOID", "usfs_WHP_mean", "tribal_any", "aiannh_name",
                                                       "ur_class", "svi_overall"])
    st["GEOID"] = st.GEOID.astype(str)
    tr = tr.merge(st, on="GEOID", how="left")
    B = np.load(CACHE / f"{region}_bldg.npz")
    bx, by, btr, bseg = B["x"].astype(np.float64), B["y"].astype(np.float64), B["tract"], B["seg"]
    nb = len(bx)
    n_t = np.bincount(btr[btr >= 0], minlength=len(tr))
    per_b = np.where(n_t > 0, tr.pop_total.values / np.maximum(n_t, 1), 0.0)
    bpop = np.where(btr >= 0, per_b[np.maximum(btr, 0)], 0.0)

    G = build_graph(c, region, "roads", None, "d", region_poly, pieces)
    batt = attach_buildings(G, bx, by, bseg)
    log("buildings attached to an edge", int((batt >= 0).sum()), "of", nb)
    # pocket of each building: non-core endpoint of its edge
    pu = np.where(batt >= 0, G["pocket"][G["eu"][np.maximum(batt, 0)]], -1)
    pv = np.where(batt >= 0, G["pocket"][G["ev"][np.maximum(batt, 0)]], -1)
    bpocket = np.where(pu >= 0, pu, pv)

    # extended graph (adds service / track / unknown / living_street)
    GE = build_graph(c, region, "roads-unfiltered", DRIVABLE + EXTRA_VEHICLE, "e", region_poly, pieces)
    # map drivable connectors -> extended connectors (same connector ids)
    m = pd.Series(np.arange(len(GE["cid"])), index=GE["cid"])
    d2e = m.reindex(G["cid"]).fillna(-1).astype(np.int64).values
    e_core_of_dnode = np.where(d2e >= 0, GE["core"][np.maximum(d2e, 0)], False)

    # ---------------- A. pockets ----------------
    has = bpocket >= 0
    pk = pd.DataFrame({"pocket": bpocket[has], "pop": bpop[has], "t": btr[has]})
    P = pk.groupby("pocket").agg(n_bldg=("pop", "size"), est_pop=("pop", "sum"),
                                 tract_idx=("t", lambda s: np.bincount(s[s >= 0]).argmax() if (s >= 0).any() else -1))
    # pocket nodes, bridge edge, extended-graph status
    nodes = pd.DataFrame({"node": np.arange(G["N"] - 1), "pocket": G["pocket"][:-1]})
    nodes = nodes[nodes.pocket >= 0]
    nodes["x"] = G["nx"][nodes.node]
    nodes["y"] = G["ny"][nodes.node]
    nodes["ecore"] = e_core_of_dnode[nodes.node]
    pg = nodes.groupby("pocket").agg(n_nodes=("node", "size"), cx=("x", "mean"), cy=("y", "mean"),
                                     ext_core_share=("ecore", "mean"))
    P = P.join(pg, how="left")
    # bridge edge attaching each pocket (edge with one end core, other end in pocket)
    u, v = G["u"], G["v"]
    pc = np.where(G["pocket"][u] >= 0, G["pocket"][u], G["pocket"][v])
    att = (G["core"][u] ^ G["core"][v]) & (pc >= 0) & ~G["virtual"]
    att_df = pd.DataFrame({"pocket": pc[att], "edge": np.where(att)[0]})
    P["n_attach_edges"] = att_df.groupby("pocket").size().reindex(P.index).fillna(0)
    # pockets whose single link to the core is the virtual OUT edge: their only exit is a road leaving the extract
    # -> unverifiable, excluded from all counts (flag kept for transparency)
    P["exit_outside_region"] = P.n_attach_edges == 0
    first_edge = att_df.groupby("pocket").edge.first().reindex(P.index)
    names = c.execute(f"SELECT s.si, s.class, r.names.primary AS name FROM d_s s JOIN '{p(region, 'roads')}' r USING(id) "
                      "ORDER BY s.si").fetchdf()
    eseg_full = np.r_[G["eseg"], np.full(int(G["virtual"].sum()), -1)]
    bs = eseg_full[first_edge.fillna(0).astype(int).values]
    P["egress_road"] = np.where(first_edge.notna() & (bs >= 0), names.name.values[np.maximum(bs, 0)], None)
    P["egress_class"] = np.where(first_edge.notna() & (bs >= 0), names["class"].values[np.maximum(bs, 0)], None)
    # road km inside pocket
    pe = np.where(G["pocket"][G["eu"]] >= 0, G["pocket"][G["eu"]], G["pocket"][G["ev"]])
    P["road_km"] = pd.Series(G["elen"][pe >= 0] / 1000).groupby(pe[pe >= 0]).sum().reindex(P.index).fillna(0)
    P["ext_second_exit"] = P.ext_core_share > 0.5  # pocket becomes part of the meshed core once service/track added
    # topology-defect test: a dead-end node of the pocket lying within NEARMISS_M of a drivable segment that is not
    # part of the pocket (undershoot / missing connector => the pocket may be an artefact of the connector graph)
    real_u, real_v = G["eu"], G["ev"]
    deg = np.bincount(np.r_[real_u, real_v], minlength=G["N"])
    isb = np.zeros(G["N"], bool)
    isb[G["bnodes"]] = True
    dn = np.where((deg == 1) & (G["pocket"] >= 0) & ~isb)[0]
    seg_pocket = pd.DataFrame({"s": np.r_[G["eseg"], G["eseg"]], "pk": np.r_[G["pocket"][real_u], G["pocket"][real_v]]})
    seg_pocket = seg_pocket[seg_pocket.pk >= 0].drop_duplicates()
    stree = shapely.STRtree(G["geoms"])
    qi, qs = stree.query(shapely.points(G["nx"][dn], G["ny"][dn]), predicate="dwithin", distance=NEARMISS_M)
    pairs = pd.DataFrame({"node": dn[qi], "s": qs, "pk": G["pocket"][dn[qi]]})
    own = pairs.merge(seg_pocket, on=["s", "pk"], how="left", indicator=True)
    pairs = own[own["_merge"] == "left_only"].copy()
    pairs["d"] = shapely.distance(shapely.points(G["nx"][pairs.node], G["ny"][pairs.node]), G["geoms"][pairs.s])
    P["dangle_gap_m"] = pairs.groupby("pk").d.min().reindex(P.index)
    P["near_miss"] = P.dangle_gap_m.notna()
    # distance of the pocket's nodes to the region outline (clip-edge caveat)
    ptree2 = shapely.STRtree(pieces)
    (qn, _), dn_ = ptree2.query_nearest(shapely.points(nodes.x.values, nodes.y.values), return_distance=True,
                                        all_matches=False)
    bd = np.full(len(nodes), np.nan)
    bd[qn] = dn_
    nodes["bd"] = bd
    P["min_dist_outline_km"] = nodes.groupby("pocket").bd.min().reindex(P.index) / 1000
    P = P.reset_index()
    P["GEOID"] = [tr.GEOID.iloc[int(k)] if k >= 0 else None for k in P.tract_idx]
    P = P.merge(tr.drop(columns=["pop_total"]), on="GEOID", how="left")
    lon, lat = xy_to4326(P.cx.values, P.cy.values)
    P["lon"], P["lat"] = lon, lat
    big = P.n_bldg >= 10
    P = P[big].copy()
    P["locality"] = nearest_locality(c, region, P.cx.values, P.cy.values)
    P["region"] = region
    P.to_parquet(OUT / f"{region}_pockets.parquet", index=False)
    v50 = (P.n_bldg >= 50) & ~P.exit_outside_region
    log("pockets >=10 bldg", len(P), ">=50 (verifiable)", int(v50.sum()), "bldg in pockets>=50",
        int(P.n_bldg[v50].sum()), "est pop", round(P.est_pop[v50].sum()), "exit-outside excluded",
        int(((P.n_bldg >= 50) & P.exit_outside_region).sum()))

    # ---------------- B. settlement clusters ----------------
    blab, cells = cluster_buildings(bx, by)
    ok = (blab >= 0) & (batt >= 0)
    att_nodes = pd.DataFrame({"cl": np.r_[blab[ok], blab[ok]],
                              "node": np.r_[G["eu"][batt[ok]], G["ev"][batt[ok]]]}).drop_duplicates()
    att_by = {k: g.node.values for k, g in att_nodes.groupby("cl")}
    cells_by = {k: (g.x.values, g.y.values) for k, g in cells.groupby("cl")}
    log("clusters", len(cells_by), "with attachments", len(att_by))
    ex = cluster_exits(G, att_by, cells_by)
    log("max-flow done (drivable)")
    # extended graph: buildings attach to the same drivable edge; map its end connectors into GE
    att_by_e = {}
    for k, nds in att_by.items():
        e_n = d2e[nds]
        att_by_e[k] = e_n[e_n >= 0]
    ex_e = cluster_exits(GE, att_by_e, cells_by)
    log("max-flow done (extended)")
    inc = blab >= 0
    C = pd.DataFrame({"cl": blab[inc], "pop": bpop[inc], "t": btr[inc], "x": bx[inc], "y": by[inc],
                      "pocket": bpocket[inc] >= 0, "att": batt[inc] >= 0})
    C = C.groupby("cl").agg(n_bldg=("pop", "size"), est_pop=("pop", "sum"), cx=("x", "mean"), cy=("y", "mean"),
                            share_in_pocket=("pocket", "mean"), share_attached=("att", "mean"),
                            tract_idx=("t", lambda s: np.bincount(s[s >= 0]).argmax() if (s >= 0).any() else -1))
    C["exits"] = pd.Series(ex).reindex(C.index).fillna(0).astype(int)
    # exits are only meaningful when most of the cluster's buildings sit within 150 m of a drivable segment;
    # otherwise the cluster's streets are missing from the drivable layer (reported separately as "unmapped")
    C["exits_valid"] = C.share_attached >= MIN_SHARE_ATTACHED
    C["exits_ext"] = pd.Series(ex_e).reindex(C.index).fillna(0).astype(int)
    C["n_cells"] = cells.groupby("cl").size().reindex(C.index)
    C = C.reset_index()
    C["GEOID"] = [tr.GEOID.iloc[int(k)] if k >= 0 else None for k in C.tract_idx]
    C = C.merge(tr.drop(columns=["pop_total"]), on="GEOID", how="left")
    lon, lat = xy_to4326(C.cx.values, C.cy.values)
    C["lon"], C["lat"] = lon, lat
    C["locality"] = nearest_locality(c, region, C.cx.values, C.cy.values)
    C["region"] = region
    C.to_parquet(OUT / f"{region}_clusters.parquet", index=False)
    V = C[C.exits_valid]
    log("clusters exits distribution (valid)", V.exits.clip(upper=4).value_counts().sort_index().to_dict(),
        "ext", V.exits_ext.clip(upper=4).value_counts().sort_index().to_dict(), "unmapped clusters",
        int((~C.exits_valid).sum()), "bldg", int(C.n_bldg[~C.exits_valid].sum()))

    # ---------------- per-tract ----------------
    T = tr[["GEOID", "pop_total"]].copy()
    T["region"] = region
    T["n_bldg"] = n_t
    psize = pd.Series(np.where(P.exit_outside_region, 0, P.n_bldg.values), index=P.pocket.values)
    bps = np.where(bpocket >= 0, psize.reindex(np.maximum(bpocket, 0)).fillna(0).values, 0)
    pext = pd.Series(P.ext_second_exit.values, index=P.pocket.values)
    bpext = np.where(bpocket >= 0, pext.reindex(np.maximum(bpocket, 0)).fillna(False).values.astype(bool), False)
    okb = btr >= 0
    def cnt(mask):
        return np.bincount(btr[okb & mask], minlength=len(tr))
    for th in (10, 50, 100):
        T[f"n_bldg_pocket{th}"] = cnt(bps >= th)
    T["n_bldg_pocket50_noext"] = cnt((bps >= 50) & ~bpext)
    cex = pd.Series(np.where(C.exits_valid, C.exits.values, -1), index=C.cl.values)
    cexe = pd.Series(np.where(C.exits_valid, C.exits_ext.values, -1), index=C.cl.values)
    cunm = pd.Series(~C.exits_valid.values, index=C.cl.values)
    bex = np.where(blab >= 0, cex.reindex(np.maximum(blab, 0)).fillna(-1).values, -1)
    bexe = np.where(blab >= 0, cexe.reindex(np.maximum(blab, 0)).fillna(-1).values, -1)
    T["n_bldg_cluster"] = cnt(blab >= 0)
    T["n_bldg_cluster_1exit"] = cnt(bex == 1)
    T["n_bldg_cluster_le2exit"] = cnt((bex >= 1) & (bex <= 2))
    T["n_bldg_cluster_1exit_ext"] = cnt(bexe == 1)
    bunm = np.where(blab >= 0, cunm.reindex(np.maximum(blab, 0)).fillna(False).values.astype(bool), False)
    T["n_bldg_cluster_unmapped"] = cnt(bunm)
    for col in ["n_bldg_pocket10", "n_bldg_pocket50", "n_bldg_pocket100", "n_bldg_pocket50_noext",
                "n_bldg_cluster_1exit", "n_bldg_cluster_le2exit", "n_bldg_cluster_1exit_ext", "n_bldg_cluster_unmapped"]:
        T[col.replace("n_bldg", "est_pop")] = T[col] * per_b
    T.to_parquet(OUT / f"{region}_tract_egress.parquet", index=False)
    # per-building pocket size (0 = not in a pocket), in the order of cache/<region>_bldg.npz: lets a point dataset
    # (e.g. CAL FIRE damage inspections, 19_dins.py) be matched to the nearest building and its pocket
    np.savez_compressed(CACHE / f"{region}_bldg_pocket.npz", bps=bps.astype(np.int32), bpocket=bpocket.astype(np.int32))
    log(region, "done", f"{time.time() - T0:.0f}s",
        T[["n_bldg", "n_bldg_pocket50", "est_pop_pocket50", "n_bldg_cluster_1exit", "est_pop_cluster_1exit"]].sum().round().to_dict())


if __name__ == "__main__":
    main(sys.argv[1])
