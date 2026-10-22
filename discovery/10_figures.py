"""Figures for FINDINGS.md (matplotlib, 150 dpi, light surface, the dataviz skill's reference palette).

  figs/pocket_rancho_tehama.png, figs/pocket_lake_california.png, figs/pocket_whitehorn_cove.png
      the pocket's roads and buildings (blue), the single link (orange), the meshed network (grey),
      the OSM-only second link if the external check found one (aqua, gates marked), tract outline, scale bar
  figs/whp_quartile_by_region.png
      share of residents in single-egress pockets (>= 50 buildings) by wildfire-hazard quartile, small multiples
      (pooled quartiles for all regions, within-region quartiles per region), 95% tract-bootstrap CIs.
      Table view: out/fig_whp_by_region.csv
Inputs: graph rebuilt as in 02_egress.py; out/<region>_pockets.parquet; out/external_validation.csv and the cached
Overpass responses in external/osm_raw/ (from 07_external_osm.py); out/tract_metrics_all.parquet.
Usage: python -I discovery/10_figures.py
"""
from __future__ import annotations

import gzip
import importlib
import json
import sys
from pathlib import Path

DISC = Path(__file__).resolve().parent
sys.path.insert(0, str(DISC))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from common import OUT, REGIONS, con, load_tracts, log, p, region_outline, to5070, wkb_to_geoms, xy_to4326, \
    xy_to5070  # noqa: E402

eg = importlib.import_module("02_egress")
ox = importlib.import_module("07_external_osm")

FIGS = DISC / "figs"
SURF, INK, INK2, MUTED, GRID, BASE = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, BLUE_L, ORANGE, AQUA = "#2a78d6", "#86b6ef", "#eb6834", "#1baf7a"
OTHER_BLDG = "#dcdbd4"
SEED, B = 20261012, 2000
plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9, "axes.edgecolor": BASE,
                     "axes.labelcolor": INK2, "xtick.color": MUTED, "ytick.color": MUTED, "text.color": INK})

MAPS = [  # (label in 07 NAMED, region, file stem, title, place)
    ("Rancho Tehama", "northern-ca", "pocket_rancho_tehama", "Rancho Tehama, Tehama County, CA"),
    ("Lake California", "northern-ca", "pocket_lake_california", "Lake California, Cottonwood, CA"),
    ("Whitehorn Cove (Wagoner County)", "eastern-ok", "pocket_whitehorn_cove", "Whitehorn Cove, Wagoner County, OK"),
]


def lines_xy(geoms):
    out = []
    for g in geoms:
        for part in shapely.get_parts(g):
            co = shapely.get_coordinates(part)
            if len(co) >= 2:
                out.append(co)
    return out


def nice_scale(width_m):
    for km in (0.5, 1, 2, 5, 10, 20):
        if km * 1000 >= width_m * 0.15:
            return km
    return 20


def pocket_map(label, region, stem, title, G, tr, tg, c, ev_df):
    named = {n[0]: n for n in ox.NAMED}
    _, _, geoid, road = named[label]
    P = pd.read_parquet(OUT / f"{region}_pockets.parquet")
    P = P[(P.GEOID == geoid) & ~P.exit_outside_region]
    if road is not None:
        P = P[P.egress_road == road]
    row = P.sort_values("n_bldg", ascending=False).iloc[0]
    pid = int(row.pocket)
    pk, eu, evv = G["pocket"], G["eu"], G["ev"]
    real = np.arange(len(eu))
    in_p = (pk[eu] == pid) & (pk[evv] == pid)
    att = ((pk[eu] == pid) & G["core"][evv]) | ((pk[evv] == pid) & G["core"][eu])
    pg = ox.edge_geoms(G, real[in_p])
    ag = ox.edge_geoms(G, real[att])
    x0, y0, x1, y1 = shapely.total_bounds(np.r_[pg, ag])
    pad = max(700.0, 0.10 * max(x1 - x0, y1 - y0))
    X0, Y0, X1, Y1 = x0 - pad, y0 - pad, x1 + pad, y1 + pad
    box = shapely.box(X0, Y0, X1, Y1)
    segs = G["seg_tree"].query(box, predicate="intersects")
    seg_mask = np.zeros(G["n"], bool)
    seg_mask[segs] = True
    win = real[seg_mask[G["eseg"]]]
    core_e = win[G["core"][eu[win]] & G["core"][evv[win]]]
    other_e = win[~(G["core"][eu[win]] & G["core"][evv[win]]) & ~in_p[win] & ~att[win]]
    cg = ox.edge_geoms(G, core_e)
    og = ox.edge_geoms(G, other_e)
    # buildings in the window (footprints), assigned exactly as in the pipeline
    lo, la = xy_to4326(np.array([X0, X1, X0, X1]), np.array([Y0, Y0, Y1, Y1]))
    w_, e_, s_, n_ = lo.min(), lo.max(), la.min(), la.max()
    bq = (f"SELECT ST_AsWKB(geometry) w, (bbox.xmin+bbox.xmax)/2 cx, (bbox.ymin+bbox.ymax)/2 cy "
          f"FROM '{p(region, 'buildings')}' WHERE bbox.xmin <= {e_} AND bbox.xmax >= {w_} "
          f"AND bbox.ymin <= {n_} AND bbox.ymax >= {s_}")
    bd = c.execute(bq).fetchdf()
    bgeo = to5070(wkb_to_geoms(bd.w))
    bx, by = xy_to5070(bd.cx.values, bd.cy.values)
    (qq, tt), _ = G["seg_tree"].query_nearest(shapely.points(bx, by), max_distance=150.0, return_distance=True,
                                              all_matches=False)
    bseg = np.full(len(bx), -1, np.int64)
    bseg[qq] = tt
    batt = eg.attach_buildings(G, bx, by, bseg)
    pu = np.where(batt >= 0, pk[eu[np.maximum(batt, 0)]], -1)
    pv = np.where(batt >= 0, pk[evv[np.maximum(batt, 0)]], -1)
    inpocket = (pu == pid) | (pv == pid)

    def polys(gs):
        out = []
        for g in gs:
            for part in shapely.get_parts(g):
                if part.geom_type == "Polygon":
                    out.append(shapely.get_coordinates(part.exterior))
        return out

    # OSM-only second link (07_external_osm.py)
    er = ev_df[(ev_df.label == label)]
    osm_lines, gates, osm_cls, osm_note = [], [], None, ""
    if len(er):
        er = er.iloc[0]
        osm_cls = er.osm_class
        if osm_cls in ("b", "c") and isinstance(er.link_osm_way_ids, str):
            raw = ox.RAW / f"{region}_p{pid}_m{int(er.margin_m)}.json.gz"
            with gzip.open(raw, "rt", encoding="utf-8") as fh:
                j = json.load(fh)
            nodes = {el["id"]: (el["lon"], el["lat"]) for el in j["elements"] if el["type"] == "node" and "lat" in el}
            bar = {el["id"] for el in j["elements"] if el["type"] == "node" and
                   el.get("tags", {}).get("barrier") in ox.BARRIERS}
            ids = er.link_path_osm_way_ids if isinstance(er.get("link_path_osm_way_ids"), str) else er.link_osm_way_ids
            want = {int(v) for v in str(ids).split(";") if v}
            for el in j["elements"]:
                if el["type"] == "way" and el["id"] in want:
                    nd = [nodes[k] for k in el["nodes"] if k in nodes]
                    if len(nd) >= 2:
                        xs, ys = xy_to5070(np.array([q[0] for q in nd]), np.array([q[1] for q in nd]))
                        osm_lines.append(np.c_[xs, ys])
            if isinstance(er.get("link_barriers_lonlat"), str) and er.link_barriers_lonlat:
                for pair in er.link_barriers_lonlat.split(";"):
                    blon, blat = (float(v) for v in pair.split(","))
                    gx, gy = xy_to5070(np.array([blon]), np.array([blat]))
                    gates.append((float(gx[0]), float(gy[0])))
        osm_note = {"a": "OSM check (2026-10-13): no second connection, single access confirmed",
                    "b": "OSM check (2026-10-13): a public second link exists in OSM",
                    "c": "OSM check (2026-10-13): second link only via a gated / private / track road"}.get(osm_cls, "")
    # tract outline
    ti = int(np.where(tr.GEOID.values == geoid)[0][0])
    t_line = shapely.intersection(shapely.boundary(tg[ti]), box)

    W, H = X1 - X0, Y1 - Y0
    fw = 8.0
    fh = min(max(fw * H / W, 5.0), 9.5) + 1.3
    fig, ax = plt.subplots(figsize=(fw, fh), dpi=150)
    fig.patch.set_facecolor(SURF)
    ax.set_facecolor(SURF)
    ax.add_collection(PolyCollection(polys(bgeo[~inpocket]), facecolors=OTHER_BLDG, edgecolors="none", zorder=1))
    ax.add_collection(PolyCollection(polys(bgeo[inpocket]), facecolors=BLUE_L, edgecolors="none", zorder=2))
    ax.add_collection(LineCollection(lines_xy(og), colors=BASE, linewidths=0.6, zorder=3))
    ax.add_collection(LineCollection(lines_xy(cg), colors=MUTED, linewidths=0.9, zorder=4))
    ax.add_collection(LineCollection(lines_xy(pg), colors=BLUE, linewidths=1.3, zorder=5))
    if osm_lines:
        ax.add_collection(LineCollection(osm_lines, colors=AQUA, linewidths=2.2, zorder=6))
    ax.add_collection(LineCollection(lines_xy(ag), colors=ORANGE, linewidths=3.2, capstyle="round", zorder=7))
    if gates:
        gx, gy = zip(*gates)
        ax.scatter(gx, gy, s=46, marker="X", color=INK, edgecolors=SURF, linewidths=1.2, zorder=8)
    for part in shapely.get_parts(t_line):
        co = shapely.get_coordinates(part)
        if len(co) >= 2:
            ax.plot(co[:, 0], co[:, 1], color=INK2, lw=0.8, ls=(0, (4, 3)), zorder=2.5)
    # direct label on the single link
    am = None
    if len(ag):
        alink = max(shapely.get_parts(shapely.line_merge(shapely.union_all(ag))), key=lambda g: g.length)
        am = (shapely.line_interpolate_point(alink, 0.5, normalized=True) if alink.geom_type == "LineString"
              else shapely.centroid(alink))
    if am is not None and not am.is_empty:
        ax.scatter([am.x], [am.y], s=70, color=ORANGE, edgecolors=SURF, linewidths=2, zorder=8)
        right = am.x > X0 + 0.6 * W
        top = am.y > Y0 + 0.7 * H
        ax.annotate(f"single link: {row.egress_road}", xy=(am.x, am.y),
                    xytext=(-16 if right else 16, -22 if top else 18), textcoords="offset points",
                    ha="right" if right else "left", fontsize=8.5, color=INK, zorder=9,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8),
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURF, ec="none", alpha=0.9))
    if osm_lines:
        mid = osm_lines[len(osm_lines) // 2]
        mp = mid[len(mid) // 2]
        ax.annotate("OSM-only second link\n(gated track)" if osm_cls == "c" else "OSM-only second link",
                    xy=(mp[0], mp[1]), xytext=(16, -26), textcoords="offset points", fontsize=8.5, color=INK,
                    zorder=9, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8),
                    bbox=dict(boxstyle="round,pad=0.25", fc=SURF, ec="none", alpha=0.9))
    # scale bar
    km = nice_scale(W)
    sx, sy = X1 - 0.04 * W - km * 1000, Y0 + 0.05 * H
    ax.add_patch(matplotlib.patches.Rectangle((sx - 0.015 * W, sy - 0.025 * H), km * 1000 + 0.03 * W, 0.075 * H,
                                              facecolor=SURF, edgecolor="none", alpha=0.9, zorder=9.5))
    ax.plot([sx, sx + km * 1000], [sy, sy], color=INK, lw=2.5, solid_capstyle="butt", zorder=10)
    ax.text(sx + km * 500, sy + 0.015 * H, f"{km:g} km", ha="center", va="bottom", fontsize=8, color=INK, zorder=10)
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_aspect("equal")
    ax.set_axis_off()
    handles = [Line2D([], [], color=BLUE, lw=1.6, label="roads inside the pocket"),
               Patch(color=BLUE_L, label=f"buildings inside the pocket ({int(inpocket.sum()):,} in view)"),
               Line2D([], [], color=ORANGE, lw=3.2, label="the single link to the meshed network"),
               Line2D([], [], color=MUTED, lw=1.0, label="meshed (2-edge-connected) network"),
               Line2D([], [], color=BASE, lw=1.0, label="other dead-end roads"),
               ]
    if not t_line.is_empty:
        handles.append(Line2D([], [], color=INK2, lw=0.8, ls=(0, (4, 3)), label=f"tract {geoid} boundary"))
    if osm_lines:
        handles.insert(3, Line2D([], [], color=AQUA, lw=2.2, label="OSM-only second link (not in Overture drivable)"))
    if gates:
        handles.insert(4, Line2D([], [], marker="X", color=INK, ls="none", markersize=7, label="gate / barrier (OSM)"))
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.01), ncol=2, frameon=False,
                    fontsize=8, handlelength=2.2, columnspacing=1.5)
    for t in leg.get_texts():
        t.set_color(INK2)
    fig.text(0.02, 0.985, f"{title}: {int(row.n_bldg):,} buildings behind one road segment", fontsize=12,
             fontweight="semibold", color=INK, va="top")
    sub = (f"Overture drivable graph (challenge extract). Tract {geoid}"
           f"{' (boundary outside view)' if t_line.is_empty else ''}; est. residents {row.est_pop:,.0f}; "
           f"tract WHP {row.usfs_WHP_mean:,.0f}; {row.road_km:.0f} km of road inside the pocket.")
    fig.text(0.02, 0.952, sub, fontsize=8.5, color=INK2, va="top")
    if osm_note:
        fig.text(0.02, 0.928, osm_note + ".", fontsize=8.5, color=INK2, va="top")
    fig.text(0.98, 0.006, "Overture Maps data via challenge extract; OSM © OpenStreetMap contributors (ODbL 1.0); "
                          "EPSG:5070", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.17)
    out = FIGS / f"{stem}.png"
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)
    log("wrote", out)
    return out


def summary_chart():
    d = pd.read_parquet(OUT / "tract_metrics_all.parquet")
    rng = np.random.default_rng(SEED)
    panels = [("All five regions (pooled quartiles)", d, "whp_q")] + \
             [(f"{r} (within-region quartiles)", d[d.region == r], "whp_q_region") for r in REGIONS]
    rows = []
    for name, sub0, col in panels:
        for q in ["Q1", "Q2", "Q3", "Q4"]:
            sub = sub0[sub0[col].astype(str) == q]
            x, y = sub.est_pop_pocket50.values, sub.pop_total.values
            idx = rng.integers(0, len(sub), size=(B, len(sub)))
            reps = x[idx].sum(1) / y[idx].sum(1)
            rows.append(dict(panel=name, quartile=q, n_tracts=len(sub), rate=x.sum() / y.sum(),
                             ci_lo=np.percentile(reps, 2.5), ci_hi=np.percentile(reps, 97.5)))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "fig_whp_by_region.csv", index=False)
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.2), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURF)
    ymax = float(T.ci_hi.max()) * 100 * 1.15
    xs = np.arange(4)
    for ax, (name, _, _) in zip(axes.ravel(), panels):
        t = T[T.panel == name]
        r = t.rate.values * 100
        lo, hi = t.ci_lo.values * 100, t.ci_hi.values * 100
        ax.set_facecolor(SURF)
        ax.vlines(xs, lo, hi, color=BLUE, lw=2, zorder=2)
        ax.scatter(xs, r, s=48, color=BLUE, edgecolors=SURF, linewidths=2, zorder=3)
        for k in (0, 3):
            ax.text(xs[k] + 0.12, r[k], f"{r[k]:.1f}%", fontsize=8, color=INK, va="center", ha="left")
        ax.set_xticks(xs, ["Q1\nlowest", "Q2", "Q3", "Q4\nhighest"])
        ax.set_xlim(-0.45, 3.65)
        ax.set_ylim(0, ymax)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=0)
        ratio = r[3] / r[0]
        ax.set_title(f"{name.split(' (')[0]}\nQ4 / Q1 = {ratio:.1f}x", fontsize=9, color=INK, loc="left")
    for ax in axes[:, 0]:
        ax.set_ylabel("% of residents in a single-egress pocket")
    fig.text(0.01, 0.985, "Residents whose only road link to the meshed network is one segment, "
                          "by wildfire-hazard quartile", fontsize=12, fontweight="semibold", color=INK, va="top")
    fig.text(0.01, 0.945, "Pockets of at least 50 buildings in the Overture drivable graph. Dots: population-weighted "
                          "share; bars: 95% tract-bootstrap CI. Quartiles of tract-mean USFS WHP (pooled for the first "
                          "panel, within each region otherwise).", fontsize=8, color=INK2, va="top", wrap=True)
    fig.subplots_adjust(left=0.07, right=0.99, top=0.82, bottom=0.08, hspace=0.55, wspace=0.12)
    out = FIGS / "whp_quartile_by_region.png"
    fig.savefig(out, dpi=150, facecolor=SURF)
    plt.close(fig)
    log("wrote", out)
    return out


def main():
    FIGS.mkdir(exist_ok=True)
    ev_df = pd.read_csv(OUT / "external_validation.csv", dtype={"GEOID": str})
    for region in ["northern-ca", "eastern-ok"]:
        maps = [m for m in MAPS if m[1] == region]
        c = con()
        tr, tg = load_tracts(region)
        region_poly, pieces = region_outline(tg)
        G = eg.build_graph(c, region, "roads", None, "d", region_poly, pieces)
        G["seg_tree"] = shapely.STRtree(G["geoms"])
        for label, _, stem, title in maps:
            pocket_map(label, region, stem, title, G, tr, tg, c, ev_df)
        del G
    summary_chart()


if __name__ == "__main__":
    main()
