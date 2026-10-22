"""External validation of single-egress pockets against current OpenStreetMap (public Overpass API).

Question per pocket: does OSM have a second road connection between the pocket and the meshed network that
Overture's drivable graph lacks?

Sample (seeded): the 16 named communities of FINDINGS.md section 5, plus a simple random sample of 40 pockets
(>= 50 buildings, verifiable exit) from ALL pockets in the highest pooled wildfire-hazard quartile (WHP Q4,
frame = 5,972 pockets), drawn with numpy default_rng(20261013).

Method per pocket (the graph is rebuilt exactly as in 02_egress.py, so pocket ids match out/<region>_pockets):
  1. Window = bounding box of the pocket's Overture edges (+ its single link) expanded by MARGIN_M.
  2. One Overpass query: all highway ways of vehicle classes in the window (incl. service and track), their
     nodes (skeleton) and barrier nodes. Response cached gzip-compressed in discovery/external/osm_raw/ together
     with a .meta.json (endpoint, query text, UTC retrieval time, osm_base timestamp, bytes).
  3. OSM graph from shared node ids. Pocket nodes: OSM nodes <= 15 m from the pocket's Overture edges and > 40 m
     from core edges and > 15 m from the known link. Core nodes: symmetric (<= 15 m from core edges).
  4. Unit-capacity max-flow (edge-disjoint paths) from pocket nodes to core nodes:
       public graph = excludes highway=track, access/motor_vehicle/vehicle/motorcar in RESTRICTED,
                      service=driveway/parking_aisle/drive-through/emergency_access, and edges at barrier nodes;
       full graph   = all vehicle ways.
     (b) public flow >= 2 : a second public link exists in OSM that the Overture drivable graph lacks;
         sub-type from the OSM ways of that link that are > 15 m from any Overture drivable segment:
         'road' (a road class), 'service' (public service road), 'topology' (no such ways: OSM has a junction
         that Overture's connector graph lacks); share of that gap covered by Overture -roads-unfiltered.
     (c) full flow >= 2 but public flow < 2 : second link only via private / gated / track / driveway.
     (a) full flow == 1 : single access in OSM as well.
     'unresolved' : no OSM path found inside the window, or a second path only leaves the window; such
     pockets are re-queried once with MARGIN2_M.
     A second path that runs along the known link (>= 80% of its length within 40 m of it) is treated as the
     same corridor (e.g. a dual carriageway), not as an independent link.

OSM data (c) OpenStreetMap contributors, ODbL 1.0 (https://www.openstreetmap.org/copyright).
Usage: python -I discovery/07_external_osm.py
Outputs: discovery/out/external_validation.csv, discovery/out/external_validation_summary.csv
"""
from __future__ import annotations

import datetime as dt
import gzip
import importlib
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DISC = Path(__file__).resolve().parent
sys.path.insert(0, str(DISC))  # run with python -I; downloads live in external/osm_raw, never imported

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402
from scipy.sparse import csr_matrix  # noqa: E402
from scipy.sparse.csgraph import maximum_flow  # noqa: E402
from shapely.ops import substring  # noqa: E402

from common import (DRIVABLE, OUT, REGIONS, con, load_tracts, log, p, region_outline, to5070,  # noqa: E402
                    wkb_to_geoms, xy_to4326, xy_to5070)

eg = importlib.import_module("02_egress")

MAIN = "https://overpass-api.de/api/interpreter"  # main public instance (tried on every even attempt)
MIRRORS = ["https://overpass.kumi.systems/api/interpreter",  # public mirrors, tried on odd attempts
           "https://overpass.private.coffee/api/interpreter"]
USER_AGENT = "bias-bounty-mapping-equity-discovery/1.0 (research script; small bbox queries)"
RAW = DISC / "external" / "osm_raw"
MARGIN_M = 750.0
MARGIN2_M = 3000.0
BUDGET_WIRE_BYTES = 25_000_000
RANDOM_BUDGET_S = 60 * 60  # wall-clock budget for downloading the random sample (named pockets always tried)
SEED = 20261013
N_RANDOM = 40
HW = ("motorway|trunk|primary|secondary|tertiary|unclassified|residential|living_street|service|track|road|"
      "motorway_link|trunk_link|primary_link|secondary_link|tertiary_link|busway")
QUERY = """[out:json][timeout:120];
way["highway"~"^({hw})$"]({s:.6f},{w:.6f},{n:.6f},{e:.6f})->.w;
.w out body qt;
node(w.w)->.n;
.n out skel qt;
node.n["barrier"];
out body qt;"""
RESTRICTED = {"private", "no", "agricultural", "forestry", "emergency", "permit", "military", "customers",
              "delivery"}
INFORMAL_SERVICE = {"driveway", "parking_aisle", "drive-through", "emergency_access"}
BARRIERS = {"gate", "lift_gate", "swing_gate", "sliding_gate", "chain", "bollard", "block", "jersey_barrier",
            "log", "kissing_gate", "turnstile", "full-height_turnstile", "barrier_board", "debris"}
NEAR_M, FAR_M, LINK_M = 15.0, 40.0, 15.0
PARALLEL_M, PARALLEL_SHARE = 40.0, 0.8
BIG = 1_000_000

NAMED = [  # label, region, GEOID, single-link road name (None = largest pocket of the tract)
    ("Rancho Tehama", "northern-ca", "06103000300", "Rancho Tehama Road"),
    ("Lake California", "northern-ca", "06103000201", "Lake California Drive"),
    ("Olympic Valley", "northern-ca", "06061022014", "Olympic Valley Road"),
    ("Westwood Avenue (Redding)", "northern-ca", "06089011002", "Westwood Avenue"),
    ("Icicle Road", "eastern-wa", "53007960201", "Icicle Road"),
    ("Stirrup Road (Leavenworth)", "eastern-wa", "53007960202", "Stirrup Road"),
    ("Casassa Road (Cle Elum)", "eastern-wa", "53037975201", "Casassa Road"),
    ("Salmon la Sac Road", "eastern-wa", "53037975101", "Salmon la Sac Road"),
    ("Seven Bays (Davenport)", "eastern-wa", "53043960200", "Seven Bays Road East"),
    ("Crescent Bar (Quincy)", "eastern-wa", "53025010500", "Crescent Bar Road Northwest"),
    ("Whitehorn Cove (Wagoner County)", "eastern-ok", "40145030201", "Whitehorn Cove Road"),
    ("Basin Road (Mannford)", "eastern-ok", "40037020706", "Basin Road"),
    ("Rimrock / Lake Montezuma", "maricopa-az", "04025001605", "North Montezuma Avenue"),
    ("Bagdad", "maricopa-az", "04025001300", "Hillside-Bagdad Highway"),
    ("Sargent (FM 457)", "south-central-tx", "48321730503", "Farm-to-Market Road 457"),
    ("Port O'Connor", "south-central-tx", "48057000501", None),
]


# ----------------------------------------------------------------------------------------------- selection
def select_pockets() -> pd.DataFrame:
    tm = pd.read_parquet(OUT / "tract_metrics_all.parquet")[["GEOID", "whp_q"]]
    P = pd.concat([pd.read_parquet(OUT / f"{r}_pockets.parquet") for r in REGIONS], ignore_index=True)
    P = P[(P.n_bldg >= 50) & ~P.exit_outside_region].merge(tm, on="GEOID", how="left")
    rows = []
    for label, region, geoid, road in NAMED:
        s = P[(P.region == region) & (P.GEOID == geoid)]
        if road is not None:
            s = s[s.egress_road == road]
        s = s.sort_values("n_bldg", ascending=False).head(1)
        assert len(s) == 1, (label, region, geoid, road)
        rows.append(s.assign(sample="named", label=label))
    named = pd.concat(rows, ignore_index=True)
    frame = P[(P.whp_q.astype(str) == "Q4")]
    key = frame.region + ":" + frame.pocket.astype(str)
    frame = frame[~key.isin(named.region + ":" + named.pocket.astype(str))]
    rng = np.random.default_rng(SEED)
    pick = rng.choice(len(frame), size=N_RANDOM, replace=False)
    rand = frame.iloc[np.sort(pick)].assign(sample="random_q4", label="")
    log("sample frame (WHP-Q4 pockets >= 50 bldg, excluding named):", len(frame))
    return pd.concat([named, rand], ignore_index=True)


# ----------------------------------------------------------------------------------------------- overpass
_spent = 0


class NotAttempted(RuntimeError):
    pass


FETCH = dict(attempts=4, allow_download=True)


def fetch(key: str, s: float, w: float, n: float, e: float, margin: float) -> tuple[dict, dict]:
    """Cached Overpass call; on repeated failure the window is split into 2 x 2 tiles. Returns (json, meta)."""
    try:
        return fetch_one(key, s, w, n, e, margin, FETCH["attempts"])
    except OverpassUnavailable:
        if not FETCH["allow_download"] or not FETCH.get("tiles", False):
            raise
        log(f"  {key}: splitting window into 2 x 2 tiles")
        ms, mw = (s + n) / 2, (w + e) / 2
        parts, metas = [], []
        for ti, (a, b, c_, d) in enumerate([(s, w, ms, mw), (s, mw, ms, e), (ms, w, n, mw), (ms, mw, n, e)]):
            j, m = fetch_one(f"{key}_t{ti}", a, b, c_, d, margin, FETCH["attempts"])
            parts.append(j)
            metas.append(m)
        seen, elems = set(), []
        for j in parts:
            for el in j["elements"]:
                k = (el["type"], el["id"], "tags" in el)
                if k not in seen:
                    seen.add(k)
                    elems.append(el)
        meta = dict(metas[0])
        meta.update(query="(2 x 2 tiles) " + " || ".join(m["query"] for m in metas), bbox=[s, w, n, e],
                    wire_bytes=sum(m["wire_bytes"] for m in metas), json_bytes=sum(m["json_bytes"] for m in metas),
                    endpoint=";".join(sorted({m["endpoint"] for m in metas})),
                    retrieved_utc=max(m["retrieved_utc"] for m in metas))
        return dict(elements=elems, osm3s=parts[0].get("osm3s", {})), meta


def fetch_one(key: str, s: float, w: float, n: float, e: float, margin: float, attempts: int):
    global _spent
    raw_f = RAW / f"{key}_m{int(margin)}.json.gz"
    meta_f = RAW / f"{key}_m{int(margin)}.meta.json"
    if raw_f.exists() and meta_f.exists():
        with gzip.open(raw_f, "rt", encoding="utf-8") as fh:
            return json.load(fh), json.loads(meta_f.read_text(encoding="utf-8"))
    if not FETCH["allow_download"]:
        raise NotAttempted("time budget for downloads exhausted")
    if _spent > BUDGET_WIRE_BYTES:
        raise NotAttempted("download byte budget exhausted")
    q = QUERY.format(hw=HW, s=s, w=w, n=n, e=e)
    data = urllib.parse.urlencode({"data": q}).encode()
    last = None
    for attempt in range(attempts):
        ep = [MAIN, MAIN, MIRRORS[0], MAIN, MIRRORS[1], MAIN][attempt % 6]
        if attempt >= 1:
            time.sleep(min(5 * 2 ** attempt, 120))  # back off: 10, 20, 40, 80, 120 s
        req = urllib.request.Request(ep, data=data, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"})
        try:
            with urllib.request.urlopen(req, timeout=120 if ep == MAIN else 60) as r:
                wire = r.read()
                enc = r.headers.get("Content-Encoding")
                status = r.status
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as ex:
            last = ex
            log(f"  overpass error ({key}, attempt {attempt}, {ep}): {ex}")
            continue
        body = gzip.decompress(wire) if enc == "gzip" else wire
        j = json.loads(body.decode("utf-8"))
        if "remark" in j and "error" in str(j["remark"]).lower():
            last = RuntimeError(j["remark"])
            log(f"  overpass remark ({key}): {j['remark']}")
            continue
        _spent += len(wire)
        meta = dict(endpoint=ep, query=q, retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
                    osm_base=j.get("osm3s", {}).get("timestamp_osm_base"), http_status=status,
                    wire_bytes=len(wire), json_bytes=len(body), bbox=[s, w, n, e], margin_m=margin)
        with gzip.open(raw_f, "wt", encoding="utf-8") as fh:
            json.dump(j, fh)
        meta_f.write_text(json.dumps(meta, indent=1), encoding="utf-8")
        log(f"  fetched {key}: {meta['wire_bytes']:,} B from {ep}, osm_base {meta['osm_base']}")
        time.sleep(2.0)  # be polite to the public instance: strictly sequential, pause between calls
        return j, meta
    raise OverpassUnavailable(f"overpass failed for {key}: {last}")


class OverpassUnavailable(RuntimeError):
    pass


# ----------------------------------------------------------------------------------------------- geometry
def edge_geoms(G, idx):
    L = G["length"][G["eseg"][idx]]
    a0 = G["ea0"][idx]
    a1 = a0 + np.where(L > 0, G["elen"][idx] / np.maximum(L, 1e-9), 0.0)
    out = []
    for k, i in enumerate(idx):
        g = G["geoms"][G["eseg"][i]]
        if g.geom_type != "LineString":  # clipped segments can be multi-part: merge, else keep whole segment
            gm = shapely.line_merge(g)
            if gm.geom_type != "LineString":
                out.append(g)
                continue
            g = gm
        lo, hi = float(min(a0[k], a1[k])), float(max(a0[k], a1[k]))
        sub = substring(g, lo, hi, normalized=True)
        out.append(sub if not sub.is_empty else g)
    return np.array(out, dtype=object)


def osm_graph(j):
    nodes = {e["id"]: (e["lon"], e["lat"]) for e in j["elements"] if e["type"] == "node" and "lat" in e}
    barrier_nodes = {e["id"] for e in j["elements"] if e["type"] == "node" and
                     e.get("tags", {}).get("barrier") in BARRIERS}
    ways = [e for e in j["elements"] if e["type"] == "way"]
    ids = np.array(sorted(nodes))
    pos = {nid: k for k, nid in enumerate(ids)}
    lon = np.array([nodes[i][0] for i in ids])
    lat = np.array([nodes[i][1] for i in ids])
    x, y = xy_to5070(lon, lat)
    eu, ev, ew = [], [], []
    winfo = []
    for wi, w in enumerate(ways):
        t = w.get("tags", {})
        hw = t.get("highway", "")
        acc = {t.get(k, "") for k in ("access", "motor_vehicle", "motorcar", "vehicle")}
        restricted = bool(acc & RESTRICTED)
        informal = hw == "service" and t.get("service", "") in INFORMAL_SERVICE
        winfo.append(dict(way_id=w["id"], highway=hw, service=t.get("service", ""), restricted=restricted,
                          informal=informal, track=hw == "track", name=t.get("name", "")))
        nd = [pos[n] for n in w.get("nodes", []) if n in pos]
        for a, b in zip(nd[:-1], nd[1:]):
            if a != b:
                eu.append(a)
                ev.append(b)
                ew.append(wi)
    eu, ev, ew = np.array(eu, np.int64), np.array(ev, np.int64), np.array(ew, np.int64)
    is_barrier = np.isin(ids, list(barrier_nodes)) if barrier_nodes else np.zeros(len(ids), bool)
    W = pd.DataFrame(winfo)
    return dict(ids=ids, x=x, y=y, lon=lon, lat=lat, eu=eu, ev=ev, ew=ew, W=W, barrier=is_barrier)


def max_flow_paths(n, eu, ev, src, dst):
    """Unit-capacity undirected max-flow between node sets; returns (flow, list of paths as edge-index lists)."""
    if len(src) == 0 or len(dst) == 0 or len(eu) == 0:
        return 0, []
    S, T = n, n + 1
    m = len(eu)
    rows = np.r_[eu, ev, np.full(len(src), S), dst]
    cols = np.r_[ev, eu, src, np.full(len(dst), T)]
    caps = np.r_[np.ones(2 * m, np.int64), np.full(len(src), BIG), np.full(len(dst), BIG)]
    M = csr_matrix((caps.astype(np.int32), (rows, cols)), shape=(n + 2, n + 2))
    M.sum_duplicates()
    M.setdiag(0)
    M.eliminate_zeros()
    res = maximum_flow(M, S, T)
    f = int(res.flow_value)
    F = res.flow.tocoo()
    pos = F.data > 0
    fu, fv, fd = F.row[pos], F.col[pos], F.data[pos]
    out = {}
    for a, b, c in zip(fu, fv, fd):
        out.setdefault(a, []).append([b, int(c)])
    # map (a,b) -> an edge index (either orientation)
    emap = {}
    for k in range(m):
        emap.setdefault((eu[k], ev[k]), k)
        emap.setdefault((ev[k], eu[k]), k)
    paths = []
    for _ in range(min(f, 4)):
        path, cur, seen = [], S, set()
        while cur != T and cur in out and len(path) < 100000:
            nxt = None
            for item in out[cur]:
                if item[1] > 0:
                    nxt = item
                    break
            if nxt is None:
                break
            nxt[1] -= 1
            b = nxt[0]
            if cur not in (S,) and b != T:
                k = emap.get((cur, b))
                if k is not None:
                    path.append(k)
            if (cur, b) in seen:
                break
            seen.add((cur, b))
            cur = b
        paths.append(path)
    return f, paths


# ----------------------------------------------------------------------------------------------- per pocket
def classify(row, G, region, c, margin):
    pid = int(row.pocket)
    pk = G["pocket"]
    real = np.arange(len(G["eu"]))
    eu_, ev_ = G["eu"], G["ev"]
    in_p = (pk[eu_] == pid) & (pk[ev_] == pid)
    att = ((pk[eu_] == pid) & G["core"][ev_]) | ((pk[ev_] == pid) & G["core"][eu_])
    n_nodes = int((pk[:-1] == pid).sum())
    assert n_nodes == int(row.n_nodes), ("pocket id mismatch", region, pid, n_nodes, row.n_nodes)
    pe = real[in_p]
    ae = real[att]
    pg = edge_geoms(G, pe) if len(pe) else np.array([], dtype=object)
    ag = edge_geoms(G, ae)
    allg = np.r_[pg, ag]
    xmin, ymin, xmax, ymax = shapely.total_bounds(allg)
    X0, Y0, X1, Y1 = xmin - margin, ymin - margin, xmax + margin, ymax + margin
    lo, la = xy_to4326(np.array([X0, X1, X0, X1]), np.array([Y0, Y0, Y1, Y1]))
    s, w, n, e = la.min(), lo.min(), la.max(), lo.max()
    key = f"{region}_p{pid}"
    j, meta = fetch(key, s, w, n, e, margin)
    O = osm_graph(j)
    # Overture core edges in the window
    box = shapely.box(X0, Y0, X1, Y1)
    segs = G["seg_tree"].query(box, predicate="intersects")
    seg_mask = np.zeros(G["n"], bool)
    seg_mask[segs] = True
    ce = real[seg_mask[G["eseg"]] & G["core"][eu_] & G["core"][ev_]]
    cg = edge_geoms(G, ce) if len(ce) else np.array([], dtype=object)
    pts = shapely.points(O["x"], O["y"])

    def dist_to(geoms):
        if len(geoms) == 0:
            return np.full(len(pts), np.inf)
        t = shapely.STRtree(geoms)
        (qi, _), d = t.query_nearest(pts, return_distance=True, all_matches=False)
        out = np.full(len(pts), np.inf)
        out[qi] = d
        return out

    d_p, d_c, d_a = dist_to(pg), dist_to(cg), dist_to(ag)
    pocket_nodes = np.where((d_p <= NEAR_M) & (d_c > FAR_M) & (d_a > LINK_M))[0]
    core_nodes = np.where((d_c <= NEAR_M) & (d_p > FAR_M) & (d_a > LINK_M))[0]
    inside = (O["lon"] >= w) & (O["lon"] <= e) & (O["lat"] >= s) & (O["lat"] <= n)
    out_nodes = np.where(~inside)[0]
    Wd = O["W"]
    ew = O["ew"]
    e_restr = Wd.restricted.values[ew] if len(ew) else np.zeros(0, bool)
    e_track = Wd.track.values[ew] if len(ew) else np.zeros(0, bool)
    e_inf = Wd.informal.values[ew] if len(ew) else np.zeros(0, bool)
    e_bar = O["barrier"][O["eu"]] | O["barrier"][O["ev"]] if len(ew) else np.zeros(0, bool)
    public = ~(e_restr | e_track | e_inf | e_bar)
    N = len(O["ids"])
    f_pub, paths_pub = max_flow_paths(N, O["eu"][public], O["ev"][public], pocket_nodes, core_nodes)
    pub_idx = np.where(public)[0]
    paths_pub = [[pub_idx[k] for k in pth] for pth in paths_pub]
    f_full, paths_full = max_flow_paths(N, O["eu"], O["ev"], pocket_nodes, core_nodes)
    f_win, _ = max_flow_paths(N, O["eu"], O["ev"], pocket_nodes, np.r_[core_nodes, out_nodes])

    # Overture geometry for characterising a second link: the region-wide drivable tree (OSM ways can run
    # outside the query window), and the -roads-unfiltered non-drivable segments around each link's extent
    drv_tree = G["seg_tree"]
    drv_list = ",".join("'" + k + "'" for k in DRIVABLE)
    unf_src = p(region, "roads-unfiltered")

    def unf_tree_for(lo0, la0, lo1, la1):
        q = (f"SELECT ST_AsWKB(geometry) w FROM '{unf_src}' WHERE class NOT IN ({drv_list}) "
             f"AND bbox.xmax >= {lo0} AND bbox.xmin <= {lo1} AND bbox.ymax >= {la0} AND bbox.ymin <= {la1}")
        u = c.execute(q).fetchdf()
        return shapely.STRtree(to5070(wkb_to_geoms(u.w))) if len(u) else None

    ag_union = shapely.union_all(ag)

    def describe(path):
        if not path:
            return None
        k = np.array(path)
        a, b = O["eu"][k], O["ev"][k]
        segs_g = shapely.linestrings(np.stack([np.c_[O["x"][a], O["y"][a]], np.c_[O["x"][b], O["y"][b]]], axis=1))
        lens = shapely.length(segs_g)
        mids = shapely.line_interpolate_point(segs_g, 0.5, normalized=True)
        near_link = shapely.distance(mids, ag_union) <= PARALLEL_M
        par_share = float(lens[near_link].sum() / max(lens.sum(), 1e-9))
        (qi, _), dd = drv_tree.query_nearest(mids, return_distance=True, all_matches=False)
        ddrv = np.full(len(mids), np.inf)
        ddrv[qi] = dd
        gap = ddrv > NEAR_M
        unf_tree = None
        if gap.any():
            gco = shapely.get_coordinates(segs_g[gap])
            glo, gla = xy_to4326(gco[:, 0], gco[:, 1])
            unf_tree = unf_tree_for(glo.min() - 0.002, gla.min() - 0.002, glo.max() + 0.002, gla.max() + 0.002)
        if unf_tree is not None and gap.any():
            (qi, _), dd = unf_tree.query_nearest(mids[gap], return_distance=True, all_matches=False)
            du = np.full(int(gap.sum()), np.inf)
            du[qi] = dd
            in_unf = float(lens[gap][du <= NEAR_M].sum() / max(lens[gap].sum(), 1e-9))
        else:
            in_unf = np.nan
        wk = ew[k]
        gw = Wd.iloc[np.unique(wk[gap])] if gap.any() else Wd.iloc[[]]
        reasons = set()
        allw = Wd.iloc[np.unique(wk)]
        if allw.restricted.any():
            reasons.add("private/no-access")
        if allw.track.any():
            reasons.add("track")
        if allw.informal.any():
            reasons.add("driveway/parking")
        if (O["barrier"][a] | O["barrier"][b]).any():
            reasons.add("gate/barrier")
        if gap.any():
            gm = mids[gap][np.argmax(lens[gap])]
            glon, glat = xy_to4326(np.array([shapely.get_x(gm)]), np.array([shapely.get_y(gm)]))
            glon, glat = float(glon[0]), float(glat[0])
        else:
            m0 = mids[np.argmin(np.where(near_link, np.inf, 0))]
            glon, glat = xy_to4326(np.array([shapely.get_x(m0)]), np.array([shapely.get_y(m0)]))
            glon, glat = float(glon[0]), float(glat[0])
        bar_nodes = np.unique(np.r_[a[O["barrier"][a]], b[O["barrier"][b]]])
        blo, bla = O["lon"][bar_nodes], O["lat"][bar_nodes]
        path_ways = allw.way_id.values
        return dict(par_share=par_share, gap_m=float(lens[gap].sum()), gap_hw=sorted(set(gw.highway)),
                    path_way_ids=";".join(str(v) for v in path_ways[:40]),
                    barriers=";".join(f"{x:.6f},{y:.6f}" for x, y in zip(blo, bla)),
                    gap_highway=";".join(sorted(set(gw.highway) | {f"service={v}" for v in gw.service if v})),
                    gap_way_ids=";".join(str(v) for v in gw.way_id.values[:8]),
                    gap_names=";".join(sorted({v for v in gw.name if v}))[:120],
                    in_overture_unfiltered=in_unf, reasons=";".join(sorted(reasons)), lon=glon, lat=glat)

    dpub = [d for d in (describe(pth) for pth in paths_pub) if d is not None]
    dful = [d for d in (describe(pth) for pth in paths_full) if d is not None]

    def second(ds):
        if len(ds) < 2:
            return None
        ds = sorted(ds, key=lambda d: d["par_share"])
        return ds[0] if ds[0]["par_share"] < PARALLEL_SHARE else "parallel"

    known = max(dful, key=lambda d: d["par_share"]) if dful else None
    sec_pub, sec_full = second(dpub), second(dful)
    cls, sub, info, note = "unresolved", "", None, ""
    if len(pocket_nodes) == 0 or len(core_nodes) == 0:
        note = f"no OSM nodes matched (pocket {len(pocket_nodes)}, core {len(core_nodes)})"
    elif f_pub >= 2 and isinstance(sec_pub, dict):
        info = sec_pub
        cls = "b"
        if info["gap_m"] <= 1.0:
            sub = "topology"
        elif set(info["gap_hw"]) <= {"service"}:
            sub = "service"
        else:
            sub = "road"
    elif f_full >= 2 and isinstance(sec_full, dict):
        info = sec_full
        cls = "c"
        sub = info["reasons"] or "restricted"
    elif f_full >= 1:
        if f_win >= 2 and f_full < 2:
            note = "second path only leaves the query window"
        else:
            cls, sub = "a", ""
            if f_full >= 2:
                note = "second path runs along the known link (same corridor)"
    else:
        note = "no OSM path between pocket and core inside the window"
    return dict(osm_class=cls, osm_subclass=sub, flow_public=f_pub, flow_full=f_full, flow_with_window_exits=f_win,
                pocket_osm_nodes=len(pocket_nodes), core_osm_nodes=len(core_nodes), osm_ways=len(Wd), osm_nodes=N,
                query_s=s, query_w=w, query_n=n, query_e=e, margin_m=margin, retrieved_utc=meta["retrieved_utc"],
                osm_base=meta["osm_base"], wire_bytes=meta["wire_bytes"], endpoint=meta["endpoint"],
                link_gap_m=None if info is None else round(info["gap_m"], 1),
                link_highway=None if info is None else info["gap_highway"],
                link_osm_way_ids=None if info is None else info["gap_way_ids"],
                link_names=None if info is None else info["gap_names"],
                link_restrictions=None if info is None else info["reasons"],
                link_in_overture_unfiltered_share=None if info is None else info["in_overture_unfiltered"],
                link_lon=None if info is None else round(info["lon"], 6),
                link_lat=None if info is None else round(info["lat"], 6),
                known_link_restrictions=None if known is None else known["reasons"],
                link_path_osm_way_ids=None if info is None else info["path_way_ids"],
                link_barriers_lonlat=None if info is None else info["barriers"], note=note)


def run_one(r, G, region, c):
    try:
        res = classify(r, G, region, c, MARGIN_M)
        if res["osm_class"] == "unresolved" and "window" in res["note"]:
            res = classify(r, G, region, c, MARGIN2_M)
    except OverpassUnavailable as ex:
        res = dict(osm_class="unresolved", osm_subclass="", note=f"Overpass unavailable: {ex}")
    except NotAttempted as ex:
        res = dict(osm_class="not_attempted", osm_subclass="", note=str(ex))
    log(f"{region} pocket {r.pocket} ({r.label or 'random'}): {res['osm_class']} {res['osm_subclass']} "
        f"flows pub/full/win {res.get('flow_public')}/{res.get('flow_full')}/{res.get('flow_with_window_exits')} "
        f"{res['note']}")
    return dict(sample=r.sample, label=r.label, region=region, GEOID=r.GEOID, pocket=int(r.pocket),
                n_bldg=int(r.n_bldg), est_pop=round(float(r.est_pop)), egress_road=r.egress_road,
                egress_class=r.egress_class, usfs_WHP_mean=round(float(r.usfs_WHP_mean), 1), whp_q=str(r.whp_q),
                overture_ext_core_share=round(float(r.ext_core_share), 3),
                overture_second_link_low_class=bool(r.ext_second_exit), lat=round(float(r.lat), 5),
                lon=round(float(r.lon), 5), **res)


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    ph = k / n
    den = 1 + z * z / n
    cen = (ph + z * z / (2 * n)) / den
    half = z * np.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return cen - half, cen + half


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    sel = select_pockets()
    if "--select-only" in sys.argv:
        print(sel[["sample", "label", "region", "GEOID", "pocket", "n_bldg", "egress_road", "whp_q"]].to_string())
        return
    if "--limit" in sys.argv:
        k = int(sys.argv[sys.argv.index("--limit") + 1])
        sel = sel.head(k)
    rows = []
    deadline = None
    for phase in ("named", "random_q4"):
        part = sel[sel["sample"] == phase]
        if phase == "random_q4":
            deadline = time.time() + RANDOM_BUDGET_S
            FETCH["attempts"], FETCH["tiles"] = 2, False
        else:
            FETCH["attempts"], FETCH["tiles"] = 4, True
        for region in REGIONS:
            sub = part[part.region == region]
            if sub.empty:
                continue
            c = con()
            tr, tg = load_tracts(region)
            region_poly, pieces = region_outline(tg)
            G = eg.build_graph(c, region, "roads", None, "d", region_poly, pieces)
            G["seg_tree"] = shapely.STRtree(G["geoms"])
            for r in sub.itertuples():
                FETCH["allow_download"] = (deadline is None or time.time() < deadline) and "--offline" not in sys.argv
                rows.append(run_one(r, G, region, c))
            del G
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "external_validation.csv", index=False)
    summ = []
    for smp, g in D.groupby("sample"):
        n_ = len(g)
        for cl in ["a", "b", "c", "unresolved", "not_attempted"]:
            k = int((g.osm_class == cl).sum())
            lo, hi = wilson(k, n_)
            summ.append(dict(sample=smp, osm_class=cl, n=k, of=n_, share=k / n_, ci_lo=lo, ci_hi=hi))
        resolved = g[g.osm_class.isin(["a", "b", "c"])]
        for cl in ["a", "b", "c"]:
            k = int((resolved.osm_class == cl).sum())
            lo, hi = wilson(k, len(resolved))
            summ.append(dict(sample=smp + " (resolved only)", osm_class=cl, n=k, of=len(resolved),
                             share=k / max(len(resolved), 1), ci_lo=lo, ci_hi=hi))
    R = D[(D["sample"] == "random_q4") & D.osm_class.isin(["a", "b", "c"])]
    if len(R):
        rng = np.random.default_rng(SEED)
        wts = R.est_pop.values.astype(float)
        idx = rng.integers(0, len(R), size=(2000, len(R)))
        for cl in ["a", "b", "c"]:
            ind = (R.osm_class.values == cl).astype(float)
            reps = (wts[idx] * ind[idx]).sum(1) / wts[idx].sum(1)
            summ.append(dict(sample="random_q4 (resolved, resident-weighted)", osm_class=cl, n=int(ind.sum()),
                             of=len(R),
                             share=(wts * ind).sum() / wts.sum(), ci_lo=np.percentile(reps, 2.5),
                             ci_hi=np.percentile(reps, 97.5)))
    S = pd.DataFrame(summ)
    S.to_csv(OUT / "external_validation_summary.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    pd.set_option("display.max_colwidth", 40)
    print(D[["sample", "label", "region", "GEOID", "n_bldg", "egress_road", "osm_class", "osm_subclass",
             "flow_public", "flow_full", "link_highway", "link_restrictions", "link_gap_m", "note"]].to_string())
    print(S.round(3).to_string())
    print("total wire bytes downloaded this run:", _spent)


if __name__ == "__main__":
    main()
