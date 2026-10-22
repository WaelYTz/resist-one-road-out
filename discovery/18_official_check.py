"""Official check of the single-egress pockets against the FHWA HPMS 2024 public road inventory.

Question per pocket: is there a second exit made of roads that a public agency officially owns? 07_external_osm.py
answered it from OSM tags (access, track, gates). Here the 'public' judgement comes from an independent official
source: the FHWA Highway Performance Monitoring System (HPMS) 2024 public release, whose road segments carry an
ownership code (1 State, 2 County, 3 Town, 4 City, 11/12 park-forest, 21/25 other agency, 26 Private, 27 Railroad,
31/32 toll, 40 other public, 50 Tribe, 60-74 Federal, 80 Other; null = unknown).

Sample: the same 56 pockets as 07 (16 named + 40 seeded random WHP-Q4 pockets, ext.select_pockets()).
Method per pocket (Overture graph rebuilt as in 02/07, so pocket ids and node matching are identical):
  1. Window = the final 07 query window (pocket edges + single link, + 750 m, or + 3000 m where 07 re-queried).
     One HPMS query per window on the state layer (envelope intersects, outSR=5070, maxAllowableOffset=2,
     geometryPrecision=0, paged by objectid). Response cached gzip-compressed with a .meta.json (URL, UTC
     retrieval time, bytes) in discovery/external/hpms_raw/.
  2. Connectivity = the cached OSM graph of 07 with the same pocket / core node matching. HPMS geometry only
     tags OSM edges; it is never used for topology (FHWA does not evaluate HPMS topology).
  3. Each OSM edge (node pair) is sampled every 5 m. Tag 'public' if >= 80% of its length lies within 15 m of an
     HPMS segment with a public ownership code (PUBLIC, and is_restricted != 1); 'private' if >= 80% within 15 m
     of ownership 26/27 (or a restricted public road); 'null' if >= 80% within 15 m of another HPMS segment
     (ownership null or 80); else 'absent' (not in HPMS). 'null' + 'absent' = 'unknown'.
     'free' edges are allowed whatever their owner: >= 80% within 15 m of the pocket's own Overture streets or of
     the known single link, plus the edges of the OSM route of the known link (the full-graph max-flow path that
     runs most along it, as in 07) that lie >= 80% within 40 m of the link. The question is about a second way
     out, not about the streets inside or the first link.
  4. Unit-capacity max-flow (edge-disjoint paths) from pocket nodes to core nodes on nested edge sets:
       G1 = free + public;  G2 = G1 + null;  G3 = G2 + absent;  full = all OSM vehicle ways (as 07).
     A second exit counts only if one path runs < 80% of its length within 40 m of the known link (as 07).
       B_off : G1 has a second exit (an officially public second exit).
       C_off : only a larger set has one; subclass 'null' (G2), 'absent' (G3) or 'private' (full graph only).
       A_off : single exit even in the full OSM graph (or the only second path runs along the known link).
       U     : no cached OSM graph, no node match, no HPMS segment in the window, or second path leaves window.
     A_off is the OSM (a) class by construction; the independent information is the B/C split.
     Diagnostic: the pair of edge-disjoint exits with the least length outside G1 (min-cost flow of 2 units on
     the full OSM graph); for its less-parallel path: non-public length and OSM classes (valid if < 80% along the
     known link). Sensitivity: C_off whose gap is <= 50 m with no OSM access restriction, track or barrier.
  5. Arizona only: ADOT ATIS Public (ADOT's 2024 HPMS submission, ownership as text) re-tags the same edges:
     public agency, open private ('PRI' unrestricted / 'OAO' open to all, code 80), gated ('OAG', 26) or unknown;
     nested sets public -> + open private -> + unknown -> + gated.
  6. C_off pockets only: USFS National Forest System roads (EDW RoadBasic) re-tag unknown edges as open NFS roads
     (openforuseto = ALL; maintenance level 3-5 passenger car, or 2 high clearance); nested flows as above.
  7. Hull check (no OSM): public HPMS features crossing the boundary of the convex hull of the pocket's streets
     buffered 100 m, more than 40 m from the known link; 'touching' if within 15 m of a pocket street.

Sources: FHWA HPMS 2024 public release, https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_{ST}_2024/
FeatureServer/0 (US federal data; FHWA notes topology not evaluated). ADOT MPD, ATIS Public,
https://services6.arcgis.com/clPWQMwZfdWn4MQZ/ArcGIS/rest/services/ADOT_ATIS_Public_view/FeatureServer/0 (as-is).
USDA Forest Service EDW, https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_RoadBasic_01/MapServer/0.
No TIGER / TIGERweb or TIGER-derived service is queried.
Usage: python -I discovery/18_official_check.py [--offline]
Outputs: discovery/out/official_validation.csv, discovery/out/official_validation_section.md
"""
from __future__ import annotations

import datetime as dt
import gzip
import importlib
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

DISC = Path(__file__).resolve().parent
sys.path.insert(0, str(DISC))  # run with python -I; downloads live in external/hpms_raw, never imported

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shapely  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import dijkstra  # noqa: E402

from common import OUT, REGIONS, con, load_tracts, log, region_outline, xy_to4326  # noqa: E402

eg = importlib.import_module("02_egress")
ext = importlib.import_module("07_external_osm")

RAW = DISC / "external" / "hpms_raw"
STATE = {"eastern-ok": "OK", "maricopa-az": "AZ", "northern-ca": "CA", "south-central-tx": "TX", "eastern-wa": "WA"}
HPMS_URL = "https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_{st}_2024/FeatureServer/0/query"
HPMS_FIELDS = ("objectid,data_year,route_id,routename,f_system,ownership,facility_type,structure_type,county_id,"
               "is_restricted")
ATIS_URL = ("https://services6.arcgis.com/clPWQMwZfdWn4MQZ/ArcGIS/rest/services/ADOT_ATIS_Public_view/"
            "FeatureServer/0/query")
ATIS_FIELDS = "OBJECTID,RouteID,RouteNameShort,Ownership,Maintenance,Jurisdiction,FunctionalSystem_Desc,County"
USFS_URL = "https://apps.fs.usda.gov/arcx/rest/services/EDW/EDW_RoadBasic_01/MapServer/0/query"
USFS_FIELDS = "objectid,id,name,oper_maint_level,openforuseto,jurisdiction,system"
USER_AGENT = "bias-bounty-mapping-equity-discovery/1.0 (research script; small bbox queries)"
PAGE = 2000
BUDGET_WIRE_BYTES = 9_000_000
PUBLIC = {1, 2, 3, 4, 11, 12, 21, 25, 31, 32, 40, 50, 60} | set(range(62, 75))
PRIVATE = {26, 27}
ROADCLS = {"motorway", "trunk", "primary", "secondary", "tertiary", "unclassified", "residential", "living_street",
           "road", "motorway_link", "trunk_link", "primary_link", "secondary_link", "tertiary_link"}
STEP_M, NEAR_M, SHARE = 5.0, 15.0, 0.8
PARALLEL_M, PARALLEL_SHARE = ext.PARALLEL_M, ext.PARALLEL_SHARE
HULL_BUF_M = 100.0
SHORT_GAP_M = 50.0
OFFLINE = "--offline" in sys.argv
MANUAL = ["Rancho Tehama", "Lake California", "Whitehorn Cove (Wagoner County)"]


# ----------------------------------------------------------------------------------------------- download
_spent = 0


class NotCached(RuntimeError):
    pass


def get_json(url: str, attempts: int = 4) -> tuple[dict, int]:
    last = None
    for attempt in range(attempts):
        if attempt:
            time.sleep(5 * 2 ** attempt)
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "gzip"})
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                wire = r.read()
                enc = r.headers.get("Content-Encoding")
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as ex:
            last = ex
            log(f"  http error (attempt {attempt}): {ex}")
            continue
        j = json.loads((gzip.decompress(wire) if enc == "gzip" else wire).decode("utf-8"))
        if "error" in j:
            last = RuntimeError(str(j["error"]))
            log(f"  service error (attempt {attempt}): {j['error']}")
            continue
        return j, len(wire)
    raise RuntimeError(f"request failed: {last}")


def fetch_layer(kind: str, key: str, url: str, fields: str, order: str, s: float, w: float, n: float, e: float):
    """Cached, paged envelope query of an ArcGIS layer. Returns (features sorted by object id, meta)."""
    global _spent
    raw_f, meta_f = RAW / f"{key}_{kind}.json.gz", RAW / f"{key}_{kind}.meta.json"
    if raw_f.exists() and meta_f.exists():
        with gzip.open(raw_f, "rt", encoding="utf-8") as fh:
            return json.load(fh), json.loads(meta_f.read_text(encoding="utf-8"))
    if OFFLINE:
        raise NotCached(f"{kind} response not cached")
    if _spent > BUDGET_WIRE_BYTES:
        raise NotCached("download byte budget exhausted")
    feats, wire, pages, urls = [], 0, 0, []
    while True:
        q = dict(geometry=f"{w:.6f},{s:.6f},{e:.6f},{n:.6f}", geometryType="esriGeometryEnvelope", inSR=4326,
                 spatialRel="esriSpatialRelIntersects", outFields=fields, returnGeometry="true", outSR=5070,
                 maxAllowableOffset=2, geometryPrecision=0, orderByFields=order, resultOffset=len(feats),
                 resultRecordCount=PAGE, f="json")
        u = url + "?" + urllib.parse.urlencode(q)
        j, b = get_json(u)
        urls.append(u)
        wire += b
        _spent += b
        pages += 1
        got = j.get("features", [])
        feats += got
        if not got or not (j.get("exceededTransferLimit") or len(got) >= PAGE):
            break
        time.sleep(1.0)
    feats.sort(key=lambda f: f["attributes"][order])
    meta = dict(kind=kind, url=url, first_query=urls[0], pages=pages, n_features=len(feats), wire_bytes=wire,
                retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"), bbox=[s, w, n, e])
    with gzip.open(raw_f, "wt", encoding="utf-8") as fh:
        json.dump(feats, fh)
    meta_f.write_text(json.dumps(meta, indent=1), encoding="utf-8")
    log(f"  fetched {key} {kind}: {len(feats)} features, {wire:,} B, {pages} page(s)")
    time.sleep(1.0)
    return feats, meta


def clean(name) -> str:
    """ADOT pads HPMS route names with '~' (e.g. '10W~SAN~JUAN~~~DR'); collapse to spaces."""
    return re.sub(r"~+", " ", str(name)).strip()


def feat_geoms(feats):
    out = []
    for f in feats:
        paths = [pth for pth in (f.get("geometry") or {}).get("paths", []) if len(pth) >= 2]
        out.append(shapely.MultiLineString(paths) if paths else None)
    return np.array(out, dtype=object)


# ----------------------------------------------------------------------------------------------- ownership
def hpms_cat(a: dict) -> str:
    o = a.get("ownership")
    if o is None:
        return "null"
    o = int(o)
    if o in PUBLIC:
        return "private" if a.get("is_restricted") == 1 else "public"
    if o in PRIVATE:
        return "private"
    return "null"  # 80 'other' and undocumented codes: ownership not established as public or private


def atis_cat(a: dict) -> str:
    t = (a.get("Ownership") or "").strip()
    if not t:
        return "unknown"
    low = t.lower()
    if "gated" in low or t.startswith("OAG"):
        return "gated"
    m = re.search(r"\((\d+)\)\s*$", t)
    code = int(m.group(1)) if m else (1 if t.startswith("DOT") else None)
    if code is None:
        return "unknown"
    if code in PUBLIC:
        return "public"
    if code == 80 or "unrestricted" in low or "open to all" in low:
        return "open_private"
    if code in PRIVATE:
        return "gated" if "restricted" in low else "open_private"
    return "unknown"


def usfs_cat(a: dict) -> str:
    if (a.get("openforuseto") or "").strip().upper() != "ALL":
        return "closed"
    ml = (a.get("oper_maint_level") or "").strip()[:1]
    return "car" if ml in ("3", "4", "5") else ("high_clearance" if ml == "2" else "closed")


# ----------------------------------------------------------------------------------------------- edge tagging
def sample_edges(O):
    a, b = O["eu"], O["ev"]
    x0, y0, x1, y1 = O["x"][a], O["y"][a], O["x"][b], O["y"][b]
    L = np.hypot(x1 - x0, y1 - y0)
    k = np.maximum(1, np.ceil(L / STEP_M)).astype(np.int64)
    eid = np.repeat(np.arange(len(a)), k)
    off = np.arange(int(k.sum())) - np.repeat(np.cumsum(k) - k, k)
    t = (off + 0.5) / np.repeat(k, k)
    px = x0[eid] + t * (x1 - x0)[eid]
    py = y0[eid] + t * (y1 - y0)[eid]
    return shapely.points(px, py), eid, k, L


def near_share(pts, eid, k, geoms, dist=NEAR_M):
    """Per edge: share of its sample points within dist of any of geoms."""
    geoms = [g for g in geoms if g is not None]
    if not geoms:
        return np.zeros(len(k))
    t = shapely.STRtree(np.array(geoms, dtype=object))
    ii, _ = t.query(pts, predicate="dwithin", distance=dist)
    hit = np.zeros(len(pts))
    hit[np.unique(ii)] = 1.0
    return np.bincount(eid, weights=hit, minlength=len(k)) / k


def tag_edges(pts, eid, k, geoms, cats, order):
    """Edge tag from per-category coverage shares; the first category in `order` reaching SHARE wins, then
    order[-1] if the union of all features covers SHARE, else 'absent'."""
    tag = np.full(len(k), "absent", dtype=object)
    sh = {c_: near_share(pts, eid, k, geoms[cats == c_]) for c_ in order}
    any_ = near_share(pts, eid, k, geoms)
    tag[any_ >= SHARE] = order[-1]
    for c_ in reversed(order):
        tag[sh[c_] >= SHARE] = c_
    return tag, any_


# ----------------------------------------------------------------------------------------------- flows
def path_info(O, kk, ag_union):
    kk = np.asarray(kk, np.int64)
    a, b = O["eu"][kk], O["ev"][kk]
    segs = shapely.linestrings(np.stack([np.c_[O["x"][a], O["y"][a]], np.c_[O["x"][b], O["y"][b]]], axis=1))
    lens = shapely.length(segs)
    mids = shapely.line_interpolate_point(segs, 0.5, normalized=True)
    par = float(lens[shapely.distance(mids, ag_union) <= PARALLEL_M].sum() / max(lens.sum(), 1e-9))
    return dict(edges=kk, lens=lens, par_share=par)


def flow_paths(O, mask, src, dst, ag_union):
    idx = np.where(mask)[0]
    f, paths = ext.max_flow_paths(len(O["ids"]), O["eu"][idx], O["ev"][idx], src, dst)
    return f, [path_info(O, idx[np.array(pth)], ag_union) for pth in paths if pth]


def second(f, ds):
    """As 07: the least-parallel path is a second exit if flow >= 2 and it runs < 80% along the known link."""
    if f < 2 or len(ds) < 2:
        return None
    d = min(ds, key=lambda q: q["par_share"])
    return d if d["par_share"] < PARALLEL_SHARE else "parallel"


def min_cost_two_exits(O, ok, src, dst, ag_union, L):
    """Two edge-disjoint pocket-to-core paths with the least total length outside `ok` (min-cost flow of 2 units,
    successive shortest paths with Johnson potentials; arc cost = length outside ok + 1e-3 x length). Returns the
    less-parallel path of the pair with its non-public length and OSM classes, or None."""
    N = len(O["ids"])
    S, T, n = N, N + 1, N + 2
    m = len(O["eu"])
    c = np.where(ok, 0.0, L) + 1e-3 * L + 1e-6
    A = pd.DataFrame(dict(fr=np.r_[O["eu"], O["ev"], np.full(len(src), S), dst],
                          to=np.r_[O["ev"], O["eu"], src, np.full(len(dst), T)],
                          c=np.r_[c, c, np.full(len(src) + len(dst), 1e-6)],
                          e=np.r_[np.arange(m), np.arange(m), np.full(len(src) + len(dst), -1)]))
    A = A[A.fr != A.to].sort_values(["fr", "to", "c", "e"]).drop_duplicates(["fr", "to"]).reset_index(drop=True)
    key = dict(zip(zip(A.fr.values.tolist(), A.to.values.tolist()), range(len(A))))

    def shortest(fr, to, w):
        M = coo_matrix((w, (fr, to)), shape=(n, n)).tocsr()
        d, pr = dijkstra(M, directed=True, indices=S, return_predecessors=True)
        if pr[T] < 0:
            return d, None
        path, cur = [], T
        while cur != S:
            path.append((int(pr[cur]), cur))
            cur = int(pr[cur])
        return d, path[::-1]

    d1, P1 = shortest(A.fr.values, A.to.values, A.c.values)
    if P1 is None:
        return None
    osm1 = [a for a in P1 if a[0] != S and a[1] != T]
    drop = np.zeros(len(A), bool)
    for a, b in osm1:  # used arc leaves the residual graph; its reverse becomes a cancellation arc
        drop[key[(a, b)]] = True
        if (b, a) in key:
            drop[key[(b, a)]] = True
    rev = pd.DataFrame([(b, a, -A.c.values[key[(a, b)]], A.e.values[key[(a, b)]]) for a, b in osm1],
                       columns=["fr", "to", "c", "e"])
    R = pd.concat([A[~drop], rev], ignore_index=True)
    reach = np.isfinite(d1)
    R = R[reach[R.fr.values] & reach[R.to.values]]
    rc = np.maximum(R.c.values + d1[R.fr.values] - d1[R.to.values], 0.0) + 1e-12
    _, P2 = shortest(R.fr.values, R.to.values, rc)
    if P2 is None:
        return None
    flow = {}
    for a, b in P1 + P2:
        if flow.get((b, a), 0) > 0:
            flow[(b, a)] -= 1
        else:
            flow[(a, b)] = flow.get((a, b), 0) + 1
    out = {}
    for (a, b), f in sorted(flow.items()):
        if f > 0:
            out.setdefault(a, []).extend([b] * f)
    paths = []
    for _ in range(2):
        cur, pth = S, []
        while cur != T and out.get(cur):
            b = out[cur].pop()
            if cur != S and b != T:
                pth.append(int(A.e.values[key[(cur, b)]]))
            cur = b
        paths.append(pth)
    ds = [path_info(O, pth, ag_union) for pth in paths if pth]
    if len(ds) < 2:
        return None
    d = min(ds, key=lambda q: q["par_share"])
    kk = d["edges"]
    bad = ~ok[kk]
    W = O["W"]
    hw = W.highway.values[O["ew"][kk]].astype(str)
    tr = ((W.restricted.values | W.informal.values | W.track.values)[O["ew"][kk]] | O["barrier"][O["eu"][kk]]
          | O["barrier"][O["ev"][kk]])
    lab = np.where(tr & ~W.track.values[O["ew"][kk]], np.char.add(hw, "(restricted)"), hw)
    by = pd.Series(d["lens"][bad]).groupby(lab[bad]).sum().sort_values(ascending=False)
    d.update(nonpublic_m=float(d["lens"][bad].sum()), nonpublic_restricted=bool((bad & tr).any()),
             nonpublic_osm=";".join(f"{k_} {v:.0f} m" for k_, v in by.items()))
    return d


def describe_path(d, tag, free, H, Hg, O):
    """Length (m) by tag on the non-free part of a path, and the HPMS routes it follows."""
    kk, lens = d["edges"], d["lens"]
    nf = ~free[kk]
    by = {t: round(float(lens[nf & (tag[kk] == t)].sum())) for t in ("public", "null", "absent", "private")}
    routes = []
    ok = np.array([g is not None for g in Hg], bool)
    if nf.any() and ok.any():
        a, b = O["eu"][kk[nf]], O["ev"][kk[nf]]
        mids = shapely.points((O["x"][a] + O["x"][b]) / 2, (O["y"][a] + O["y"][b]) / 2)
        (qi, ti), dd = shapely.STRtree(Hg[ok]).query_nearest(mids, return_distance=True, all_matches=False)
        for h in np.where(ok)[0][ti[dd <= NEAR_M]]:
            at = H[h]["attributes"]
            r = f"{clean(at.get('routename') or at.get('route_id'))} (own {at.get('ownership')})"
            if r not in routes:
                routes.append(r)
    return by, ";".join(routes[:6])


def hull_crossings(pg, ag_union, Hg, cats, H):
    """Public HPMS routes crossing the boundary of the pocket hull + 100 m away from the known link."""
    if len(pg) == 0:
        return []
    pocket_union = shapely.union_all(pg)
    ring = shapely.boundary(shapely.convex_hull(pocket_union).buffer(HULL_BUF_M))
    rows = []
    for i in np.where(cats == "public")[0]:
        g = Hg[i]
        if g is None or not shapely.intersects(g, ring):
            continue
        cr = [q for q in shapely.get_parts(shapely.intersection(g, ring))
              if q.geom_type == "Point" and shapely.distance(q, ag_union) > PARALLEL_M]
        if not cr:
            continue
        at = H[i]["attributes"]
        lo, la = xy_to4326(np.array([cr[0].x]), np.array([cr[0].y]))
        rows.append(dict(route=clean(at.get("routename") or at.get("route_id")), ownership=at.get("ownership"),
                         f_system=at.get("f_system"), crossings=len(cr), lon=round(float(lo[0]), 5),
                         lat=round(float(la[0]), 5), dist_m=float(shapely.distance(g, pocket_union)),
                         dist_link_m=float(shapely.distance(g, ag_union))))
    if not rows:
        return []
    R = pd.DataFrame(rows).sort_values(["dist_m", "route"])
    out = []
    for (route, own), g in R.groupby(["route", "ownership"], sort=False, dropna=False):
        f0 = g.iloc[0]
        out.append(dict(route=route, ownership=own, f_system=f0.f_system, segments=len(g), crossings=int(g.crossings.sum()),
                        lat=f0.lat, lon=f0.lon, dist_to_pocket_street_m=round(float(g.dist_m.min()), 1),
                        dist_to_link_m=round(float(g.dist_link_m.min()), 1), touching=bool(g.dist_m.min() <= NEAR_M)))
    return out


# ----------------------------------------------------------------------------------------------- per pocket
def check_pocket(row, G, region, c):
    pid = int(row.pocket)
    pk = G["pocket"]
    real = np.arange(len(G["eu"]))
    eu_, ev_ = G["eu"], G["ev"]
    in_p = (pk[eu_] == pid) & (pk[ev_] == pid)
    att = ((pk[eu_] == pid) & G["core"][ev_]) | ((pk[ev_] == pid) & G["core"][eu_])
    assert int((pk[:-1] == pid).sum()) == int(row.n_nodes), ("pocket id mismatch", region, pid)
    pe, ae = real[in_p], real[att]
    pg = ext.edge_geoms(G, pe) if len(pe) else np.array([], dtype=object)
    ag = ext.edge_geoms(G, ae)
    margin = float(row.margin_m) if pd.notna(row.margin_m) else ext.MARGIN_M
    xmin, ymin, xmax, ymax = shapely.total_bounds(np.r_[pg, ag])
    X0, Y0, X1, Y1 = xmin - margin, ymin - margin, xmax + margin, ymax + margin
    lo, la = xy_to4326(np.array([X0, X1, X0, X1]), np.array([Y0, Y0, Y1, Y1]))
    s, w, n, e = la.min(), lo.min(), la.max(), lo.max()
    if pd.notna(row.query_s):
        assert np.allclose([s, w, n, e], [row.query_s, row.query_w, row.query_n, row.query_e], atol=1e-9), \
            ("window differs from 07", region, pid)
    key = f"{region}_p{pid}"
    res = dict(query_s=s, query_w=w, query_n=n, query_e=e, margin_m=margin)
    ag_union = shapely.union_all(ag)
    # --- HPMS
    try:
        H, hm = fetch_layer("hpms", key, HPMS_URL.format(st=STATE[region]), HPMS_FIELDS, "objectid", s, w, n, e)
    except NotCached as ex:
        return dict(res, off_class="U", off_subclass="", note=f"HPMS: {ex}")
    Hg = feat_geoms(H)
    hc = np.array([hpms_cat(f["attributes"]) for f in H], dtype=object)
    own = pd.Series([f["attributes"].get("ownership") for f in H], dtype="float")
    res.update(hpms_features=len(H), hpms_public_features=int((hc == "public").sum()),
               hpms_private_features=int((hc == "private").sum()), hpms_null_ownership_features=int(own.isna().sum()),
               hpms_restricted_features=int(sum(f["attributes"].get("is_restricted") == 1 for f in H)),
               hpms_data_year=";".join(sorted({str(f["attributes"].get("data_year")) for f in H})),
               hpms_retrieved_utc=hm["retrieved_utc"], hpms_wire_bytes=hm["wire_bytes"], hpms_pages=hm["pages"])
    hull = hull_crossings(pg, ag_union, Hg, hc, H)
    res.update(hull_public_routes=len(hull), hull_public_touching_routes=sum(h["touching"] for h in hull),
               hull_public_route_names=";".join(h["route"] for h in hull)[:200],
               hull_touching_route_names=";".join(h["route"] for h in hull if h["touching"])[:200])
    res["_hull"] = hull
    res["pocket_street_km"] = round(float(shapely.length(shapely.union_all(pg))) / 1000, 1) if len(pg) else 0.0
    res["link_m"] = round(float(shapely.length(ag_union)))
    # --- ADOT ATIS (Arizona)
    A = Ag = ac = None
    if region == "maricopa-az":
        try:
            A, am = fetch_layer("atis", key, ATIS_URL, ATIS_FIELDS, "OBJECTID", s, w, n, e)
            Ag = feat_geoms(A)
            ac = np.array([atis_cat(f["attributes"]) for f in A], dtype=object)
            res.update(atis_features=len(A), atis_gated_features=int((ac == "gated").sum()),
                       atis_retrieved_utc=am["retrieved_utc"], atis_wire_bytes=am["wire_bytes"])
        except NotCached as ex:
            res["atis_note"] = str(ex)
    # --- OSM graph from the 07 cache (no download)
    ext.FETCH["allow_download"] = False
    try:
        j, _ = ext.fetch(key, s, w, n, e, margin)
    except (ext.NotAttempted, ext.OverpassUnavailable) as ex:
        return dict(res, off_class="U", off_subclass="", note=f"no cached OSM graph ({ex})")
    O = ext.osm_graph(j)
    box = shapely.box(X0, Y0, X1, Y1)
    seg_mask = np.zeros(G["n"], bool)
    seg_mask[G["seg_tree"].query(box, predicate="intersects")] = True
    ce = real[seg_mask[G["eseg"]] & G["core"][eu_] & G["core"][ev_]]
    cg = ext.edge_geoms(G, ce) if len(ce) else np.array([], dtype=object)
    pts_n = shapely.points(O["x"], O["y"])

    def dist_to(geoms):
        if len(geoms) == 0:
            return np.full(len(pts_n), np.inf)
        (qi, _), d = shapely.STRtree(geoms).query_nearest(pts_n, return_distance=True, all_matches=False)
        out = np.full(len(pts_n), np.inf)
        out[qi] = d
        return out

    d_p, d_c, d_a = dist_to(pg), dist_to(cg), dist_to(ag)
    src = np.where((d_p <= ext.NEAR_M) & (d_c > ext.FAR_M) & (d_a > ext.LINK_M))[0]
    dst = np.where((d_c <= ext.NEAR_M) & (d_p > ext.FAR_M) & (d_a > ext.LINK_M))[0]
    inside = (O["lon"] >= w) & (O["lon"] <= e) & (O["lat"] >= s) & (O["lat"] <= n)
    out_nodes = np.where(~inside)[0]
    res.update(pocket_osm_nodes=len(src), core_osm_nodes=len(dst), osm_edges=len(O["eu"]))
    if len(src) == 0 or len(dst) == 0:
        return dict(res, off_class="U", off_subclass="", note="no OSM node match")
    if len(H) == 0:
        return dict(res, off_class="U", off_subclass="", note="no HPMS segment in window")
    pts, eid, k, L = sample_edges(O)
    tag, any_ = tag_edges(pts, eid, k, Hg, hc, ["public", "private", "null"])
    m_all = np.ones(len(tag), bool)
    ff, dsf = flow_paths(O, m_all, src, dst, ag_union)
    sf = second(ff, dsf)
    known = max(dsf, key=lambda q: q["par_share"]) if dsf else None
    link_route = np.zeros(len(tag), bool)
    if known is not None:
        link_route[known["edges"]] = True
    near_own = near_share(pts, eid, k, list(np.r_[pg, ag])) >= SHARE
    corr40 = near_share(pts, eid, k, list(ag), dist=PARALLEL_M) >= SHARE  # the known link's corridor (07: 40 m)
    link_free = link_route & corr40
    free = near_own | link_free
    on_link = (near_share(pts, eid, k, list(ag)) >= SHARE) | link_free
    pocket_e = near_own & ~on_link
    roadcls = np.isin(O["W"].highway.values[O["ew"]], list(ROADCLS)) & ~O["W"].restricted.values[O["ew"]]
    Lk, Lp, Lr = max(L[on_link].sum(), 1e-9), max(L[pocket_e].sum(), 1e-9), max(L[roadcls].sum(), 1e-9)
    res.update(window_osm_km=round(float(L.sum()) / 1000, 1),
               window_osm_road_km=round(float(L[roadcls].sum()) / 1000, 1),
               window_road_hpms_share=round(float(L[roadcls & (any_ >= SHARE)].sum() / Lr), 3),
               window_road_public_share=round(float(L[roadcls & (tag == "public")].sum() / Lr), 3),
               link_public_share=round(float(L[on_link & (tag == "public")].sum() / Lk), 3),
               link_hpms_share=round(float(L[on_link & (any_ >= SHARE)].sum() / Lk), 3),
               pocket_street_hpms_share=round(float(L[pocket_e & (any_ >= SHARE)].sum() / Lp), 3),
               pocket_street_public_share=round(float(L[pocket_e & (tag == "public")].sum() / Lp), 3),
               pocket_street_private_share=round(float(L[pocket_e & (tag == "private")].sum() / Lp), 3))
    lk = []
    if on_link.any() and any(g is not None for g in Hg):
        mids = shapely.points((O["x"][O["eu"][on_link]] + O["x"][O["ev"][on_link]]) / 2,
                              (O["y"][O["eu"][on_link]] + O["y"][O["ev"][on_link]]) / 2)
        okh = np.array([g is not None for g in Hg], bool)
        (_, ti), dd = shapely.STRtree(Hg[okh]).query_nearest(mids, return_distance=True, all_matches=False)
        for h in np.unique(np.where(okh)[0][ti[dd <= NEAR_M]]):
            a_ = H[h]["attributes"]
            lk.append(f"{clean(a_.get('routename') or a_.get('route_id'))} (own {a_.get('ownership')})")
    res["link_hpms"] = ";".join(sorted(set(lk))) if lk else "not in HPMS"
    g1 = free | (tag == "public")
    g2 = g1 | (tag == "null")
    g3 = g2 | (tag == "absent")
    f1, ds1 = flow_paths(O, g1, src, dst, ag_union)
    f2, ds2 = flow_paths(O, g2, src, dst, ag_union)
    f3, ds3 = flow_paths(O, g3, src, dst, ag_union)
    s1, s2, s3 = second(f1, ds1), second(f2, ds2), second(f3, ds3)
    fw, _ = ext.max_flow_paths(len(O["ids"]), O["eu"], O["ev"], src, np.r_[dst, out_nodes])
    res.update(flow_g1_public=f1, flow_g2_null=f2, flow_g3_nonprivate=f3, flow_full=ff, flow_window_exits=fw)
    cls, sub, info, note = "U", "", None, ""
    if isinstance(s1, dict):
        cls, info = "B_off", s1
    elif isinstance(s2, dict):
        cls, sub, info = "C_off", "null", s2
    elif isinstance(s3, dict):
        cls, sub, info = "C_off", "absent", s3
    elif isinstance(sf, dict):
        cls, sub, info = "C_off", "private", sf
    elif ff >= 1:
        if fw >= 2 and ff < 2:
            note = "second path only leaves the query window"
        else:
            cls = "A_off"
            note = "second path runs along the known link" if ff >= 2 else ""
    else:
        note = "no OSM path between pocket and core inside the window"
    if info is not None:
        by, routes = describe_path(info, tag, free, H, Hg, O)
        res.update(second_path_public_m=by["public"], second_path_null_m=by["null"],
                   second_path_absent_m=by["absent"], second_path_private_m=by["private"],
                   second_path_hpms_routes=routes)
    if cls in ("B_off", "C_off"):
        mn = min_cost_two_exits(O, g1, src, dst, ag_union, L)
        if mn is not None:
            res.update(min_nonpublic_m=round(mn["nonpublic_m"], 1), min_nonpublic_par_share=round(mn["par_share"], 3),
                       min_nonpublic_restricted=mn["nonpublic_restricted"], min_nonpublic_osm=mn["nonpublic_osm"])
    # --- ATIS re-tagging (Arizona)
    if A is not None:
        if len(A):
            atag, _ = tag_edges(pts, eid, k, Ag, ac, ["public", "gated", "open_private", "unknown"])
            atag = np.where(atag == "absent", "unknown", atag)
            a1 = free | (atag == "public")
            a2 = a1 | (atag == "open_private")
            a3 = a2 | (atag == "unknown")
            sa = [second(*flow_paths(O, m_, src, dst, ag_union)) for m_ in (a1, a2, a3)]
            if isinstance(sa[0], dict):
                acl = "B_atis"
            elif isinstance(sa[1], dict):
                acl = "B_atis_open_private"
            elif isinstance(sa[2], dict):
                acl = "C_atis_unknown"
            elif isinstance(sf, dict):
                acl = "C_atis_gated_only"
            else:
                acl = cls if cls in ("A_off", "U") else "A_off"
            res.update(atis_class=acl,
                       atis_window_road_gated_share=round(float(L[roadcls & (atag == "gated")].sum() / Lr), 3))
        else:
            res.update(atis_class="U")
    # --- USFS National Forest System roads (C_off pockets only)
    if cls == "C_off":
        try:
            F, _ = fetch_layer("usfs", key, USFS_URL, USFS_FIELDS, "objectid", s, w, n, e)
        except NotCached as ex:
            F = None
            res["usfs_note"] = str(ex)
        if F is not None:
            Fg = feat_geoms(F)
            fc = np.array([usfs_cat(f["attributes"]) for f in F], dtype=object)
            res["usfs_features"] = len(F)
            if len(F):
                ftag, _ = tag_edges(pts, eid, k, Fg, fc, ["car", "high_clearance", "closed"])
                unk = (tag == "null") | (tag == "absent")
                u1 = g1 | (unk & (ftag == "car"))
                u2 = u1 | (unk & (ftag == "high_clearance"))
                su = [second(*flow_paths(O, m_, src, dst, ag_union)) for m_ in (u1, u2)]
                res["usfs_class"] = ("B_usfs_car" if isinstance(su[0], dict) else
                                     "B_usfs_high_clearance" if isinstance(su[1], dict) else "C_off")
            else:
                res["usfs_class"] = "no NFS road in window"
    return dict(res, off_class=cls, off_subclass=sub, note=note)


# ----------------------------------------------------------------------------------------------- summary
def wilson(k, n, z=1.96):
    return ext.wilson(k, n, z)


def kish_wilson(ind, wts):
    wts = np.asarray(wts, float)
    ind = np.asarray(ind, float)
    p_ = (wts * ind).sum() / wts.sum()
    neff = wts.sum() ** 2 / (wts ** 2).sum()
    lo, hi = wilson(p_ * neff, neff)
    return p_, lo, hi, neff


def kappa(a, b, labels):
    M = pd.crosstab(pd.Categorical(a, labels), pd.Categorical(b, labels), dropna=False).values.astype(float)
    n_ = M.sum()
    po = np.trace(M) / n_
    pe = (M.sum(0) * M.sum(1)).sum() / n_ ** 2
    return (po - pe) / (1 - pe) if pe < 1 else np.nan, int(n_)


def pct(x):
    return f"{100 * x:.1f}%"


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    V = pd.read_csv(OUT / "external_validation.csv")
    sel = ext.select_pockets()
    assert sorted(V.region + ":" + V.pocket.astype(str)) == sorted(sel.region + ":" + sel.pocket.astype(str)), \
        "sample differs from 07"
    sel = sel.merge(V[["region", "pocket", "osm_class", "osm_subclass", "margin_m", "query_s", "query_w", "query_n",
                       "query_e"]], on=["region", "pocket"], how="left")
    rows, hulls = [], {}
    for region in REGIONS:
        sub = sel[sel.region == region]
        if sub.empty:
            continue
        c = con()
        tr, tg = load_tracts(region)
        region_poly, pieces = region_outline(tg)
        G = eg.build_graph(c, region, "roads", None, "d", region_poly, pieces)
        G["seg_tree"] = shapely.STRtree(G["geoms"])
        for r in sub.itertuples():
            res = check_pocket(r, G, region, c)
            hull = res.pop("_hull", [])
            if r.label:
                hulls[r.label] = hull
            log(f"{region} p{r.pocket} ({r.label or 'random'}): OSM {r.osm_class} -> {res['off_class']} "
                f"{res['off_subclass']} flows g1/g2/g3/full {res.get('flow_g1_public')}/{res.get('flow_g2_null')}/"
                f"{res.get('flow_g3_nonprivate')}/{res.get('flow_full')} min-nonpublic {res.get('min_nonpublic_m')} "
                f"[{res.get('min_nonpublic_osm', '')}] atis {res.get('atis_class', '')} usfs {res.get('usfs_class', '')}"
                f" {res['note']}")
            rows.append(dict(sample=r.sample, label=r.label, region=region, GEOID=r.GEOID, pocket=int(r.pocket),
                             n_bldg=int(r.n_bldg), est_pop=round(float(r.est_pop)), egress_road=r.egress_road,
                             whp_q=str(r.whp_q), lat=round(float(r.lat), 5), lon=round(float(r.lon), 5),
                             osm_class=r.osm_class, osm_subclass=r.osm_subclass, **res))
        del G
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "official_validation.csv", index=False)
    write_section(D, hulls)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_colwidth", 50)
    print(D[["sample", "label", "region", "pocket", "est_pop", "osm_class", "osm_subclass", "off_class",
             "off_subclass", "flow_g1_public", "flow_full", "min_nonpublic_m", "window_road_hpms_share",
             "hull_public_touching_routes", "note"]].to_string())
    print("total wire bytes downloaded this run:", _spent)


def key_rows(D, label, cls_col):
    out = []
    for smp, name in (("random_q4", "Random WHP-Q4 sample"), ("named", "Named communities")):
        g = D[(D["sample"] == smp) & (D.off_class != "U")]
        ind = (~g[cls_col].astype(str).str.startswith("B")).values
        lo, hi = wilson(ind.sum(), len(g))
        pw, wlo, whi, neff = kish_wilson(ind, g.est_pop.values)
        out.append(f"| {name}{label} | {int(ind.sum())} of {len(g)} | {pct(ind.mean())} ({pct(lo)}–{pct(hi)}) | "
                   f"{pct(pw)} ({pct(wlo)}–{pct(whi)}; n_eff {neff:.0f}) |")
    return out


def write_section(D, hulls):
    metas = [json.loads(f.read_text(encoding="utf-8")) for f in sorted(RAW.glob("*.meta.json"))]
    by_kind = {k_: [m for m in metas if m["kind"] == k_] for k_ in ("hpms", "atis", "usfs")}
    hm = by_kind["hpms"]
    t0 = min(m["retrieved_utc"] for m in metas)[:16].replace("T", " ")
    t1 = max(m["retrieved_utc"] for m in metas)[:16].replace("T", " ")
    mb = {k_: sum(m["wire_bytes"] for m in v) / 1e6 for k_, v in by_kind.items()}
    hn = sum(m["n_features"] for m in hm)
    years = sorted({y for v in D.hpms_data_year.dropna() for y in str(v).split(";") if y != "None"})
    res = D[D.off_class != "U"]
    both = res[res.osm_class.isin(["a", "b", "c"])]
    k3, n3 = kappa(both.osm_class.map({"a": "A", "b": "B", "c": "C"}), both.off_class.str[0], ["A", "B", "C"])
    k2, _ = kappa(np.where(both.osm_class == "b", "pub", "no"), np.where(both.off_class == "B_off", "pub", "no"),
                  ["pub", "no"])
    bc = both[both.osm_class.isin(["b", "c"])]
    kbc, nbc = kappa(bc.osm_class.map({"b": "B", "c": "C"}), bc.off_class.str[0], ["B", "C"])
    n_bC = int(((both.osm_class == "b") & (both.off_class == "C_off")).sum())
    n_cB = int(((both.osm_class == "c") & (both.off_class == "B_off")).sum())
    agree = int((both.osm_class.map({"a": "A", "b": "B", "c": "C"}) == both.off_class.str[0]).sum())
    comp = res.groupby("region").window_road_hpms_share.median()
    L = []
    L.append("#### Official check (FHWA HPMS 2024)\n")
    L.append(f"- **Source.** FHWA Highway Performance Monitoring System, 2024 public release (field `data_year` = "
             f"{', '.join(years)}): `https://geo.dot.gov/server/rest/services/Hosted/HPMS_FULL_{{ST}}_2024/FeatureServer/0` "
             f"for ST = TX, CA, WA, AZ, OK. {len(hm)} envelope queries ({hn:,} segments), {mb['hpms']:.2f} MB gzip. "
             f"US federal data; FHWA notes that topology is not evaluated, so HPMS is used only for ownership, never "
             f"for connectivity. Arizona cross-check: ADOT MPD, ATIS Public (ADOT's 2024 HPMS submission with "
             f"ownership as text; as-is disclaimer), {len(by_kind['atis'])} queries, {mb['atis']:.2f} MB. Forest "
             f"roads: USDA Forest Service EDW RoadBasic (`https://apps.fs.usda.gov/arcx/rest/services/EDW/"
             f"EDW_RoadBasic_01/MapServer/0`), {len(by_kind['usfs'])} queries, {mb['usfs']:.2f} MB. All retrieved "
             f"{t0}Z to {t1}Z (UTC) and cached with metadata in `discovery/external/hpms_raw/`; script "
             f"`18_official_check.py` (`--offline` re-runs from the cache).")
    L.append("- **Method.** Same 56 pockets, windows, cached OSM graph and pocket/core node matching as above. Every "
             "OSM edge is sampled every 5 m and tagged *public* if at least 80% of its length lies within 15 m of an "
             "HPMS segment with a public ownership code (1–4, 11, 12, 21, 25, 31, 32, 40, 50, 60, 62–74), *private* "
             "for 26/27, and *unknown* otherwise (ownership null or 80, or not in HPMS). The pocket's own streets and "
             "the known single link count whatever their owner. **B_off**: a max-flow second exit exists using only "
             "public edges outside the pocket, not along the known link (< 80% of its length within 40 m of it). "
             "**C_off**: two exits only once unknown or private edges are allowed. **A_off**: a single exit even in "
             "the full OSM graph. **U**: no cached OSM graph, or no node match. A_off equals OSM (a) by construction; "
             "the independent information is whether the second exit is officially public.")
    tab = pd.crosstab(D.osm_class.fillna("unresolved"), D.off_class)
    cols = [c_ for c_ in ["A_off", "B_off", "C_off", "U"] if c_ in tab.columns]
    tab = tab.reindex(columns=cols, fill_value=0)
    lab = {"a": "(a) single access in OSM", "b": "(b) public second link in OSM", "c": "(c) only private/gated/track",
           "unresolved": "OSM unresolved"}
    L.append("")
    L.append("| OSM class (all 56) | " + " | ".join(cols) + " |")
    L.append("|---|" + "---|" * len(cols))
    for idx in ["a", "b", "c", "unresolved"]:
        if idx in tab.index:
            L.append(f"| {lab[idx]} | " + " | ".join(str(int(v)) for v in tab.loc[idx].values) + " |")
    L.append("")
    L.append(f"- **Agreement.** {agree} of {n3} pockets resolved in both agree on A/B/C. Cohen's kappa is {k3:.2f} "
             f"for 3 classes and {k2:.2f} for 'public second exit, yes/no'. On the {nbc} pockets where OSM has a second "
             f"exit (b vs c, the part not shared by construction), kappa is {kbc:.2f}. Disagreements: {n_bC} pockets "
             f"where OSM calls the second link public but HPMS does not list it as public (b to C_off), and {n_cB} "
             f"the other way (c to B_off).")
    L.append("")
    L.append("| No official public second exit (A_off + C_off) | Pockets | Share (95% Wilson CI) | Resident-weighted "
             "(Wilson CI, Kish n_eff) |")
    L.append("|---|---|---|---|")
    L += key_rows(D, "", "off_class")
    sens_notes = []
    D2 = D.copy()
    short = ((D2.off_class == "C_off") & (D2.min_nonpublic_m <= SHORT_GAP_M)
             & (D2.min_nonpublic_par_share < PARALLEL_SHARE) & (D2.min_nonpublic_restricted == False))  # noqa: E712
    D2["sens_gap"] = np.where(short, "B_gap", D2.off_class)
    if short.any():
        L += key_rows(D2, f", open gaps ≤ {SHORT_GAP_M:.0f} m counted as public", "sens_gap")
    else:
        sens_notes.append(f"counting open gaps of ≤ {SHORT_GAP_M:.0f} m as public")
    if "usfs_class" in D2:
        usb = D2.usfs_class.astype(str).str.startswith("B")
        D2["sens_usfs"] = np.where(usb, "B_usfs", D2.off_class)
        if usb.any():
            L += key_rows(D2, ", + open NFS roads (USFS)", "sens_usfs")
        else:
            sens_notes.append("counting open National Forest System roads (USFS, maintenance level 2–5) as public")
    L.append("")
    R = D[(D["sample"] == "random_q4") & D.osm_class.isin(["a", "b", "c"]) & (D.off_class != "U")]
    k_osm = int(R.osm_class.isin(["a", "c"]).sum())
    lo_, hi_ = wilson(k_osm, len(R))
    L.append(f"- **Reading.** For comparison, OSM tags alone give {k_osm} of {len(R)} ({pct(k_osm / len(R))}, "
             f"{pct(lo_)}–{pct(hi_)}) random pockets without a public second road. The official inventory is "
             f"stricter, because public service roads in OSM are not in HPMS. "
             + (("The result does not change when " + " or ".join(sens_notes) + ".") if sens_notes else ""))
    C = D[D.off_class == "C_off"]
    cs = C.off_subclass.value_counts()
    gaps = C.min_nonpublic_m.dropna()
    nres = int((C.min_nonpublic_restricted == True).sum())  # noqa: E712
    L.append(f"- **C_off detail ({len(C)} pockets).** The second exit appears once these are allowed: roads absent "
             f"from HPMS in {int(cs.get('absent', 0))} pockets, roads with null or other HPMS ownership in "
             f"{int(cs.get('null', 0))}, HPMS private roads in {int(cs.get('private', 0))}. Taking the cheapest pair "
             f"of exits (min-cost flow), the second exit's non-public stretch has a median of {gaps.median():.0f} m "
             f"(range {gaps.min():.0f}–{gaps.max():.0f} m; {int((gaps > 500).sum())} pockets over 500 m). In {nres} "
             f"of them, OSM also tags that stretch as a track, private, a driveway or gated. "
             f"{int(short.sum())} pockets have an open gap of ≤ {SHORT_GAP_M:.0f} m, which could be a matching artefact.")
    bC = C[C.osm_class == "b"].dropna(subset=["min_nonpublic_osm"])
    if len(bC):
        L.append(f"  - The {len(bC)} OSM-(b) pockets that become C_off: their second link is public in OSM but not in "
                 "HPMS. " + "; ".join(f"{r.label or r.egress_road} ({r.min_nonpublic_osm.replace(';', ' + ')})"
                                      for r in bC.itertuples()) + ".")
    for smp in ("named", "random_q4"):
        g = D[(D["sample"] == smp) & (D.off_class == "B_off")]
        if len(g):
            L.append(f"- **B_off ({smp}).** " + "; ".join(
                f"{(r.label or r.egress_road)} ({r.region}): {r.second_path_hpms_routes or 'via pocket streets'}"
                for r in g.itertuples()))
    L.append("- **HPMS completeness varies by state.** This is the median share of OSM road-class length "
             "(residential and above, not access-restricted) in the windows that lies within 15 m of any HPMS "
             "segment: " + ", ".join(f"{STATE[r_]} {v:.0%}" for r_, v in comp.sort_values(ascending=False).items())
             + ". Where coverage is low, 'absent from HPMS' includes unlisted public local roads, not only private "
               "ones.")
    az = D[(D.region == "maricopa-az") & D.get("atis_class", pd.Series(index=D.index, dtype=object)).notna()]
    if len(az):
        vc = az.atis_class.value_counts()
        gated = az[az.atis_class == "C_atis_gated_only"]
        nwin = int((D[D.region == "maricopa-az"].atis_gated_features.fillna(0) > 0).sum())
        L.append(f"- **Arizona, ADOT ATIS ownership text.** {len(az)} resolved AZ pockets: " +
                 ", ".join(f"{k_} {int(v)}" for k_, v in vc.sort_index().items()) +
                 f". Gated HOA segments ('OAG-Owners Association - Gated') appear in {nwin} of "
                 f"{int((D.region == 'maricopa-az').sum())} AZ windows. Pockets whose only second exit runs "
                 f"through them: {len(gated)}, " +
                 (", ".join(f"{r.egress_road} (pocket {r.pocket}, OSM {r.osm_class})" for r in gated.itertuples())
                  or "none") + ".")
    if "usfs_class" in D:
        u = D[D.usfs_class.notna()]
        vc = u.usfs_class.value_counts()
        wn = u[u.usfs_features.fillna(0) > 0]
        ub = u[u.usfs_class.astype(str).str.startswith("B")]
        L.append(f"- **Forest roads (USFS).** National Forest System roads lie in {len(wn)} of the {len(u)} C_off "
                 f"windows ({', '.join(str(r.label or r.egress_road) for r in wn.itertuples())}). " +
                 ("When open NFS roads count as public, these pockets gain a second exit: " +
                  "; ".join(f"{r.label or r.egress_road} ({r.usfs_class})" for r in ub.itertuples()) + "."
                  if len(ub) else "Counting open NFS roads (maintenance level 2–5) as public gives none of them a "
                                  "second exit."))
    L.append("- **Manual check.** These are public HPMS routes crossing the boundary of the pocket hull (+ 100 m), "
             "more than 40 m from the known link. They come from HPMS geometry only, with no OSM. F = HPMS "
             "functional system. A route a few metres from the link is the road the link joins.")
    for lb in MANUAL:
        hh = hulls.get(lb, [])
        r = D[D.label == lb].iloc[0]
        head = (f"  - {lb}: {r.pocket_street_km:.0f} km of pocket streets, {r.pocket_street_hpms_share:.0%} of them "
                f"in HPMS; HPMS along the {r.link_m:.0f} m link: {r.link_hpms}")
        if not hh:
            L.append(f"{head}. Crossing routes: none.")
        else:
            L.append(f"{head}. Crossing routes: " + "; ".join(
                f"{h['route']} (own {h['ownership']}, F{h['f_system']}, {h['dist_to_link_m']:.0f} m from the link, "
                f"{h['dist_to_pocket_street_m']:.0f} m from the nearest pocket street"
                + (", touching" if h["touching"] else "") + ")" for h in hh) +
                (". None comes within 15 m of a pocket street." if not any(h["touching"] for h in hh) else "."))
    nr = int(D.hpms_restricted_features.fillna(0).sum())
    L.append("- **Limitations.** HPMS lists public roads only. Private roads, driveways and most tracks are absent, "
             "and so are many public local roads where state coverage is low (above). So 'unknown' mixes private "
             "roads with unlisted public ones, and C_off is an upper bound on 'no public second exit'. Ownership is "
             f"null on {D.hpms_null_ownership_features.sum() / max(D.hpms_features.sum(), 1):.0%} of the segments "
             f"in the windows. `is_restricted` is set on {nr} of them. The vintage is 2024 (OSM 2026, Overture "
             "2025). HPMS ownership says nothing about gates or seasonal closures, so a public road may still be "
             "gated. Matching uses 15 m buffers on HPMS geometry simplified to 2 m, so offsets beyond 15 m read as "
             "unknown.")
    (OUT / "official_validation_section.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    log("wrote", OUT / "official_validation_section.md")


if __name__ == "__main__":
    main()
