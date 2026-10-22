"""Shared helpers for the road-topology discovery (Overture-only, challenge files only)."""
from __future__ import annotations

import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyproj
import shapely
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
DISC = ROOT / "discovery"
OUT = DISC / "out"
CACHE = DISC / "cache"
OUT.mkdir(exist_ok=True, parents=True)
CACHE.mkdir(exist_ok=True, parents=True)

REGIONS = ["eastern-ok", "maricopa-az", "northern-ca", "south-central-tx", "eastern-wa"]

DRIVABLE = ["motorway", "trunk", "primary", "secondary", "tertiary", "residential", "unclassified"]
# classes that only exist in the -unfiltered layer and can carry a car (service roads, tracks, unknown roads)
EXTRA_VEHICLE = ["service", "track", "unknown", "living_street"]
# non-vehicle classes
NONVEH = ["footway", "path", "cycleway", "steps", "pedestrian", "bridleway"]

EDGE_BUF_M = 100.0  # components with a segment within this distance of the region outline are "edge" components

_TR = pyproj.Transformer.from_crs("EPSG:4326", "EPSG:5070", always_xy=True)
_TR_INV = pyproj.Transformer.from_crs("EPSG:5070", "EPSG:4326", always_xy=True)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def con():
    c = duckdb.connect()
    c.execute("LOAD spatial;")
    c.execute("SET threads TO 8;")
    return c


def p(region: str, kind: str) -> str:
    if kind in ("tracts",):
        return str(ROOT / "data" / "strata" / region / f"{region}-census-tracts.parquet")
    if kind == "strata":
        return str(ROOT / "data" / "strata" / region / f"{region}-strata-tract-table.parquet")
    return str(ROOT / "data" / "reference" / region / f"{region}-overture-{kind}.parquet")


def to5070(geoms):
    return shapely.transform(geoms, lambda xy: np.column_stack(_TR.transform(xy[:, 0], xy[:, 1])))


def xy_to5070(x, y):
    return _TR.transform(np.asarray(x), np.asarray(y))


def xy_to4326(x, y):
    return _TR_INV.transform(np.asarray(x), np.asarray(y))


def wkb_to_geoms(series):
    return shapely.from_wkb([bytes(v) for v in series])


def load_tracts(region: str):
    c = con()
    df = c.execute(
        f"SELECT GEOID, pop_total, ST_AsWKB(geometry) w FROM '{p(region, 'tracts')}' ORDER BY GEOID"
    ).fetchdf()
    g = to5070(wkb_to_geoms(df.w))
    df = df.drop(columns="w")
    df["GEOID"] = df["GEOID"].astype(str)
    return df, g


def region_outline(tract_geoms):
    """Region polygon (dissolve of its tracts) and its outline chopped into short pieces for STRtree queries."""
    try:
        reg = shapely.coverage_union_all(tract_geoms)
    except Exception:  # noqa: BLE001
        reg = shapely.union_all(tract_geoms)
    reg = shapely.make_valid(reg)
    rings = shapely.get_parts(shapely.boundary(reg))
    pieces = []
    for r in rings:
        co = shapely.get_coordinates(r)
        for i in range(0, len(co) - 1, 40):
            seg = co[i : min(i + 41, len(co))]
            if len(seg) >= 2:
                pieces.append(shapely.linestrings(seg))
    return reg, np.array(pieces, dtype=object)


def assign_points_to_tracts(x, y, tract_geoms, chunk=1_000_000):
    """Index of containing tract for each point (-1 if none). Deterministic: first (lowest) tract index."""
    tree = shapely.STRtree(tract_geoms)
    out = np.full(len(x), -1, dtype=np.int32)
    for s in range(0, len(x), chunk):
        pts = shapely.points(x[s : s + chunk], y[s : s + chunk])
        ii, jj = tree.query(pts, predicate="intersects")
        order = np.lexsort((jj, ii))
        ii, jj = ii[order], jj[order]
        first = np.r_[True, ii[1:] != ii[:-1]]
        out[s + ii[first]] = jj[first]
    return out


def components_from_connectors(seg_idx: np.ndarray, conn_codes: np.ndarray, n_seg: int):
    """Connected components of segments, where two segments are linked if they share a connector id.
    seg_idx / conn_codes: one row per (segment, connector) incidence."""
    n_conn = int(conn_codes.max()) + 1 if len(conn_codes) else 0
    # bipartite graph segments (0..n_seg-1) + connectors (n_seg..)
    rows = seg_idx
    cols = conn_codes + n_seg
    n = n_seg + n_conn
    A = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n)).tocsr()
    ncomp, lab = connected_components(A, directed=False)
    return lab[:n_seg]


def graph_codes(c, region, layer, classes, tag):
    """Create temp tables {tag}_s(id, class, si) and return incidence arrays (si, ci, at, n_seg)."""
    where = ""
    if classes is not None:
        where = "WHERE class IN (" + ",".join(f"'{k}'" for k in classes) + ")"
    src = p(region, layer)
    c.execute(f"CREATE OR REPLACE TEMP TABLE {tag}_s AS SELECT id, class, (row_number() OVER (ORDER BY id) - 1)::INTEGER si "
              f"FROM '{src}' {where}")
    c.execute(f"CREATE OR REPLACE TEMP TABLE {tag}_inc AS SELECT id, u.connector_id cid, u['at'] AS pos FROM "
              f"(SELECT id, unnest(connectors) u FROM '{src}' {where})")
    c.execute(f"CREATE OR REPLACE TEMP TABLE {tag}_c AS SELECT cid, (row_number() OVER (ORDER BY cid) - 1)::INTEGER ci "
              f"FROM (SELECT DISTINCT cid FROM {tag}_inc)")
    a = c.execute(f"SELECT s.si, cc.ci, i.pos FROM {tag}_inc i JOIN {tag}_s s USING(id) JOIN {tag}_c cc USING(cid) "
                  f"ORDER BY s.si, i.pos, cc.ci").fetchnumpy()
    n = c.execute(f"SELECT count(*) FROM {tag}_s").fetchone()[0]
    return a["si"].astype(np.int64), a["ci"].astype(np.int64), a["pos"].astype(np.float64), int(n)


def edge_flags(geoms, region_poly, pieces):
    tree = shapely.STRtree(geoms)
    _, jj = tree.query(pieces, predicate="dwithin", distance=EDGE_BUF_M)
    near = np.zeros(len(geoms), bool)
    near[np.unique(jj)] = True
    mid = shapely.line_interpolate_point(geoms, 0.5, normalized=True)
    shapely.prepare(region_poly)
    inside = shapely.contains_xy(region_poly, shapely.get_x(mid), shapely.get_y(mid))
    return near | ~inside, tree


def nearest_locality(c, region, x, y, max_m=8000.0):
    """Locality string (addresses[1].locality) of the nearest Overture POI to each (x, y) in EPSG:5070."""
    pois = c.execute(f"SELECT ST_X(geometry) x, ST_Y(geometry) y, addresses[1].locality AS locality FROM '{p(region, 'pois')}' "
                     "WHERE addresses[1].locality IS NOT NULL ORDER BY id").fetchdf()
    px, py = xy_to5070(pois.x.values, pois.y.values)
    tree = shapely.STRtree(shapely.points(px, py))
    (q, t), d = tree.query_nearest(shapely.points(np.asarray(x), np.asarray(y)), max_distance=max_m, return_distance=True)
    out = pd.Series(pois["locality"].values[t]).groupby(q).first().reindex(np.arange(len(x)))
    return out.values


def find_bridges(n_nodes: int, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Iterative Tarjan bridge finding on an undirected multigraph given by edge lists (parallel edges allowed,
    self-loops ignored). Returns boolean mask over edges."""
    m = len(u)
    src = np.concatenate([u, v])
    dst = np.concatenate([v, u])
    eid = np.concatenate([np.arange(m), np.arange(m)])
    keep = src != dst
    src, dst, eid = src[keep], dst[keep], eid[keep]
    order = np.argsort(src, kind="stable")
    dst = dst[order].tolist()
    eid = eid[order].tolist()
    indptr = np.searchsorted(src[order], np.arange(n_nodes + 1)).tolist()
    disc = [-1] * n_nodes
    low = [0] * n_nodes
    is_bridge = np.zeros(m, bool)
    t = 0
    for root in range(n_nodes):
        if disc[root] != -1 or indptr[root] == indptr[root + 1]:
            continue
        disc[root] = low[root] = t
        t += 1
        stack = [(root, -1, indptr[root])]
        while stack:
            node, pe, ptr = stack[-1]
            if ptr < indptr[node + 1]:
                stack[-1] = (node, pe, ptr + 1)
                w = dst[ptr]
                e = eid[ptr]
                if e == pe:
                    continue
                if disc[w] == -1:
                    disc[w] = low[w] = t
                    t += 1
                    stack.append((w, e, indptr[w]))
                elif disc[w] < low[node]:
                    low[node] = disc[w]
            else:
                stack.pop()
                if stack:
                    parent = stack[-1][0]
                    if low[node] < low[parent]:
                        low[parent] = low[node]
                    if low[node] > disc[parent]:
                        is_bridge[pe] = True
    return is_bridge
