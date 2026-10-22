"""Overview maps, one per region: share of residents in single-egress pockets per tract, with WHP-Q4 tracts outlined.

Inputs: out/tract_metrics_all.parquet (from 03_stratify) and the challenge tract polygons
(data/strata/<region>/<region>-census-tracts.parquet). Output: figs/region_maps.png
Usage: python -I discovery/13_region_maps.py
"""
from pathlib import Path

import duckdb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import PatchCollection  # noqa: E402
from matplotlib.patches import Polygon as MplPolygon  # noqa: E402
from shapely import wkb  # noqa: E402

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "strata"
REGIONS = ["northern-ca", "eastern-wa", "maricopa-az", "eastern-ok", "south-central-tx"]
BINS = [0, 0.001, 0.05, 0.20, 0.50, 1.01]
COLORS = ["#f2f2f2", "#fcd9b6", "#f59f63", "#d6542a", "#8c1d0c"]
LABELS = ["0%", "0–5%", "5–20%", "20–50%", "≥ 50%"]


def patches(geom):
    polys = getattr(geom, "geoms", [geom])
    return [MplPolygon(np.asarray(p.exterior.coords), closed=True) for p in polys if p.geom_type == "Polygon"]


def main():
    m = pd.read_parquet(HERE / "out" / "tract_metrics_all.parquet")
    m["share"] = (m["est_pop_pocket50"] / m["pop_total"].clip(lower=1)).clip(0, 1)
    con = duckdb.connect()
    con.sql("LOAD spatial;")
    fig, axes = plt.subplots(2, 3, figsize=(11, 7.2), dpi=150)
    for ax, r in zip(axes.flat, REGIONS):
        p = (DATA / r / f"{r}-census-tracts.parquet").as_posix()
        g = con.sql(f"SELECT GEOID, ST_AsWKB(geometry::GEOMETRY) AS w FROM '{p}'").df()
        g = g.merge(m[["GEOID", "share", "whp_q"]], on="GEOID", how="left")
        for q4 in (False, True):
            sub = g[(g["whp_q"] == "Q4") == q4]
            ps, cs = [], []
            for row in sub.itertuples():
                k = int(np.digitize(0 if pd.isna(row.share) else row.share, BINS) - 1)
                for pt in patches(wkb.loads(bytes(row.w))):
                    ps.append(pt)
                    cs.append(COLORS[max(0, min(k, 4))])
            ax.add_collection(PatchCollection(ps, facecolor=cs, edgecolor="#222" if q4 else "#bbb",
                                              linewidths=0.6 if q4 else 0.15, zorder=2 if q4 else 1))
        ax.autoscale_view()
        ax.set_aspect("equal")
        ax.set_title(r, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    ax = axes.flat[-1]
    ax.axis("off")
    for c, lab in zip(COLORS, LABELS):
        ax.add_patch(plt.Rectangle((0, 0), 0, 0, fc=c, ec="#999", label=lab))
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, fc="white", ec="#222", lw=1.2, label="WHP Q4 tract (outlined)"))
    ax.legend(loc="center", frameon=False, fontsize=9, title="residents in single-egress pockets", title_fontsize=9)
    fig.suptitle("Where single-egress pockets are: share of tract residents behind one road segment", fontsize=11)
    fig.tight_layout()
    out = HERE / "figs" / "region_maps.png"
    fig.savefig(out)
    print("->", out)


if __name__ == "__main__":
    main()
