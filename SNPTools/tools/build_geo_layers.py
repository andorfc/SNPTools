#!/usr/bin/env python3
"""build_geo_layers.py -- SNPGeo base-map layers from Natural Earth (public domain).

    python3 tools/build_geo_layers.py --ne0 ne_110m_admin_0_countries.geojson \
        --out data/geo/countries.geo.json \
        [--ne1 ne_10m_admin_1_states_provinces.geojson --out1 data/geo/admin1_na.geo.json]

Admin-1 layer (optional, --ne1): Natural Earth 10 m admin-1 restricted to
USA, CAN and MEX (--admin1-countries), each ring simplified with
Douglas-Peucker (--tol degrees, default 0.04), rings whose bounding box is
smaller than --min-extent degrees dropped (small islands), Aleutian rings east
of the antimeridian dropped, coordinates rounded to 2 dp; properties
{iso3, code, name} where code = ISO 3166-2 (iso_3166_2, e.g. US-IA), the key
SNPGEO_REGIONS[...].admin1Code carries.

Recipe (identical to fusarium build_snpgeo_data.py / build_maize_snpgeo_data.py):
Natural Earth 110 m admin-0, coordinates rounded to 2 decimal places,
consecutive duplicate points dropped, Antarctica dropped, properties reduced to
{iso3, name}. iso3 = first of ISO_A3_EH, ISO_A3, ADM0_A3, SOV_A3 that is not -99.

Pinned source: nvkelso/natural-earth-vector tag v5.1.2
(commit f1890d9f152c896d250a77557a5751a93d494776), geojson/ directory.
See data/geo/PROVENANCE.md. Python 3.8+ standard library only.
"""
import argparse, json, os


def iso3_of(props):
    for k in ("ISO_A3_EH", "ISO_A3", "ADM0_A3", "SOV_A3"):
        v = (props.get(k) or "").strip()
        if v and v != "-99":
            return v
    return None


def round_ring(ring, nd):
    out = []
    for pt in ring:
        p = [round(pt[0], nd), round(pt[1], nd)]
        if not out or out[-1] != p:
            out.append(p)
    return out if len(out) >= 4 else None


def simplify_geom(geom, nd, ring_fn=None):
    ring_fn = ring_fn or (lambda r: round_ring(r, nd))
    t = geom.get("type")
    if t == "Polygon":
        rings = [r for r in (ring_fn(r) for r in geom["coordinates"]) if r]
        return {"type": "Polygon", "coordinates": rings} if rings else None
    if t == "MultiPolygon":
        polys = []
        for poly in geom["coordinates"]:
            rings = [r for r in (ring_fn(r) for r in poly) if r]
            if rings:
                polys.append(rings)
        return {"type": "MultiPolygon", "coordinates": polys} if polys else None
    return None


def write_fc(feats, out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump({"type": "FeatureCollection", "features": feats}, fh, separators=(",", ":"))


def build_countries(ne, out_path, nd=2, drop=("ATA",)):
    feats, skipped = [], []
    for f in ne["features"]:
        p = f["properties"]
        iso = iso3_of(p)
        name = (p.get("NAME") or p.get("ADMIN") or "").strip()
        if not iso or iso in drop:
            skipped.append(name or iso); continue
        g = simplify_geom(f.get("geometry") or {}, nd)
        if not g:
            skipped.append(name); continue
        feats.append({"type": "Feature", "properties": {"iso3": iso, "name": name}, "geometry": g})
    write_fc(feats, out_path)
    return feats, skipped


def _perp(p, a, b):
    (x, y), (x1, y1), (x2, y2) = p, a, b
    dx, dy = x2 - x1, y2 - y1
    if dx == 0 and dy == 0:
        return ((x - x1) ** 2 + (y - y1) ** 2) ** 0.5
    t = max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / (dx * dx + dy * dy)))
    return ((x - x1 - t * dx) ** 2 + (y - y1 - t * dy) ** 2) ** 0.5


def douglas_peucker(pts, tol):
    """Iterative Douglas-Peucker on an open polyline; keeps both end points."""
    if len(pts) < 3:
        return pts[:]
    keep = [False] * len(pts); keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        i, j = stack.pop()
        dmax, idx = 0.0, None
        for k in range(i + 1, j):
            d = _perp(pts[k], pts[i], pts[j])
            if d > dmax:
                dmax, idx = d, k
        if idx is not None and dmax > tol:
            keep[idx] = True
            stack += [(i, idx), (idx, j)]
    return [p for p, k in zip(pts, keep) if k]


def simplify_ring(ring, tol, nd, min_extent, drop_east_of=None):
    xs = [p[0] for p in ring]; ys = [p[1] for p in ring]
    if max(xs) - min(xs) < min_extent and max(ys) - min(ys) < min_extent:
        return None
    if drop_east_of is not None and min(xs) > drop_east_of:
        return None
    closed = ring[:-1] if ring[0] == ring[-1] else ring[:]
    if len(closed) < 3:
        return None
    # split the closed ring at its farthest point so DP keeps its shape
    far = max(range(len(closed)), key=lambda k: (closed[k][0] - closed[0][0]) ** 2 + (closed[k][1] - closed[0][1]) ** 2)
    a = douglas_peucker(closed[:far + 1], tol)
    b = douglas_peucker(closed[far:] + [closed[0]], tol)
    out = round_ring(a[:-1] + b, nd)
    if out and out[0] != out[-1]:
        out.append(out[0])
    return out if out and len(out) >= 4 else None


def build_admin1(ne1, out_path, countries=("USA", "CAN", "MEX"), tol=0.04, nd=2, min_extent=0.15):
    feats, skipped = [], []
    for f in ne1["features"]:
        p = f["properties"]
        iso = (p.get("adm0_a3") or "").strip()
        if iso not in countries:
            continue
        code = (p.get("iso_3166_2") or "").strip()
        name = (p.get("name") or p.get("name_en") or "").strip()
        east = 0.0 if iso == "USA" else None      # Aleutians west of 180 deg would wrap
        g = simplify_geom(f.get("geometry") or {}, nd,
                          ring_fn=lambda r: simplify_ring(r, tol, nd, min_extent, east))
        if not g or not code:
            skipped.append(code or name); continue
        feats.append({"type": "Feature", "properties": {"iso3": iso, "code": code, "name": name}, "geometry": g})
    write_fc(feats, out_path)
    return feats, skipped


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ne0", required=True, help="ne_110m_admin_0_countries.geojson")
    ap.add_argument("--out", default="data/geo/countries.geo.json")
    ap.add_argument("--ne1", help="ne_10m_admin_1_states_provinces.geojson (optional admin-1 layer)")
    ap.add_argument("--out1", default="data/geo/admin1_na.geo.json")
    ap.add_argument("--admin1-countries", default="USA,CAN,MEX")
    ap.add_argument("--tol", type=float, default=0.04)
    ap.add_argument("--min-extent", type=float, default=0.15)
    a = ap.parse_args()
    feats, skipped = build_countries(json.load(open(a.ne0, encoding="utf-8")), a.out)
    print(f"{a.out}: {len(feats)} features; skipped {skipped}")
    if a.ne1:
        f1, s1 = build_admin1(json.load(open(a.ne1, encoding="utf-8")), a.out1,
                              tuple(a.admin1_countries.split(",")), a.tol, 2, a.min_extent)
        print(f"{a.out1}: {len(f1)} features; skipped {s1}")


if __name__ == "__main__":
    main()
