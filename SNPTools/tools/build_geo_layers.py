#!/usr/bin/env python3
"""build_geo_layers.py -- SNPGeo base-map layers from Natural Earth (public domain).

    python3 tools/build_geo_layers.py --ne0 ne_110m_admin_0_countries.geojson \
        --out data/geo/countries.geo.json

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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ne0", required=True, help="ne_110m_admin_0_countries.geojson")
    ap.add_argument("--out", default="data/geo/countries.geo.json")
    a = ap.parse_args()
    feats, skipped = build_countries(json.load(open(a.ne0, encoding="utf-8")), a.out)
    print(f"{a.out}: {len(feats)} features; skipped {skipped}")


if __name__ == "__main__":
    main()
