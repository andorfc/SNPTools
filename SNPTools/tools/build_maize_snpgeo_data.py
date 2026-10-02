#!/usr/bin/env python3
"""
build_maize_snpgeo_data.py -- SNPGeo assets for the maize GRIN-linked v1 set.
Port of fusarium SNPTools/build_snpgeo_data.py; output formats are unchanged.

Outputs (relative to --outdir, mirroring SNPTools/)
---------------------------------------------------
  data/geo/countries.geo.json   Natural Earth admin-0 polygons, trimmed + rounded,
                                properties {iso3, name}  (identical recipe to fusarium)
  js/snpgeo.regions.js          window.SNPGEO_COUNTRY_ISO  country spelling -> ISO3
                                window.SNPGEO_REGIONS      sample_id -> geo record

Record keys (the fusarium JS reads iso3, country, state, county; the rest ride along):
  iso3, country        ISO3 from GRIN passport (T5 country_iso3), validated against
                       Natural Earth; country = Natural Earth display name
  state                GRIN admin-1 name (T5 admin1_name)
  admin1Code           ISO 3166-2 code (T5 admin1_iso3166_2), for a later admin-1 map
  location             GRIN collection-site text, when present
  locality             GRIN origin statement verbatim ("Developed - Iowa, United States")
  originType           Developed / Donated / Collected
  lat, lon             ONLY GRIN passport coordinates, and only when both parse as
                       in-range decimal degrees (fusarium rule: never invent a point).
                       The admin-1 / country label points in T5 (origin_*_derived)
                       are NOT emitted as coordinates.
  grin, acqYear        GRIN accession number; year of GRIN acquisition (not a
                       collection year, hence not the fusarium `year` key)
  src                  't5' (GRIN passport) or 'crosswalk' (class C: GRIN origin text
                       recorded in the Grzybowski crosswalk; no T5 row)

Usage:
  python3 build_maize_snpgeo_data.py --ne ne_110m_admin_0_countries.geojson \
      --samples data/zmgrin2026_samples.tsv --t5 T5.tsv --crosswalk crosswalk.tsv \
      --outdir out/SNPTools
"""
import argparse, csv, json, os, re, sys
from collections import Counter, OrderedDict

# ---------------------------------------------------------------------------
# Geometry helpers: copied verbatim in behaviour from fusarium build_snpgeo_data.py
# ---------------------------------------------------------------------------
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


def simplify_geom(geom, nd):
    t = geom.get("type")
    if t == "Polygon":
        rings = [r for r in (round_ring(r, nd) for r in geom["coordinates"]) if r]
        return {"type": "Polygon", "coordinates": rings} if rings else None
    if t == "MultiPolygon":
        polys = []
        for poly in geom["coordinates"]:
            rings = [r for r in (round_ring(r, nd) for r in poly) if r]
            if rings:
                polys.append(rings)
        return {"type": "MultiPolygon", "coordinates": polys} if polys else None
    return None


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
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump({"type": "FeatureCollection", "features": feats}, fh, separators=(",", ":"))
    return feats, skipped


def ne_name_index(ne):
    idx = {}
    for f in ne["features"]:
        p = f["properties"]
        iso = iso3_of(p)
        if not iso:
            continue
        for k in ("NAME", "ADMIN", "NAME_LONG", "NAME_EN", "GEOUNIT"):
            v = (p.get(k) or "").strip().lower()
            if v:
                idx.setdefault(v, iso)
    return idx


# GRIN origin spellings that do not match Natural Earth verbatim (maize-specific).
COUNTRY_ALIAS = {
    "united states": "united states of america",
    "virgin islands (british)": "british virgin islands",
    "soviet union": None,          # historical; no polygon  (T5 carries iso3 '108' here)
    "former soviet union": None,
    "unknown": None,
}
ISO3_OK = re.compile(r"^[A-Z]{3}$")
ORIGIN_PREFIX = re.compile(r"^\s*(Developed|Donated|Collected|Unknown)\s*[\u2013-]\s*", re.I)


def resolve_country(name, idx):
    key = (name or "").strip().lower()
    if not key:
        return None
    if key in COUNTRY_ALIAS:
        if COUNTRY_ALIAS[key] is None:
            return None
        key = COUNTRY_ALIAS[key]
    return idx.get(key)


def origin_parts(origin):
    o = origin or ""
    m = ORIGIN_PREFIX.match(o)
    kind = m.group(1).capitalize() if m else ""
    parts = [x.strip() for x in ORIGIN_PREFIX.sub("", o).split(",") if x.strip()]
    country = parts[-1] if parts else ""
    admin1 = ", ".join(parts[:-1]) if len(parts) > 1 else ""
    if country.lower() == "soviet union" and admin1.lower() == "former":
        country, admin1 = "Former Soviet Union", ""
    return kind, admin1, country


BLANK = {".", "-", "--", "na", "n/a", "nan", "none", "null", "unknown", ""}


def clean(v):
    s = ("" if v is None else str(v)).replace("\u00a0", " ").strip()
    return "" if s.lower() in BLANK else s


def dec_degrees(v, lim):
    s = clean(v)
    if not s or not re.fullmatch(r"-?\d{1,3}(?:\.\d+)?", s):
        return None
    x = float(s)
    return x if -lim <= x <= lim else None


def read_tsv(p):
    with open(p, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--t5", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--catalog", help="js/zmgrin.catalog.js, for join validation")
    ap.add_argument("--precision", type=int, default=2)
    a = ap.parse_args()

    ne = json.load(open(a.ne, encoding="utf-8"))
    geo_out = os.path.join(a.outdir, "data", "geo", "countries.geo.json")
    feats, skipped = build_countries(ne, geo_out, nd=a.precision)
    idx = ne_name_index(ne)
    iso_to_name = {f["properties"]["iso3"]: f["properties"]["name"] for f in feats}
    print(f"[countries] {len(feats)} features -> {geo_out} ({os.path.getsize(geo_out)/1024:.0f} KB, "
          f"{a.precision} dp); skipped: {', '.join(sorted(map(str, skipped)))}")

    samples = read_tsv(a.samples)
    t5 = {r["accession_number"]: r for r in read_tsv(a.t5)}
    xw = {r["VCFname"]: r for r in read_tsv(a.crosswalk)}

    regions, spellings = OrderedDict(), {}
    stats = Counter()
    nopoly, bad_iso = Counter(), Counter()
    for s in samples:
        sid, grin = s["sample_id"], s["grin_accession"]
        p5 = t5.get(grin)
        rec = OrderedDict()
        if p5:
            rec["src"] = "t5"
            kind, adm1_txt, ctry_txt = origin_parts(p5.get("origin"))
            iso = clean(p5.get("country_iso3")).upper()
            if iso and not ISO3_OK.match(iso):
                bad_iso[f"{iso} ({ctry_txt})"] += 1
                iso = resolve_country(ctry_txt, idx) or ""
            state = clean(p5.get("admin1_name")) or adm1_txt
            adm1code = clean(p5.get("admin1_iso3166_2"))
            location = clean(p5.get("collection_site_grin_export"))
            lat = dec_degrees(p5.get("passport_latitude"), 90)
            lon = dec_degrees(p5.get("passport_longitude"), 180)
            acq = clean(p5.get("acquisition_date"))[:4]
            locality = clean(p5.get("origin"))
        else:
            rec["src"] = "crosswalk"
            xr = xw.get(s["vcf_name"], {})
            kind, state, ctry_txt = origin_parts(xr.get("grin_origin", ""))
            iso = resolve_country(ctry_txt, idx) or ""
            adm1code = location = acq = ""
            lat = lon = None
            locality = clean(xr.get("grin_origin"))
        if iso:
            rec["iso3"] = iso
            rec["country"] = iso_to_name.get(iso) or ctry_txt or iso
            if iso not in iso_to_name:
                nopoly[f"{iso} ({ctry_txt})"] += 1
            if ctry_txt:
                spellings[ctry_txt] = iso
        elif ctry_txt:
            rec["country"] = ctry_txt
            stats["no_iso3"] += 1
        if state: rec["state"] = state
        if adm1code: rec["admin1Code"] = adm1code
        if location: rec["location"] = location
        if kind: rec["originType"] = kind
        if locality: rec["locality"] = locality
        if lat is not None and lon is not None:
            rec["lat"], rec["lon"] = lat, lon
            stats["coords"] += 1
        if grin: rec["grin"] = grin
        if acq.isdigit(): rec["acqYear"] = acq
        regions[sid] = rec
        stats[rec["src"]] += 1

    js_out = os.path.join(a.outdir, "js", "snpgeo.regions.js")
    os.makedirs(os.path.dirname(js_out), exist_ok=True)
    with open(js_out, "w", encoding="utf-8") as fh:
        fh.write("/* AUTO-GENERATED by build_maize_snpgeo_data.py -- do not edit by hand.\n"
                 "   Per-sample geographic metadata for SNPGeo (maize GRIN-linked v1), keyed by the\n"
                 "   release VCF sample column (sample_id, e.g. ZmG_B73). Same format as the fusarium\n"
                 "   snpgeo.regions.js: iso3, country (Natural Earth display name), state (GRIN admin-1),\n"
                 "   plus admin1Code, location, originType, locality (GRIN origin text), lat/lon (GRIN\n"
                 "   passport coordinates only), grin, acqYear, src. No coordinates are inferred. */\n")
        fh.write("window.SNPGEO_COUNTRY_ISO = ")
        json.dump(spellings, fh, separators=(",", ":"), sort_keys=True, ensure_ascii=False)
        fh.write(";\nwindow.SNPGEO_REGIONS = ")
        json.dump(regions, fh, separators=(",", ":"), sort_keys=True, ensure_ascii=False)
        fh.write(";\n")

    n = len(regions)
    with_iso = sum(1 for r in regions.values() if "iso3" in r)
    print(f"[regions]   {n} samples -> {js_out} ({os.path.getsize(js_out)/1024:.0f} KB); "
          f"src t5={stats['t5']} crosswalk={stats['crosswalk']}")
    print(f"[origin]    {with_iso}/{n} samples resolved to ISO3; "
          f"{sum(1 for r in regions.values() if r.get('iso3') in iso_to_name)}/{n} have a polygon")
    print(f"[coords]    {stats['coords']} samples with GRIN passport lat/lon")
    print(f"[state]     {sum(1 for r in regions.values() if 'state' in r)}/{n} with admin-1 name; "
          f"{sum(1 for r in regions.values() if 'admin1Code' in r)}/{n} with ISO 3166-2 code")
    print("[countries] samples per ISO3: " + ", ".join(f"{k}={v}" for k, v in
          Counter(r.get("iso3", "none") for r in regions.values()).most_common()))
    if bad_iso:
        print(f"!! non-ISO3 country codes in T5 (resolved from origin text instead): {dict(bad_iso)}")
    if nopoly:
        print(f"!! ISO3 codes with no Natural Earth 110m polygon (listed in tables, not on the map): {dict(nopoly)}")
    if stats["no_iso3"]:
        print(f"!! {stats['no_iso3']} samples with a country name but no ISO3 (historical/unmappable)")

    if a.catalog and os.path.exists(a.catalog):
        src = open(a.catalog, encoding="utf-8").read()
        m = re.search(r"window\.SNP_CATALOG\.families\[\"([^\"]+)\"\]\s*=\s*", src)
        fam, _ = json.JSONDecoder().raw_decode(src[m.end():])
        ids = [x["id"] for p in fam["projects"] for g in p["groups"] for x in g["accessions"]]
        hit = set(ids) & set(regions)
        print(f"[join]      {len(hit)}/{len(ids)} catalogue samples have a region record; "
              f"{len(set(regions) - set(ids))} region records not in the catalogue")


if __name__ == "__main__":
    main()
