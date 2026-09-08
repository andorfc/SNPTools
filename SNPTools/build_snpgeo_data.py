#!/usr/bin/env python3
"""
build_snpgeo_data.py — generate the static assets SNPGeo needs.

Outputs
-------
  data/geo/countries.geo.json   world country polygons, trimmed + rounded.
                                properties: {iso3, name}. Shipped in-repo so the
                                map has no CDN dependency at runtime.
  js/snpgeo.regions.js          window.SNPGEO_REGIONS      isolate id -> sub-country
                                                            metadata (state/county/...)
                                window.SNPGEO_COUNTRY_ISO  country name -> ISO3

Country of origin is NOT emitted here: it already ships in js/strains.catalog.js
as window.FUSARIUM_STRAINS[].country, which stays the single source of truth.
This script only adds the sub-country detail the catalogue lacks.

Source-data notes (why this file is defensive)
----------------------------------------------
* "GPS1 /Latitude" is not reliably a latitude. Of 211 non-empty values only 82
  parse as decimal degrees; the rest are free-text localities or macro-regions
  ("Far East", "Poitou Charentes", "Cook County, Minnesota") and a handful of
  ragged DMS strings ("08°01'10''", "8'36'5''", "08\"42\"26") with no hemisphere.
  A coordinate is therefore emitted ONLY when latitude and longitude BOTH parse
  as in-range decimal degrees; every other value is kept verbatim as `locality`
  text and surfaced in the region table. Nothing is inferred or reprojected.
* Natural Earth stores ISO_A3 = "-99" for several sovereign states (France,
  Norway); ISO_A3_EH carries the real code. The name index covers NAME and ADMIN
  so both "Serbia" and "Republic of Serbia" resolve.
* A few isolate ids repeat across metadata rows (replicate entries). Duplicates
  are MERGED field-by-field preferring non-empty values, rather than last-wins.
* Catalogue ids are zero-padded inconsistently vs the metadata ("06163" vs
  "6163"), so each record is also keyed under a zero-stripped alias.

Usage
-----
  python3 build_snpgeo_data.py \
      --ne /tmp/ne110m.geojson \
      --meta ~/Desktop/snpfold_validation/manual_geographic_snptools_sept3_1229pm.txt
"""
import argparse, csv, json, os, re, sys
from collections import Counter, OrderedDict, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ID_COL = "New SNPversity ID"


# ----------------------------------------------------------------- geometry --
def iso3_of(props):
    """Natural Earth uses "-99" for some sovereign states; ISO_A3_EH is real."""
    for k in ("ISO_A3_EH", "ISO_A3", "ADM0_A3", "SOV_A3"):
        v = (props.get(k) or "").strip()
        if v and v != "-99":
            return v
    return None


def round_ring(ring, nd):
    """Round a coordinate ring, dropping consecutive duplicates it creates."""
    out = []
    for pt in ring:
        p = [round(pt[0], nd), round(pt[1], nd)]
        if not out or out[-1] != p:
            out.append(p)
    return out if len(out) >= 4 else None      # a closed ring needs >= 4 points


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


def build_countries(ne_path, out_path, nd=2, drop=("ATA",)):
    ne = json.load(open(ne_path, encoding="utf-8"))
    feats, skipped = [], []
    for f in ne["features"]:
        p = f["properties"]
        iso = iso3_of(p)
        name = (p.get("NAME") or p.get("ADMIN") or "").strip()
        if not iso or iso in drop:
            skipped.append(name or iso)
            continue
        g = simplify_geom(f.get("geometry") or {}, nd)
        if not g:
            skipped.append(name)
            continue
        feats.append({"type": "Feature",
                      "properties": {"iso3": iso, "name": name},
                      "geometry": g})
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump({"type": "FeatureCollection", "features": feats}, fh,
                  separators=(",", ":"))
    return feats, skipped


def ne_name_index(ne_path):
    """Index NAME *and* ADMIN so 'Serbia' and 'Republic of Serbia' both hit."""
    ne = json.load(open(ne_path, encoding="utf-8"))
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


# Metadata country spellings that do NOT match Natural Earth verbatim.
COUNTRY_ALIAS = {
    "usa": "united states of america",
    "united states": "united states of america",
    "austrailia": "australia",          # typo in the source file
    "czech republic": "czechia",
    "unkown": None,                      # typo of "unknown" -> no geometry
    "unknown": None,                     # explicit unknown origin -> no geometry
}


def resolve_country(name, idx):
    c = (name or "").strip()
    if not c:
        return None
    key = c.lower()
    if key in COUNTRY_ALIAS:
        target = COUNTRY_ALIAS[key]
        if target is None:
            return None                  # deliberately unmapped
        key = target
    return idx.get(key)


# --------------------------------------------------------------- field prep --
BLANK = {".", "-", "--", "na", "n/a", "nan", "none", "null", "unknown", "unkown", ""}


def clean(v):
    if v is None:
        return ""
    s = str(v).replace("\u00a0", " ").strip()
    return "" if s.lower() in BLANK else s


def clean_region(v):
    """Like clean(), but for REGION_COLS (state/county/location/...) only:
    an explicit "Unknown"/"unkown" entry is preserved as visible text
    ("Unknown", canonicalized) instead of being collapsed to blank, so a
    deliberately-labeled record is distinguishable from a genuinely empty
    cell in the region-breakdown table. True placeholder tokens (".", "-",
    "n/a", "", ...) still collapse to blank. Country's "unknown -> no
    geometry" handling is untouched -- that logic lives in COUNTRY_ALIAS /
    resolve_country and never calls this function.
    """
    if v is None:
        return ""
    s = str(v).replace("\u00a0", " ").strip()
    low = s.lower()
    if low in ("unknown", "unkown"):
        return "Unknown"
    return "" if low in BLANK else s


def dec_degrees(v):
    """Decimal degrees only. Returns None for DMS / free text / out-of-range."""
    s = clean(v)
    if not s or not re.fullmatch(r"-?\d{1,3}(?:\.\d+)?", s):
        return None
    x = float(s)
    return x if -180.0 <= x <= 180.0 else None


def id_aliases(iid):
    """('06163',) -> {'06163','6163'} so padded/unpadded ids join."""
    out = {iid}
    stripped = iid.lstrip("0")
    if stripped and stripped != iid:
        out.add(stripped)
    return out


# Sub-country columns lifted from the metadata TSV -> compact JS keys.
# "GPS1 /Latitude" / "GPS2/Longitude" are handled separately (see dec_degrees).
REGION_COLS = OrderedDict([
    ("State",             "state"),
    ("County",            "county"),
    ("Location",          "location"),
    ("Year",              "year"),
    ("Cultivar",          "cultivar"),
    ("Host Crop;Variety", "hostvar"),
    ("Additional note",   "note"),
])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ne", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--strains", default=os.path.join(HERE, "js", "strains.catalog.js"))
    ap.add_argument("--precision", type=int, default=2)
    args = ap.parse_args()

    # ---- geometry ----------------------------------------------------------
    geo_out = os.path.join(HERE, "data", "geo", "countries.geo.json")
    feats, skipped = build_countries(args.ne, geo_out, nd=args.precision)
    idx = ne_name_index(args.ne)
    print(f"[countries] {len(feats)} features -> {geo_out} "
          f"({os.path.getsize(geo_out)/1024:.0f} KB, {args.precision} dp)")
    if skipped:
        print(f"[countries] intentionally skipped: {', '.join(sorted(map(str, skipped)))}")

    # ---- catalogue (for validation only) -----------------------------------
    strains, strain_ids = None, None
    if os.path.exists(args.strains):
        src = open(args.strains, encoding="utf-8").read()
        m = re.search(r"window\.FUSARIUM_STRAINS\s*=\s*(\[.*?\]);", src, re.S)
        if m:
            strains = json.loads(m.group(1))
            strain_ids = set()
            for s in strains:
                for f in ("id", "baseId"):
                    if s.get(f):
                        strain_ids |= id_aliases(str(s[f]))

    # ---- metadata ----------------------------------------------------------
    with open(args.meta, encoding="utf-8", errors="replace") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows or ID_COL not in rows[0]:
        sys.exit(f"!! metadata is missing the {ID_COL!r} column")

    grouped, blank_id = defaultdict(list), 0
    for r in rows:
        iid = clean(r.get(ID_COL))
        if not iid:
            blank_id += 1
            continue
        grouped[iid].append(r)

    # iso3 -> canonical country display name (from Natural Earth), so typos and
    # abbreviations in the source ("Austrailia", "USA") show a clean label.
    iso_to_name = {f["properties"]["iso3"]: f["properties"]["name"] for f in feats}

    regions, coord_ok, coord_text, merged, no_origin = {}, 0, 0, [], 0
    for iid, group in grouped.items():
        if len(group) > 1:
            merged.append(iid)
        rec = {}
        # country of origin, resolved to ISO3 so the map can shade it directly.
        # The metadata "Country" column is authoritative (user-designated); the
        # strain catalogue, when present, only overrides genuinely blank values.
        raw_country = ""
        for r in group:
            v = clean(r.get("Country"))
            if v:
                raw_country = v
                break
        iso = resolve_country(raw_country, idx)
        if iso:
            rec["iso3"] = iso
            rec["country"] = iso_to_name.get(iso, raw_country)
        else:
            no_origin += 1
            if raw_country:            # a spelling we could not map to a polygon
                rec["country"] = raw_country
        # merge duplicates field-by-field, first non-empty wins
        for src_col, key in REGION_COLS.items():
            for r in group:
                v = clean_region(r.get(src_col))
                if v:
                    rec[key] = v
                    break
        # coordinates: emit only when BOTH axes parse as decimal degrees
        lat = lon = None
        for r in group:
            la, lo = dec_degrees(r.get("GPS1 /Latitude")), dec_degrees(r.get("GPS2/Longitude"))
            if la is not None and lo is not None:
                lat, lon = la, lo
                break
        if lat is not None:
            rec["lat"], rec["lon"] = lat, lon
            coord_ok += 1
        else:
            # keep the unparseable value as locality text — never invent a point
            for r in group:
                v = clean(r.get("GPS1 /Latitude"))
                if v and dec_degrees(v) is None:
                    rec["locality"] = v
                    coord_text += 1
                    break
        regions[iid] = rec

    # ---- optional catalogue override for isolates with no mapped origin -----
    cat_filled = 0
    if strains:
        cat_country = {}
        for s in strains:
            c = clean(s.get("country"))
            if c:
                for a in id_aliases(str(s.get("id", ""))):
                    cat_country[a] = c
                if s.get("baseId"):
                    for a in id_aliases(str(s["baseId"])):
                        cat_country.setdefault(a, c)
        for iid, rec in regions.items():
            if "iso3" not in rec:
                c = next((cat_country[a] for a in id_aliases(iid) if a in cat_country), None)
                iso = resolve_country(c, idx) if c else None
                if iso:
                    rec["iso3"] = iso
                    rec["country"] = iso_to_name.get(iso, c)
                    cat_filled += 1

    # ---- country -> ISO3 ---------------------------------------------------
    names = {clean(r.get("Country")) for r in rows}
    if strains:
        names |= {clean(s.get("country")) for s in strains}
    country_iso, unresolved = {}, []
    for n in sorted(x for x in names if x):
        iso = resolve_country(n, idx)
        if iso:
            country_iso[n] = iso
        elif n.lower() not in COUNTRY_ALIAS:
            unresolved.append(n)

    # ---- emit --------------------------------------------------------------
    js_out = os.path.join(HERE, "js", "snpgeo.regions.js")
    with open(js_out, "w", encoding="utf-8") as fh:
        fh.write("/* AUTO-GENERATED by build_snpgeo_data.py — do not edit by hand.\n"
                 "   Per-isolate geographic metadata for SNPGeo, keyed by SNPversity isolate id.\n"
                 "   Each record may carry: iso3, country (canonical display name), state,\n"
                 "   county, location, year, cultivar, hostvar, note, and either lat/lon\n"
                 "   (both axes parsed as decimal degrees) or `locality` (free text from the\n"
                 "   GPS column that is not a coordinate). Country is resolved from the\n"
                 "   metadata 'Country' column (authoritative); the strain catalogue only\n"
                 "   fills genuinely blank origins. This makes SNPGeo self-sufficient and\n"
                 "   independent of which branch (maize/fusarium) is checked out. */\n")
        fh.write("window.SNPGEO_COUNTRY_ISO = ")
        json.dump(country_iso, fh, separators=(",", ":"), sort_keys=True)
        fh.write(";\n")
        fh.write("window.SNPGEO_REGIONS = ")
        json.dump(regions, fh, separators=(",", ":"), sort_keys=True)
        fh.write(";\n")
    print(f"[regions]   {len(regions)} isolates -> {js_out} "
          f"({os.path.getsize(js_out)/1024:.0f} KB)")
    print(f"[regions]   {len(country_iso)} country spellings resolved to ISO3")
    print(f"[coords]    {coord_ok} usable lat/lon pairs; "
          f"{coord_text} GPS values kept as locality text")
    with_iso = sum(1 for r in regions.values() if "iso3" in r)
    print(f"[origin]    {with_iso}/{len(regions)} isolates mapped to a country polygon"
          + (f" ({cat_filled} via catalogue fallback)" if strains and cat_filled else ""))

    # ---- validation --------------------------------------------------------
    if unresolved:
        print(f"!! UNRESOLVED countries (would render as no-data): {unresolved}")
    if blank_id:
        print(f"!! {blank_id} rows had a blank {ID_COL!r} and were dropped")
    if merged:
        print(f"[regions]   merged duplicate rows for: {sorted(merged)}")

    if strain_ids is not None:
        hit = [i for i in regions if id_aliases(i) & strain_ids]
        print(f"[join]      {len(hit)}/{len(regions)} metadata ids match the catalogue")
        extra = sorted(i for i in regions if not (id_aliases(i) & strain_ids))
        if extra:
            print(f"!! in metadata but not in catalogue: {extra}")
        region_keys = set()
        for i in regions:
            region_keys |= id_aliases(i)
        nogeo = sorted({str(s["id"]) for s in strains
                        if not (id_aliases(str(s["id"])) & region_keys)
                        and not (s.get("baseId") and id_aliases(str(s["baseId"])) & region_keys)})
        if nogeo:
            print(f"!! catalogue isolates with no sub-country metadata ({len(nogeo)}): {nogeo}")
        filled = Counter()
        for i in hit:
            filled.update(regions[i].keys())
        print("[join]      field coverage among matched isolates:")
        for k in list(REGION_COLS.values()) + ["locality", "lat"]:
            print(f"              {k:>9}: {filled.get(k, 0)}/{len(hit)}")


if __name__ == "__main__":
    main()
