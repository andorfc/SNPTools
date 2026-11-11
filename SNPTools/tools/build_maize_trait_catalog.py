#!/usr/bin/env python3
"""
build_maize_trait_catalog.py -- SNPTrait data for the maize GRIN-linked release (v1.4).

Outputs (paths relative to --outdir, which mirrors the SNPTools/ tree)
------------------------------------------------------------------------
  js/zmgrin.catalog.js
      * window.SNP_CATALOG.families[<family>]   -- merged INTO the catalog that
        js/accessions.catalog.js already defines (load this file after it and
        before data.js). Same families -> projects -> groups -> accessions shape
        as fusarium js/strains.catalog.js; each accession carries the base keys
        (id, run, founder, label, rep, reps, namFounder) plus the flat, scalar
        metadata strings SNPTrait facets on.
      * window.SNPTRAIT_SCHEMA[<family>]        -- a schema object with exactly
        the shape of TRAIT_SCHEMA.graminearum / .vert in fusarium js/snptrait.js.
  data/traits/<family>.traits.json
      Full per-sample GRIN trait summary (all T4b traits, all scale types) plus
      the trait dictionary, for a later numeric-range facet / trait column view.
      Not read by the current fusarium JS.

Sources: sample manifest (build_maize_samples.py), T5 (passport), T4b (trait
summary by accession), T4c (trait dictionary), crosswalk (Grzybowski
Subpopulation / ImprovementEra for every VCF column, incl. class C).

Faceting rules (the fusarium SNPTrait matches facets by exact string equality):
  * one scalar string per facet per sample; missing -> key omitted (JS shows 'Unknown')
  * panel membership is multi-valued (B73 is NAM + Ames282 + WiDiv), so it is
    exposed as three yes/no facets plus a single-valued `panel` used for grouping
    (precedence NAM > Ames282 > WiDiv > Other GRIN); every sample sits in exactly
    one project, so Data.accessionsFor() never duplicates a sample
  * coded GRIN descriptors are decoded to their labels ("C" -> "Dent")
  * numeric GRIN descriptors are binned into tertiles computed over the v1
    samples that have a value ("Low (<= 214 g)", ...), because the current
    SNPTrait has no numeric-range facet

Usage:
  python3 build_maize_trait_catalog.py --samples data/zmgrin2026_samples.tsv \
     --t5 T5.tsv --t4b T4b.tsv --t4c T4c.tsv --crosswalk crosswalk.tsv --outdir out/SNPTools
"""
import argparse, csv, datetime, json, math, os, re, sys
from collections import Counter, OrderedDict, defaultdict

FAMILY = "zmgrin2026"
FAMILY_LABEL = "MaizeGDB GRIN-linked 2026 (release v1.4)"
PANEL_ORDER = ["NAM", "Ames282", "WiDiv", "Other GRIN"]
PANEL_META = {
    "NAM":        ("NAM founders + B73",          "#7c3aed", "Nested Association Mapping founders (and the B73 reference)."),
    "Ames282":    ("Ames 282 association panel",  "#2563eb", "Goodman/Buckler 282 panel members not already listed under NAM."),
    "WiDiv":      ("Wisconsin Diversity panel",   "#1f8a4c", "WiDiv members not already listed under NAM or Ames282."),
    "Other GRIN": ("Other GRIN accessions",       "#999999", "Grzybowski 2023 lines linked to a GRIN accession outside the three panels (class C; provisional IDs)."),
}
# coded / nominal descriptors surfaced as facets  (trait_code -> (facet key, label))
CODED_FACETS = OrderedDict([
    ("KERNEL-TYPE",  ("kernelType", "Kernel type")),
    ("COB-COLOR",    ("cobColor",   "Cob colour")),
    ("KERNEL-COLOR", ("kernelColor","Kernel colour")),
    ("STEWARTS-WILT",("stewartsWilt","Stewart's wilt score")),
])
# numeric descriptors binned into tertiles
NUMERIC_FACETS = OrderedDict([
    ("KERNEL-WEIGHT-1000",  ("kw1000",      "1000-kernel weight")),
    ("PLANT-HEIGHT",        ("plantHeight", "Plant height")),
    ("GDU-SILK-F_LAT40-49", ("gduSilk",     "GDU to silk (lat 40-49)")),
    ("EAR-LENGTH",          ("earLength",   "Ear length")),
])
# KERNEL-ROW-NUMBER is deliberately NOT faceted: 485 of 841 v1 values are 0 (not a
# possible row count), i.e. the GRIN export mixes 0-as-missing with real counts.
# It stays in data/traits/<family>.traits.json unchanged; see the design doc.


# T4c carries unit "%" for the GDU descriptors, which are growing-degree units
# (degrees F x days, 86/50 method) -- override for display only.
UNIT_OVERRIDE = {"GDU-SILK-F_LAT40-49": "GDU F"}
ORIGIN_PREFIX = re.compile(r"^\s*(Developed|Donated|Collected|Unknown)\s*[\u2013-]\s*", re.I)


def origin_parts(origin):
    """'Developed – Iowa, United States' -> ('Developed', 'Iowa', 'United States')."""
    o = origin or ""
    m = ORIGIN_PREFIX.match(o)
    kind = m.group(1).capitalize() if m else ""
    parts = [x.strip() for x in ORIGIN_PREFIX.sub("", o).split(",") if x.strip()]
    country = parts[-1] if parts else ""
    admin1 = ", ".join(parts[:-1]) if len(parts) > 1 else ""
    return kind, admin1, country


def read_tsv(p):
    with open(p, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def truthy(v):
    return str(v).strip().lower() == "true"


def num(v):
    try:
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def parse_codes(defn):
    """'1=WHITE; 2=RED' -> {'1':'White','2':'Red'}  (labels shortened, title-cased)."""
    out = {}
    for part in re.split(r";\s*(?=[A-Za-z0-9]{1,4}=)", defn or ""):
        m = re.match(r"\s*([A-Za-z0-9]{1,4})=(.+)", part)
        if not m:
            continue
        lab = m.group(2).strip()
        lab = re.sub(r"^\(\d+\)\s*", "", lab)          # "(1) Little or no spread." -> ...
        lab = lab.split(".")[0].strip() or lab
        if lab.isupper():
            lab = lab.capitalize()
        out[m.group(1)] = lab
    # ordinal score scales with long descriptions -> "4 · Abundant spread ..." for every code
    if out and all(k.isdigit() for k in out) and any(len(v) > 24 for v in out.values()):
        out = {k: f"{k} \u00b7 " + (v[:32].rstrip() + "\u2026" if len(v) > 32 else v) for k, v in out.items()}
    return out


def fmt(x):
    return f"{x:.0f}" if abs(x) >= 100 else (f"{x:.1f}" if abs(x) >= 10 else f"{x:.2f}")


def tertile_labeler(values, unit):
    vs = sorted(values)
    if len(vs) < 6:
        return None, None
    q1, q2 = vs[len(vs) // 3], vs[(2 * len(vs)) // 3]
    u = f" {unit}" if unit else ""
    labs = [f"Low (<= {fmt(q1)}{u})", f"Mid ({fmt(q1)}-{fmt(q2)}{u})", f"High (> {fmt(q2)}{u})"]
    def lab(x):
        return labs[0] if x <= q1 else (labs[1] if x <= q2 else labs[2])
    return lab, {"q1": q1, "q2": q2, "labels": labs}


def clean(v):
    s = (v or "").strip()
    return "" if s.lower() in {"", ".", "na", "n/a", "nan", "none", "null", "unknown"} else s


def main():
    ap = argparse.ArgumentParser()
    for k in ("samples", "t5", "t4b", "t4c", "crosswalk"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--family", default=FAMILY)
    ap.add_argument("--bins-from", help="existing <family>.traits.json: reuse its facetBins edges and labels, so "
                    "facet labels stay stable when samples are added (e.g. v1.3 -> v1.4)")
    a = ap.parse_args()

    samples = read_tsv(a.samples)
    t5 = {r["accession_number"]: r for r in read_tsv(a.t5)}
    t4c = {r["trait_code"]: r for r in read_tsv(a.t4c)}
    xw = {r["VCFname"]: r for r in read_tsv(a.crosswalk)}
    t4b = defaultdict(dict)
    for r in read_tsv(a.t4b):
        t4b[r["accession_number"]][r["trait_code"]] = r

    grins = {s["grin_accession"] for s in samples}
    # ---- trait dictionary (only traits observed on the v1 set) --------------------
    used = sorted({c for g in grins for c in t4b.get(g, {})})
    codes = {c: parse_codes(t4c[c]["code_definitions"]) for c in used if c in t4c}
    dictionary = OrderedDict()
    for c in used:
        d = t4c.get(c)
        if not d:
            print(f"!! trait {c} missing from T4c"); continue
        dictionary[c] = {"name": d["trait_name"], "class": d["trait_class"], "category": d["category"],
                         "scale": d["scale_type"], "unit": d["unit"], "url": d["grin_descriptor_url"]}
        if codes.get(c):
            dictionary[c]["codes"] = codes[c]

    # ---- numeric tertiles over v1 samples (one value per GRIN accession) ----------
    binners, bins_meta = {}, {}
    pinned = json.load(open(a.bins_from, encoding="utf-8")).get("facetBins", {}) if a.bins_from else {}
    for c, (key, _) in NUMERIC_FACETS.items():
        vals = [num(t4b[g][c]["mean"]) for g in grins if c in t4b.get(g, {})]
        vals = [v for v in vals if v is not None]
        if c in pinned:
            q1, q2, labs = pinned[c]["q1"], pinned[c]["q2"], pinned[c]["labels"]
            binners[c] = (lambda q1, q2, labs: lambda x: labs[0] if x <= q1 else (labs[1] if x <= q2 else labs[2]))(q1, q2, labs)
            bins_meta[c] = {"q1": q1, "q2": q2, "labels": labs, "n": len(vals),
                            "pinnedFrom": os.path.basename(a.bins_from), "pinnedN": pinned[c].get("n")}
            continue
        lab, meta = tertile_labeler(vals, UNIT_OVERRIDE.get(c, t4c.get(c, {}).get("unit", "")))
        if lab:
            binners[c], bins_meta[c] = lab, dict(meta, n=len(vals))

    # ---- per-sample accession objects --------------------------------------------
    by_panel = defaultdict(list)
    values = OrderedDict()
    fill = Counter()
    for s in samples:
        sid, grin = s["sample_id"], s["grin_accession"]
        p5 = t5.get(grin, {})
        xr = xw.get(s["vcf_name"], {})
        tr = t4b.get(grin, {})
        name = clean(p5.get("accession_name")) or s["accession_name"] or s["vcf_name"]
        acc = OrderedDict([
            ("id", sid), ("run", s["vcf_name"]), ("founder", s["base_id"]),
            ("label", f"{name} ({grin})" if grin else name),
            ("rep", int(s["rep"])), ("reps", int(s["reps"])),
            ("namFounder", s["nam_founder"] or None),
        ])
        meta = OrderedDict()
        meta["strain"] = name                              # SNPTrait display name
        meta["grin"] = grin
        meta["genotypeId"] = s["genotype_id"]
        meta["panel"] = s["primary_panel"]
        meta["inNAM"] = "yes" if truthy(s["pop_NAM"]) else "no"
        meta["inAmes282"] = "yes" if truthy(s["pop_Ames282"]) else "no"
        meta["inWiDiv"] = "yes" if truthy(s["pop_WiDiv"]) else "no"
        # passport (T5) -- class C has no T5 row; fall back to the crosswalk GRIN origin text
        if p5:
            _, _, ctry = origin_parts(clean(p5.get("origin")))
            meta["country"] = ctry or p5.get("country_iso3", "")
            meta["state"] = clean(p5.get("admin1_name"))
            meta["originType"] = clean(p5.get("origin_type"))
            meta["improvement"] = clean(p5.get("improvement_status"))
            meta["pedigree"] = clean(p5.get("pedigree"))
        else:
            kind, adm1, ctry = origin_parts(xr.get("grin_origin", ""))
            meta["country"], meta["state"], meta["originType"] = ctry, adm1, kind
        meta["subpop"] = clean(xr.get("Subpopulation"))
        era = clean(xr.get("ImprovementEra"))
        meta["era"] = {"expvp": "ExPVP"}.get(era.lower(), era)
        meta["species"] = clean(xr.get("Species"))
        for c, (key, _) in CODED_FACETS.items():
            v = tr.get(c, {}).get("mode_value", "")
            if v:
                meta[key] = codes.get(c, {}).get(v, v)
        for c, (key, _) in NUMERIC_FACETS.items():
            x = num(tr.get(c, {}).get("mean"))
            if x is not None and c in binners:
                meta[key] = binners[c](x)
        meta["nTraits"] = str(len(tr))
        for k, v in list(meta.items()):
            if v in ("", None):
                del meta[k]                                 # JS renders absent keys as 'Unknown'
            else:
                fill[k] += 1
        acc.update(meta)
        by_panel[s["primary_panel"]].append(acc)
        # full trait values for the JSON side-file
        tv = OrderedDict()
        for c in sorted(tr):
            r = tr[c]
            if r["scale_type"] == "numeric":
                tv[c] = [int(r["n_obs"] or 0), num(r["mean"]), num(r["sd"]), num(r["min"]), num(r["max"])]
            else:
                tv[c] = [int(r["n_obs"] or 0), r["mode_value"]]
        values[sid] = {"grin": grin, "traits": tv}

    projects = []
    for i, pnl in enumerate([p for p in PANEL_ORDER if by_panel.get(p)]):
        title, color, desc = PANEL_META[pnl]
        rows = sorted(by_panel[pnl], key=lambda r: (r["founder"].lower(), r["rep"]))
        projects.append({
            "id": f"{a.family}_{re.sub(r'[^a-z0-9]', '', pnl.lower())}", "title": title,
            "bioprojects": [], "ncbiUrl": "", "referenceUrl": "https://doi.org/10.1111/tpj.16123",
            "description": desc, "statedTotal": len(rows), "color": color, "count": len(rows),
            "namFounders": sorted({r["namFounder"] for r in rows if r["namFounder"]}),
            "groups": [{"name": None, "accessions": rows}],
        })
    family = {"label": FAMILY_LABEL, "count": sum(p["count"] for p in projects),
              "namFounders": sorted({n for p in projects for n in p["namFounders"]}),
              "projects": projects}

    # ---- SNPTrait schema: exact TRAIT_SCHEMA shape ------------------------------
    facets = [["panel", "Primary panel"], ["inNAM", "In NAM"], ["inAmes282", "In Ames 282"],
              ["inWiDiv", "In WiDiv"], ["subpop", "Subpopulation (Grzybowski 2023)"],
              ["country", "Country of origin"], ["state", "State / province"],
              ["improvement", "Improvement status (GRIN)"], ["originType", "Origin type (GRIN)"],
              ["era", "Improvement era"]]
    facets += [[k, l] for _, (k, l) in CODED_FACETS.items()]
    facets += [[k, l] for c, (k, l) in NUMERIC_FACETS.items() if c in binners]
    schema = {
        "groupBy": "panel", "groupLabel": "germplasm panel",
        "groupOrder": [p for p in PANEL_ORDER],
        "groupColors": {p: PANEL_META[p][1] for p in PANEL_ORDER},
        "groupNames": {p: PANEL_META[p][0] for p in PANEL_ORDER},
        "facets": facets,
        "columns": [["grin", "GRIN accession"], ["strain", "Name"], ["subpop", "Subpopulation"],
                    ["country", "Country"], ["state", "State"], ["kernelType", "Kernel type"]],
        "search": ["id", "strain", "grin", "pedigree", "state", "country", "subpop", "genotypeId"],
    }

    # ---- emit -------------------------------------------------------------------
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    js_path = os.path.join(a.outdir, "js", "zmgrin.catalog.js")
    os.makedirs(os.path.dirname(js_path), exist_ok=True)
    with open(js_path, "w", encoding="utf-8") as fh:
        fh.write("/* AUTO-GENERATED by build_maize_trait_catalog.py -- do not edit by hand.\n"
                 f"   Adds family '{a.family}' (maize GRIN-linked release v1.4, Grzybowski et al. 2023 sites) to\n"
                 "   window.SNP_CATALOG and its SNPTrait schema to window.SNPTRAIT_SCHEMA.\n"
                 "   Load AFTER js/accessions.catalog.js and BEFORE js/data.js.\n"
                 "   Accession id = release VCF sample column (sample_id); run = Grzybowski VCFname;\n"
                 "   founder = base id (replicate columns share it). */\n")
        fh.write("window.SNP_CATALOG = window.SNP_CATALOG || {generated:null, families:{}};\n")
        fh.write(f"window.SNP_CATALOG.families[{json.dumps(a.family)}] = ")
        json.dump(family, fh, separators=(",", ":"), ensure_ascii=False)
        fh.write(";\n")
        fh.write("window.SNPTRAIT_SCHEMA = window.SNPTRAIT_SCHEMA || {};\n")
        fh.write(f"window.SNPTRAIT_SCHEMA[{json.dumps(a.family)}] = ")
        json.dump(schema, fh, separators=(",", ":"), ensure_ascii=False)
        fh.write(";\n")
    tr_path = os.path.join(a.outdir, "data", "traits", f"{a.family}.traits.json")
    os.makedirs(os.path.dirname(tr_path), exist_ok=True)
    with open(tr_path, "w", encoding="utf-8") as fh:
        json.dump({"generated": now, "family": a.family,
                   "valueLayout": {"numeric": ["n_obs", "mean", "sd", "min", "max"],
                                   "ordinal_or_nominal_coded": ["n_obs", "mode_code"],
                                   "text": ["n_obs", "mode_text"]},
                   "facetBins": bins_meta, "dictionary": dictionary, "samples": values},
                  fh, separators=(",", ":"), ensure_ascii=False)

    n = family["count"]
    print(f"[catalog]  family {a.family}: {n} samples in {len(projects)} projects -> {js_path} "
          f"({os.path.getsize(js_path)/1024:.0f} KB)")
    for p in projects:
        print(f"             {p['title']:30s} {p['count']:4d}")
    print(f"[schema]   {len(facets)} facets, {len(schema['columns'])} columns")
    print("[coverage] facet fill (samples with a value):")
    for k, _ in facets:
        print(f"             {k:>12}: {fill.get(k, 0)}/{n}")
    for c, m in bins_meta.items():
        print(f"[bins]     {c}: n={m['n']} q1={m['q1']:.4g} q2={m['q2']:.4g}")
    print(f"[traits]   {len(dictionary)} traits, {sum(len(v['traits']) for v in values.values())} "
          f"sample x trait values -> {tr_path} ({os.path.getsize(tr_path)/1024:.0f} KB)")
    notr = [s for s, v in values.items() if not v["traits"]]
    if notr:
        print(f"!! {len(notr)} samples without any GRIN trait record (e.g. {notr[:5]})")


if __name__ == "__main__":
    main()
