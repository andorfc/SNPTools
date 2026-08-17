#!/usr/bin/env python3
"""
build_strains_catalog.py — compile the Fusarium strain metadata + HDF5 accession
columns into the data files the SNPTools front-end consumes:

  js/strains.catalog.js   window.SNP_CATALOG  (families -> projects -> groups -> accessions)
                          window.FUSARIUM_STRAINS  (full graminearum metadata for SNPTrait)

Sources:
  strains.js                 graminearum strain metadata (id, population, species, host,
                             chemotype, country)  -> 512 isolates, grouped by population
  <family>_acc.json          accession column names read from each verticillioides HDF5
                             (no metadata available yet -> flagged in MISSING_DATA.md)

Hierarchy:  reference family -> population/collection project -> group -> accession
"""
import re, json, datetime, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
STRAINS_JS = sys.argv[1] if len(sys.argv) > 1 else "strains.js"
VERT7600_ACC = sys.argv[2] if len(sys.argv) > 2 else "/tmp/vert7600_acc.json"
VERTMRC_ACC  = sys.argv[3] if len(sys.argv) > 3 else "/tmp/vertMRC826_acc.json"
OUT_JS = os.path.join(HERE, "js", "strains.catalog.js")

# population display + colours (from the original Strain Selector, selector.js)
POP_ORDER  = ["NA1", "NA2", "NA3", "Admixture", "Outgroups", "Unknown"]
POP_NAMES  = {"NA1":"North American 1 (NA1)", "NA2":"North American 2 (NA2)",
              "NA3":"North American 3 (NA3)", "Admixture":"Admixture",
              "Outgroups":"Outgroups", "Unknown":"Unknown"}
POP_COLORS = {"NA1":"#56B4E9","NA2":"#e0c200","NA3":"#009E73",
              "Admixture":"#E69F00","Outgroups":"#CC79A7","Unknown":"#999999"}

def parse_strains_js(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    m = re.search(r'const\s+data\s*=\s*(\[\s*[\s\S]*?\]);', txt)
    if not m:
        raise RuntimeError("could not find `const data = [...]` in " + path)
    body = m.group(1)
    body = re.sub(r'(\b\w+)\s*:', r'"\1":', body)      # quote bare keys (already quoted stay valid)
    body = re.sub(r'""(\w+)""', r'"\1"', body)          # undo double-quoting of already-quoted keys
    body = re.sub(r',\s*([}\]])', r'\1', body)          # drop trailing commas
    return json.loads(body)

from collections import Counter

# VCF spec columns / pandas junk that are NOT real accessions (must never appear
# as isolates in the picker).
NONACC = {"CHROM","POS","REF","ALT","QUAL","INFO","FILTER","FORMAT","ID"}
def is_accession(x):
    x = str(x)
    return x not in NONACC and not x.startswith("Unnamed")

strains = parse_strains_js(STRAINS_JS)
for s in strains:
    s["population"] = (s.get("population") or "Unknown").strip()
    s["species"]    = re.sub(r"\s+", " ", (s.get("species") or "unknown")).strip()
    s["host"]       = (s.get("host") or "unknown").strip()
    s["chemotype"]  = str(s.get("chemotype") or "unknown").strip()
    s["country"]    = (s.get("country") or "unknown").strip()

# Replicate disambiguation: the HDF5 stores duplicate sample columns with a
# pandas-style ".N" suffix (ET-2022-49-1, ET-2022-49-1.1). The strain metadata
# lists both under the same base id, so mirror that suffixing here IN METADATA
# ORDER so each replicate maps to its real HDF5 column and both are selectable.
# `baseId` keeps the clean display name; `rep`/`reps` drive the "· rN" chip badge.
_id_total = Counter(str(s["id"]) for s in strains)
_id_seen  = Counter()
for s in strains:
    base = str(s["id"]); k = _id_seen[base]; _id_seen[base] += 1
    s["baseId"] = base
    s["rep"]    = k + 1
    s["reps"]   = _id_total[base]
    s["id"]     = base if k == 0 else f"{base}.{k}"

# ---- graminearum family: one project per population ----
def acc_obj(sid, meta=None):
    # base display name: from metadata when present (keeps "· rN" replicates on
    # one clean name), else the raw column id (verticillioides).
    base = str(meta["baseId"]) if (meta and meta.get("baseId")) else str(sid)
    o = {"id": str(meta["id"]) if meta else str(sid),
         "run": base, "founder": base, "label": base,
         "rep": meta.get("rep", 1) if meta else 1,
         "reps": meta.get("reps", 1) if meta else 1,
         "namFounder": None}
    if meta:
        o.update({"population": meta["population"], "species": meta["species"],
                  "host": meta["host"], "chemotype": meta["chemotype"],
                  "country": meta["country"]})
    return o

gram_projects = []
by_pop = {p: [] for p in POP_ORDER}
for s in strains:
    by_pop.setdefault(s["population"], []).append(s)
for i, pop in enumerate([p for p in POP_ORDER if by_pop.get(p)] +
                        [p for p in by_pop if p not in POP_ORDER]):
    rows = sorted(by_pop[pop], key=lambda r: (str(r["id"])))
    gram_projects.append({
        "id": "gram_" + pop.lower(),
        "title": POP_NAMES.get(pop, pop),
        "bioprojects": [], "ncbiUrl": "", "referenceUrl": "",
        "description": f"{pop} population-structure cluster",
        "statedTotal": len(rows), "color": POP_COLORS.get(pop, "#888"),
        "count": len(rows), "namFounders": [],
        "groups": [{"name": None, "accessions": [acc_obj(r["id"], r) for r in rows]}],
    })

def vert_family(acc_path, label):
    # NOTE: this emits verticillioides accessions WITHOUT the SNPTrait metadata
    # (strain/region/substrate/geo/host from Fvert_MetaData_2024-05-20.xlsx) — that
    # enrichment is a SEPARATE merge step. Do NOT overwrite a metadata-enriched
    # js/strains.catalog.js with a raw rerun of this script or SNPTrait vert facets
    # will go blank. Re-run the vert-metadata merge after this, or patch in place.
    acc = json.load(open(acc_path)) if os.path.exists(acc_path) else []
    acc = [a for a in acc if is_accession(a)]          # drop VCF spec cols / Unnamed junk
    acc = sorted(acc, key=lambda x: (len(str(x)), str(x)))
    return {
        "label": label,
        "count": len(acc),
        "namFounders": [],
        "projects": [{
            "id": "all", "title": "All isolates", "bioprojects": [],
            "ncbiUrl": "", "referenceUrl": "",
            "description": "Accession columns read from the HDF5 store (metadata pending).",
            "statedTotal": len(acc), "color": "#2563eb", "count": len(acc),
            "namFounders": [],
            "groups": [{"name": None, "accessions": [acc_obj(a) for a in acc]}],
        }],
    }

catalog = {
    "generated": datetime.datetime.utcnow().isoformat() + "Z",
    "families": {
        "graminearum": {"label": "F. graminearum 2025",
                        "count": sum(p["count"] for p in gram_projects),
                        "namFounders": [], "projects": gram_projects},
        "vert7600":    vert_family(VERT7600_ACC, "F. verticillioides 7600"),
        "vertMRC826":  vert_family(VERTMRC_ACC, "F. verticillioides MRC826"),
    },
}

compact = json.dumps(catalog, separators=(",", ":"))
strains_compact = json.dumps(strains, separators=(",", ":"))
os.makedirs(os.path.dirname(OUT_JS), exist_ok=True)
with open(OUT_JS, "w", encoding="utf-8") as f:
    f.write("/* AUTO-GENERATED by build_strains_catalog.py — do not edit by hand.\n")
    f.write("   Load BEFORE data.js.  window.SNP_CATALOG feeds the SNPVersity picker;\n")
    f.write("   window.FUSARIUM_STRAINS feeds the SNPTrait metadata Strain Selector. */\n")
    f.write("window.SNP_CATALOG = " + compact + ";\n")
    f.write("window.FUSARIUM_STRAINS = " + strains_compact + ";\n")

print("graminearum isolates:", sum(p["count"] for p in gram_projects))
for p in gram_projects:
    print(f"  {p['title']:28s} {p['count']:4d}")
print("vert7600 accessions:", catalog["families"]["vert7600"]["count"])
print("vertMRC826 accessions:", catalog["families"]["vertMRC826"]["count"])
print("wrote", OUT_JS, os.path.getsize(OUT_JS), "bytes")
