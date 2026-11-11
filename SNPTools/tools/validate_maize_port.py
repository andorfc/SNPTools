#!/usr/bin/env python3
"""
validate_maize_port.py -- check the maize builder outputs against the fusarium
example files (structure) and against an independent recomputation (numbers).

  1. Structure: the maize SNP_CATALOG family, SNPTRAIT_SCHEMA, SNPGEO_REGIONS,
     SNPGEO_COUNTRY_ISO and countries.geo.json have the same container shapes,
     key sets and value types as fusarium js/strains.catalog.js,
     TRAIT_SCHEMA (js/snptrait.js) and js/snpgeo.regions.js.
  2. Integrity: every catalogue id is unique, appears in exactly one project, has a
     region record, and matches a sample column of the test VCFs.
  3. Numbers: per-site, per-country counts produced by the fusarium snpgeo.js
     aggregateGeoData() (+ proposed patch) in smoke_results.json are recomputed
     here directly from the VCF text and SNPGEO_REGIONS; SNPTrait facet counts and
     one compound filter are recomputed from the sample manifest + catalogue.

usage: python3 validate_maize_port.py <repoRoot> <portOut> <smoke_results.json> <vcf.gz>...
"""
import gzip, json, re, sys
from collections import Counter, defaultdict

REPO, OUT, SMOKE, *VCFS = sys.argv[1:]
fails, notes = [], []


def check(cond, msg):
    (notes if cond else fails).append(("PASS " if cond else "FAIL ") + msg)


def js_value(path, anchor):
    src = open(path, encoding="utf-8").read()
    m = re.search(anchor, src)
    if not m:
        raise SystemExit(f"anchor {anchor!r} not in {path}")
    obj, _ = json.JSONDecoder().raw_decode(src[m.end():].lstrip())
    return obj


def trait_schema_js(path):
    """Extract TRAIT_SCHEMA key names/types from the JS object literal (not JSON)."""
    src = open(path, encoding="utf-8").read()
    block = src[src.index("graminearum: {"):src.index("vert: {")]
    block = block[block.index("{") + 1:]
    prev = None
    while prev != block:                      # drop nested {...} / [...] contents
        prev, block = block, re.sub(r"\{[^{}]*\}|\[[^\[\]]*\]", "0", block)
    return set(re.findall(r"(?:^|[,{])\s*(\w+)\s*:", block, re.M))


def tname(v):
    return {dict: "object", list: "array", str: "string", int: "number", float: "number",
            bool: "boolean", type(None): "null"}[type(v)]


# ---------------- load ----------------
f_cat = js_value(f"{REPO}/fusarium/js/strains.catalog.js", r"window\.SNP_CATALOG\s*=")
f_reg = js_value(f"{REPO}/fusarium/js/snpgeo.regions.js", r"window\.SNPGEO_REGIONS\s*=")
f_iso = js_value(f"{REPO}/fusarium/js/snpgeo.regions.js", r"window\.SNPGEO_COUNTRY_ISO\s*=")
m_fam = js_value(f"{OUT}/SNPTools/js/zmgrin.catalog.js", r"window\.SNP_CATALOG\.families\[\"zmgrin2026\"\]\s*=")
m_sch = js_value(f"{OUT}/SNPTools/js/zmgrin.catalog.js", r"window\.SNPTRAIT_SCHEMA\[\"zmgrin2026\"\]\s*=")
m_reg = js_value(f"{OUT}/SNPTools/js/snpgeo.regions.js", r"window\.SNPGEO_REGIONS\s*=")
m_iso = js_value(f"{OUT}/SNPTools/js/snpgeo.regions.js", r"window\.SNPGEO_COUNTRY_ISO\s*=")
m_geo = json.load(open(f"{OUT}/SNPTools/data/geo/countries.geo.json"))
smoke = json.load(open(SMOKE))
manifest = [dict(zip(h, l.rstrip("\n").split("\t"))) for h in [None] for l in []]
with open(f"{OUT}/SNPTools/data/zmgrin2026_samples.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    manifest = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh]

# ---------------- 1. structure ----------------
f_fam = f_cat["families"]["graminearum"]
for lvl, fv, mv in [("family", f_fam, m_fam),
                    ("project", f_fam["projects"][0], m_fam["projects"][0]),
                    ("group", f_fam["projects"][0]["groups"][0], m_fam["projects"][0]["groups"][0])]:
    fk = {k: tname(v) for k, v in fv.items()}
    mk = {k: tname(v) for k, v in mv.items()}
    check(fk.keys() == mk.keys(), f"catalog {lvl} keys identical to fusarium ({sorted(fk)})")
    diff = {k for k in fk if k in mk and fk[k] != mk[k] and "null" not in (fk[k], mk[k])}
    check(not diff, f"catalog {lvl} value types match fusarium" + (f" (diff: {diff})" if diff else ""))
BASE = ["id", "run", "founder", "label", "rep", "reps", "namFounder"]
f_acc = [a for p in f_fam["projects"] for g in p["groups"] for a in g["accessions"]]
m_acc = [a for p in m_fam["projects"] for g in p["groups"] for a in g["accessions"]]
check(all(list(a)[:7] == BASE for a in f_acc), "fusarium accessions start with the 7 base keys (reference)")
check(all(list(a)[:7] == BASE for a in m_acc), "maize accessions start with the same 7 base keys in the same order")
check(all(isinstance(a["rep"], int) and isinstance(a["reps"], int) for a in m_acc), "rep/reps are integers")
check(all(isinstance(v, str) for a in m_acc for k, v in a.items() if k not in ("rep", "reps", "namFounder")),
      "all maize metadata values are scalar strings (SNPTrait facet matching is string equality)")
fsk = trait_schema_js(f"{REPO}/fusarium/js/snptrait.js")
check(set(m_sch) == fsk, f"SNPTRAIT_SCHEMA keys == fusarium TRAIT_SCHEMA keys ({sorted(fsk)})")
check(all(isinstance(f, list) and len(f) == 2 for f in m_sch["facets"] + m_sch["columns"]), "facets/columns are [key,label] pairs")
acc_keys = set().union(*[a.keys() for a in m_acc])
missing = [k for k in [m_sch["groupBy"]] + [f[0] for f in m_sch["facets"]] + [c[0] for c in m_sch["columns"]] + m_sch["search"]
           if k not in acc_keys]
check(not missing, "every schema key exists on the maize accessions" + (f" (missing {missing})" if missing else ""))
f_rkeys = set().union(*[r.keys() for r in f_reg.values()])
m_rkeys = set().union(*[r.keys() for r in m_reg.values()])
check({"iso3", "country", "state"} <= m_rkeys, "region records carry the keys snpgeo.js reads (iso3, country, state)")
notes.append(f"INFO region keys shared with fusarium: {sorted(m_rkeys & f_rkeys)}; maize-only: {sorted(m_rkeys - f_rkeys)}")
for k in m_rkeys & f_rkeys:
    ft = {tname(r[k]) for r in f_reg.values() if k in r}
    mt = {tname(r[k]) for r in m_reg.values() if k in r}
    check(mt <= ft, f"region key '{k}' type {mt} consistent with fusarium {ft}")
check(all(isinstance(k, str) and re.fullmatch(r"[A-Z]{3}", v) for k, v in m_iso.items()), "SNPGEO_COUNTRY_ISO is name -> ISO3")
check(m_geo["type"] == "FeatureCollection" and all(set(f["properties"]) == {"iso3", "name"} for f in m_geo["features"]),
      "countries.geo.json: FeatureCollection with properties {iso3,name}")
polys = {f["properties"]["iso3"] for f in m_geo["features"]}
check(all(r["iso3"] in polys for r in m_reg.values() if "iso3" in r), "every region iso3 has a polygon")
check(all(("lat" in r) == ("lon" in r) for r in m_reg.values()), "lat/lon emitted only as pairs")
check(all(isinstance(r.get("lat", 0.0), float) and -90 <= r.get("lat", 0) <= 90 and -180 <= r.get("lon", 0) <= 180
          for r in m_reg.values()), "lat/lon are in-range floats")

# ---------------- 2. integrity ----------------
ids = [a["id"] for a in m_acc]
check(len(ids) == len(set(ids)) == m_fam["count"] == len(manifest), f"{len(ids)} unique ids, == family count == manifest rows")
check(set(ids) == set(m_reg), "catalogue ids == SNPGEO_REGIONS keys")
check(sum(p["count"] for p in m_fam["projects"]) == len(ids), "project counts partition the family (no sample in two projects)")
vcf_samples = []
for v in VCFS:
    with gzip.open(v, "rt") as fh:
        for l in fh:
            if l.startswith("#CHROM"):
                vcf_samples.append(l.rstrip("\n").split("\t")[9:]); break
for v, s in zip(VCFS, vcf_samples):
    check(set(s) == set(ids), f"{v.split('/')[-1]}: {len(s)} VCF sample columns == catalogue ids")
nf = sorted({a["namFounder"] for a in m_acc if a["namFounder"]})
check(m_fam["namFounders"] == nf, f"family.namFounders lists the {len(nf)} tagged founders")

# ---------------- 3. numbers ----------------
prop = smoke["maize_proposed"]
man = {r["sample_id"]: r for r in manifest}
acc_by = {a["id"]: a for a in m_acc}
fc = prop["trait"]["facetCounts"]
check(fc["panel"] == dict(Counter(r["primary_panel"] for r in manifest)), "SNPTrait panel facet counts == manifest primary_panel counts")
for k in ("inNAM", "inAmes282", "inWiDiv"):
    col = "pop_" + k[2:]
    exp = Counter("yes" if man[i][col] == "True" else "no" for i in ids)
    check(fc[k] == dict(exp), f"SNPTrait {k} facet counts == manifest {col}")
for k in ("subpop", "country", "kernelType", "kw1000"):
    exp = Counter(acc_by[i].get(k, "Unknown") for i in ids)
    check(fc[k] == dict(exp), f"SNPTrait {k} facet counts == catalogue values (Unknown for absent)")
exp_f = sorted(i for i in ids if acc_by[i].get("subpop") == "SS" and acc_by[i].get("inAmes282") == "yes"
               and acc_by[i].get("kernelType") == "Dent")
check(prop["trait"]["filterSS_Ames_Dent"] == exp_f, f"compound filter SS & Ames282 & Dent: {len(exp_f)} samples, same ids")

CODE = {"0|0": 0, "0/0": 0, "0|1": 1, "1|0": 1, "0/1": 1, "1/0": 1, "1|1": 2, "1/1": 2}
for v, g in zip(VCFS, prop["geo"]):
    with gzip.open(v, "rt") as fh:
        lines = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("##")]
    cols, body = lines[0][9:], lines[1:]
    iso = {s: m_reg[s].get("iso3", "???") for s in cols}
    totals = Counter(r.get("iso3", "???") for r in m_reg.values())
    check(len(body) == g["nRows"], f"{v.split('/')[-1]}: {len(body)} sites parsed by main parseVcf")
    bad = 0
    het_total = 0
    for t, site in zip(body, g["sites"]):
        exp = defaultdict(lambda: Counter())
        for s, gt in zip(cols, t[9:]):
            c = CODE.get(gt.split(":")[0], 3)
            e = exp[iso[s]]
            if c == 3:
                e["missing"] += 1; continue
            e["called"] += 1; e["calledAlleles"] += 2; e["altAlleles"] += c
            if c: e["count"] += 1
            else: e["ref"] += 1
            if c == 1: e["het"] += 1; het_total += 1
        for k3, e in exp.items():
            got = site["c"].get(k3, {})
            for f in ("called", "count", "ref", "het", "altAlleles", "calledAlleles"):
                if got.get(f, 0) != e.get(f, 0):
                    bad += 1
            if got.get("total") != totals[k3]:
                bad += 1
        if int(t[1]) != site["pos"]:
            bad += 1
    check(bad == 0, f"{v.split('/')[-1]}: per-site x per-country counts match an independent recomputation "
                    f"({len(body)} sites x {len(totals)} country buckets; {het_total} het calls)")
un = smoke["maize_unpatched"]
check(all(sum(x["called"] for x in s["c"].values()) == 0 for gg in un["geo"] for s in gg["sites"]),
      "REPRODUCED: unpatched snpgeo.js classifies every Int8Array genotype as missing")
check(un["trait"]["schemaGroupBy"] == "country" and not un["trait"]["hasMeta"],
      "REPRODUCED: unpatched snptrait.js falls through to the F. verticillioides schema and sees no metadata")
check(smoke["control"]["vertGeoTotalAll"] == 511 and smoke["control"]["vertGeoCalledAll"] == 0,
      "REPRODUCED (fusarium branch): verticillioides SNPGeo totals come from the 511 graminearum region records")
check(str(prop["scoreModelsCall"]).startswith("Data.scoreModels is not a function"),
      "REPRODUCED: main Data has no scoreModels (geoRenderTable would throw)")

print("\n".join(notes + fails))
print(f"\n{len([n for n in notes if n.startswith('PASS')])} passed, {len(fails)} failed")
sys.exit(1 if fails else 0)
