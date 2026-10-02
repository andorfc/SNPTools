#!/usr/bin/env python3
"""
build_maize_samples.py -- derive the SNPTools sample manifest for the maize
GRIN-linked v1 release from the project tables.

One row per VCF sample column of the v1 set (crosswalk classes A, B and C of
grz2023_grin_crosswalk.tsv). This manifest is the single join key used by
build_maize_trait_catalog.py and build_maize_snpgeo_data.py:

    sample_id        SNPTools / release VCF column name (e.g. ZmG_B73)
    base_id          sample_id without the replicate suffix
    rep, reps        replicate index / count for base_id (fusarium convention)
    vcf_name         Grzybowski et al. 2023 VCF column (source genotype)
    genotype_id      ZmG:<key> from T0 (blank for class C)
    grin_accession   GRIN accession_number (e.g. "PI 550473")
    grin_dbid        GRIN germplasm dbid (blank when not in T0)
    dataset_class    crosswalk class (A/B/C)
    pop_NAM, pop_Ames282, pop_WiDiv   True/False
    primary_panel    NAM > Ames282 > WiDiv > Other GRIN (partition for SNPTrait)
    nam_founder      founder name for NAM founders (drives SNPVersity default selection)
    id_provisional   True when sample_id was minted here rather than from genotype_id

ID policy (provisional -- see the design doc, section 3):
  * sample_id = genotype_id with ':' -> '_'   (ZmG:B73 -> ZmG_B73)
  * class C rows have no genotype_id; they get ZmG_<NAMEKEY> where NAMEKEY is the
    upper-cased Uniform_Name with non-alphanumerics removed, flagged provisional.
  * When several VCF columns map to one genotype, the column whose normalised
    VCFname equals the genotype key keeps base_id; the others get base_id.1,
    base_id.2 ... in crosswalk row order (the same ".N" replicate convention as
    fusarium build_strains_catalog.py). Pass --release-manifest to override
    with the IDs actually written into the release VCFs.

Usage:
  python3 build_maize_samples.py --t0 T0_population_members.tsv \
      --crosswalk grz2023_grin_crosswalk.tsv --out data/zmgrin2026_samples.tsv
"""
import argparse, csv, re, sys
from collections import Counter, defaultdict

V1_CLASSES = ("A_", "B_", "C_")
PANEL_ORDER = ["NAM", "Ames282", "WiDiv", "Other GRIN"]


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def namekey(s):
    return re.sub(r"[^A-Z0-9]", "", (s or "").upper())


def truthy(v):
    return str(v).strip().lower() == "true"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t0", required=True)
    ap.add_argument("--crosswalk", required=True)
    ap.add_argument("--release-manifest", help="optional TSV with vcf_name, sample_id columns "
                    "(IDs as written into the release VCF); overrides the derived IDs")
    ap.add_argument("--extra-samples", help="optional TSV with vcf_name, genotype_id columns for lines "
                    "called outside the Grzybowski set (e.g. the 8 unmatched lines); GRIN accession and "
                    "panel flags are taken from T0 (primary accession), dataset_class = E_called_separately")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    t0 = read_tsv(a.t0)
    xw = read_tsv(a.crosswalk)
    t0_by_geno = defaultdict(list)
    for r in t0:
        t0_by_geno[r["genotype_id"]].append(r)
    primary_acc = {g: next((r for r in rs if truthy(r["is_primary_accession_of_genotype"])), rs[0])
                   for g, rs in t0_by_geno.items()}

    rows = [r for r in xw if r["dataset_class"].startswith(V1_CLASSES)]
    rows.sort(key=lambda r: int(r["grz_row"]))
    if a.extra_samples:
        seen = {r["VCFname"] for r in rows}
        for k, e in enumerate(read_tsv(a.extra_samples)):
            g = e["genotype_id"].strip()
            pa = primary_acc.get(g)
            if not pa:
                sys.exit(f"!! extra sample {e['vcf_name']}: genotype_id {g!r} not in T0")
            if e["vcf_name"] in seen:
                sys.exit(f"!! extra sample {e['vcf_name']} already in the crosswalk v1 set")
            rows.append({"grz_row": str(100000 + k), "VCFname": e["vcf_name"], "Uniform_Name": pa["accession_name"],
                         "genotype_id": g, "grin_accession": pa["accession_number"],
                         "dataset_class": "E_called_separately", "match_tier": "t0_primary",
                         "pop_NAM": pa["pop_NAM"], "pop_Ames282": pa["pop_Ames282"], "pop_WiDiv": pa["pop_WiDiv"]})
    if not rows:
        sys.exit("!! no v1 rows (classes A/B/C) in the crosswalk")

    # ---- base ids ------------------------------------------------------------
    for r in rows:
        g = r["genotype_id"].strip()
        if g:
            if not g.startswith("ZmG:"):
                sys.exit(f"!! unexpected genotype_id {g!r}")
            r["_base"], r["_prov"] = g.replace(":", "_", 1), False
        else:
            r["_base"], r["_prov"] = "ZmG_" + namekey(r["Uniform_Name"] or r["VCFname"]), True

    # minted class-C ids must not collide with a real genotype id
    real = {r["_base"] for r in rows if not r["_prov"]} | {g.replace(":", "_", 1) for g in t0_by_geno}
    for r in rows:
        if r["_prov"] and r["_base"] in real:
            r["_base"] += "_C"                     # disambiguate, still provisional

    # ---- replicate suffixes --------------------------------------------------
    by_base = defaultdict(list)
    for r in rows:
        by_base[r["_base"]].append(r)
    for base, grp in by_base.items():
        key = base.split("_", 1)[1]
        grp.sort(key=lambda r: (namekey(r["VCFname"]) != key, int(r["grz_row"])))
        for k, r in enumerate(grp):
            r["_rep"], r["_reps"] = k + 1, len(grp)
            r["_sid"] = base if k == 0 else f"{base}.{k}"

    if a.release_manifest:
        rel = {m["vcf_name"]: m["sample_id"] for m in read_tsv(a.release_manifest)}
        miss = [r["VCFname"] for r in rows if r["VCFname"] not in rel]
        if miss:
            sys.exit(f"!! {len(miss)} v1 VCF columns absent from the release manifest, e.g. {miss[:5]}")
        for r in rows:
            r["_sid"] = rel[r["VCFname"]]

    sids = Counter(r["_sid"] for r in rows)
    dup = [s for s, n in sids.items() if n > 1]
    if dup:
        sys.exit(f"!! duplicate sample_id(s): {dup[:10]}")

    # ---- emit -----------------------------------------------------------------
    cols = ["sample_id", "base_id", "rep", "reps", "vcf_name", "genotype_id", "grin_accession",
            "grin_dbid", "accession_name", "dataset_class", "match_tier", "pop_NAM", "pop_Ames282",
            "pop_WiDiv", "primary_panel", "nam_founder", "id_provisional"]
    out = []
    for r in rows:
        g = r["genotype_id"].strip()
        pa = primary_acc.get(g, {})
        pops = {p: truthy(r.get("pop_" + p, "")) for p in ("NAM", "Ames282", "WiDiv")}
        panel = next((p for p in ("NAM", "Ames282", "WiDiv") if pops[p]), "Other GRIN")
        role = pa.get("NAM_role", "")
        nam_founder = (r["Uniform_Name"] or pa.get("accession_name", "")) if (pops["NAM"] and role.startswith(("founder", "reference"))) else ""
        grin = r["grin_accession"].strip()
        if pa and grin and pa.get("accession_number") != grin:
            print(f"!! {r['VCFname']}: crosswalk GRIN {grin} != T0 primary {pa.get('accession_number')}")
        out.append({
            "sample_id": r["_sid"], "base_id": r["_base"], "rep": r["_rep"], "reps": r["_reps"],
            "vcf_name": r["VCFname"], "genotype_id": g, "grin_accession": grin,
            "grin_dbid": pa.get("grin_germplasm_dbid", ""), "accession_name": pa.get("accession_name") or r["Uniform_Name"],
            "dataset_class": r["dataset_class"], "match_tier": r["match_tier"],
            "pop_NAM": pops["NAM"], "pop_Ames282": pops["Ames282"], "pop_WiDiv": pops["WiDiv"],
            "primary_panel": panel, "nam_founder": nam_founder, "id_provisional": r["_prov"],
        })
    with open(a.out, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)

    reps = [o for o in out if int(o["reps"]) > 1]
    print(f"[samples]  {len(out)} v1 sample columns -> {a.out}")
    print(f"[samples]  classes: {dict(Counter(o['dataset_class'][:1] for o in out))}")
    print(f"[samples]  panels : {dict(Counter(o['primary_panel'] for o in out))}")
    print(f"[samples]  {len({o['base_id'] for o in out})} distinct base ids; "
          f"{len({o['base_id'] for o in reps})} base ids carry {len(reps)} replicate columns")
    print(f"[samples]  {sum(o['id_provisional'] for o in out)} provisional (class C) ids; "
          f"{sum(1 for o in out if o['nam_founder'])} NAM founder/reference columns")


if __name__ == "__main__":
    main()
