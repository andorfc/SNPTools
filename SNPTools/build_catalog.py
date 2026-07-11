#!/usr/bin/env python3
"""Compile the legacy hard-coded accession checkboxes in index2026_v20.html into
clean, maintainable data files:

  data/accessions.tsv   one row per accession  (the master flat file)
  data/projects.tsv     one row per bioproject section (display metadata)
  js/accessions.catalog.json   compiled hierarchy for the UI (fetchable)
  js/accessions.catalog.js     same, as window.SNP_CATALOG (drop-in <script>)

Hierarchy:  dataset family -> bioproject section -> group -> accession
"""
import re, json, datetime, os
from collections import defaultdict, OrderedDict
from bs4 import BeautifulSoup

SRC = "/mnt/user-data/uploads/index2026_v20.html"
OUTJS = "/home/claude/out/js"
OUTDATA = "/home/claude/out/data"
os.makedirs(OUTJS, exist_ok=True)
os.makedirs(OUTDATA, exist_ok=True)

CLASS_TO_FAMILY = {
    "genotypes26":  "mgdb2026",
    "genotypes":    "mgdb2024",
    "genotypes_S":  "schnable2023",
    "genotypesNAM": "nam2021",
}
FAMILY_LABEL = {
    "mgdb2026":     "MaizeGDB 2026",
    "mgdb2024":     "MaizeGDB 2024",
    "schnable2023": "Schnable 2023",
    "nam2021":      "NAM 2021",
}
# a palette cycled per-project so the table's project bar varies by bioproject
PALETTE = ["#2563eb","#0891b2","#7c3aed","#db2777","#ca8a04","#16a34a",
           "#dc2626","#0d9488","#9333ea","#ea580c","#4f46e5","#65a30d"]

PRJ_RE = re.compile(r"PRJ[A-Z]{2}\d+")
RUN_RE = re.compile(r"^(?:SRR|ERR|DRR|CRR|CRX|ERS|SRS|DRS|ERX|SRX)\d+$", re.I)
REP_RE = re.compile(r"(?:v5)?rep(\d+)$", re.I)

def clean(s): return re.sub(r"\s+", " ", (s or "")).strip()
def cls_of(inp):
    c = inp.get("class") or []
    return " ".join(c) if isinstance(c, list) else c

def parse_id(value, founder_attr):
    parts = value.split("_")
    if len(parts) > 1 and RUN_RE.match(parts[-1]):
        run, label = parts[-1], "_".join(parts[:-1])
    else:
        run, label = value, value
    rep = 1
    rm = REP_RE.search(label)
    if rm:
        rep = int(rm.group(1)); label = label[:rm.start()].rstrip("_-. ")
    founder = founder_attr or label or value
    return run, founder, rep

soup = BeautifulSoup(open(SRC, encoding="utf-8", errors="replace").read(), "lxml")

# map each panel element -> its accordion button (title/metadata)
panel_meta = {}
for btn in soup.find_all("button", class_="accordion"):
    panel = btn.find_next_sibling(lambda t: t.name == "div" and "panel" in (t.get("class") or []))
    if panel is None:
        continue
    title = clean(btn.get_text(" "))
    a = btn.find("a", href=True)
    intro = clean(panel.get_text(" "))
    tm = re.search(r"(\d[\d,]*)\s+total accessions", intro)
    ref = None
    for aa in panel.find_all("a", href=True):
        ref = aa["href"]; break
    panel_meta[id(panel)] = {
        "panel": panel, "title": title,
        "bioprojects": PRJ_RE.findall(title),
        "ncbi_url": a["href"] if a else "",
        "ref_url": ref or "",
        "description": intro[:160],
        "stated_total": int(tm.group(1).replace(",", "")) if tm else None,
    }

def nearest_panel(el):
    for p in el.parents:
        if p.name == "div" and "panel" in (p.get("class") or []):
            return p
    return None

def group_for(inp, panel):
    """the most recent 'Group: X' text before this input, within its accordion
    section (scan back until the section's accordion button, so groups declared
    in an outer panel still apply to inputs in nested inner panels)."""
    for prev in inp.previous_elements:
        nm = getattr(prev, "name", None)
        if nm == "button" and "accordion" in (prev.get("class") or []):
            break
        if nm is None:
            m = re.search(r"Group:\s*(.+)", clean(str(prev)))
            if m:
                return clean(re.sub(r"</?b>", "", m.group(1)))
    return None

# ---- single pass over checkboxes; nearest-panel assignment; dedupe by (family,id) ----
# projects keyed by panel id, preserving document order per family
proj_order = defaultdict(list)          # family -> [panel_id ...]
proj_by_pid = {}                        # panel_id -> project dict
seen = set()                            # (family, id)
rows = []                               # flat accession rows

for inp in soup.find_all("input", attrs={"type": "checkbox"}):
    fam = CLASS_TO_FAMILY.get(cls_of(inp).strip())
    val = inp.get("value", "")
    if not fam or val in ("", "skip"):
        continue
    if (fam, val) in seen:
        continue
    seen.add((fam, val))
    panel = nearest_panel(inp)
    pid = id(panel) if panel is not None else None
    meta = panel_meta.get(pid, {"title": FAMILY_LABEL.get(fam, fam), "bioprojects": [],
                                "ncbi_url": "", "ref_url": "", "description": "", "stated_total": None})
    if pid not in proj_by_pid:
        proj_order[fam].append(pid)
        proj_by_pid[pid] = {**meta, "family": fam, "accs": []}
    grp = group_for(inp, panel) if panel is not None else None
    # label = visible text right after the checkbox
    sib = inp.next_sibling
    label = clean(sib if isinstance(sib, str) else "") or val
    founder_attr = inp.get("data-founder")
    run, founder, rep = parse_id(val, founder_attr)
    rec = {"id": val, "label": label, "group": grp, "founder": founder,
           "rep": rep, "run": run, "nam_founder": founder_attr}
    proj_by_pid[pid]["accs"].append(rec)
    rows.append({"family": fam, "panel_id": pid, **rec,
                 "project_title": meta["title"], "bioprojects": ";".join(meta["bioprojects"])})

# assign stable project ids + colors, compute reps-per-founder per family
family_founder_counts = defaultdict(lambda: defaultdict(int))
for r in rows:
    family_founder_counts[r["family"]][r["founder"]] += 1

catalog = {"generated": datetime.datetime.utcnow().isoformat() + "Z", "families": {}}
pid_to_projectid = {}
for fam, pids in proj_order.items():
    projects = []
    for i, pid in enumerate(pids):
        p = proj_by_pid[pid]
        projectid = f"{fam}_p{i+1:02d}"
        pid_to_projectid[pid] = projectid
        color = PALETTE[i % len(PALETTE)]
        # bucket accessions into groups (preserve first-seen order; None -> ungrouped)
        groups = OrderedDict()
        nam_founders = []
        for a in p["accs"]:
            a["reps"] = family_founder_counts[fam][a["founder"]]
            if a["nam_founder"] and a["nam_founder"] not in nam_founders:
                nam_founders.append(a["nam_founder"])
            gname = a["group"] or ""
            groups.setdefault(gname, []).append({
                "id": a["id"], "label": a["label"], "founder": a["founder"],
                "rep": a["rep"], "run": a["run"], "reps": a["reps"],
                "namFounder": a["nam_founder"],
            })
        projects.append({
            "id": projectid, "title": p["title"], "bioprojects": p["bioprojects"],
            "ncbiUrl": p["ncbi_url"], "referenceUrl": p["ref_url"],
            "description": p["description"], "statedTotal": p["stated_total"],
            "color": color, "count": len(p["accs"]),
            "namFounders": sorted(nam_founders),
            "groups": [{"name": (g or None), "accessions": v} for g, v in groups.items()],
        })
    fam_nam = sorted({a for pr in projects for a in pr["namFounders"]})
    catalog["families"][fam] = {"label": FAMILY_LABEL.get(fam, fam),
                                "count": sum(pr["count"] for pr in projects),
                                "namFounders": fam_nam, "projects": projects}

# ---- write master TSVs (editable source of truth) ----
with open(f"{OUTDATA}/accessions.tsv", "w", encoding="utf-8") as f:
    f.write("dataset\tproject_id\tproject_title\tgroup\tid\tlabel\tfounder\trep\trun\tnam_founder\tbioprojects\n")
    for r in rows:
        f.write("\t".join(str(x) for x in [
            r["family"], pid_to_projectid[r["panel_id"]], r["project_title"],
            r["group"] or "", r["id"], r["label"], r["founder"], r["rep"], r["run"],
            r["nam_founder"] or "", r["bioprojects"],
        ]) + "\n")

with open(f"{OUTDATA}/projects.tsv", "w", encoding="utf-8") as f:
    f.write("dataset\tproject_id\ttitle\tbioprojects\tncbi_url\treference_url\tstated_total\tn_accessions\tgroups\n")
    for fam in proj_order:
        for pr in catalog["families"][fam]["projects"]:
            gnames = ";".join(g["name"] for g in pr["groups"] if g["name"])
            f.write("\t".join(str(x) for x in [
                fam, pr["id"], pr["title"], ";".join(pr["bioprojects"]),
                pr["ncbiUrl"], pr["referenceUrl"], pr["statedTotal"] if pr["statedTotal"] is not None else "",
                pr["count"], gnames,
            ]) + "\n")

# ---- write compiled catalog (json + drop-in js) ----
compact = json.dumps(catalog, separators=(",", ":"))
with open(f"{OUTJS}/accessions.catalog.json", "w", encoding="utf-8") as f:
    f.write(compact)
with open(f"{OUTJS}/accessions.catalog.js", "w", encoding="utf-8") as f:
    f.write("/* AUTO-GENERATED from data/accessions.tsv + data/projects.tsv by build_catalog.py.\n")
    f.write("   Load this BEFORE data.js, or fetch accessions.catalog.json instead. */\n")
    f.write("window.SNP_CATALOG = " + compact + ";\n")

# ---- report ----
print("family        projects  accessions  (grep truth)")
truth = {"mgdb2026":2710,"mgdb2024":1498,"schnable2023":1515,"nam2021":27}
for fam in proj_order:
    c = catalog["families"][fam]
    ok = "OK" if c["count"] == truth.get(fam) else f"!! expected {truth.get(fam)}"
    print(f"  {fam:13s} {len(c['projects']):5d}   {c['count']:8d}   {ok}")
print("\nNAM founders (mgdb2026):", len(catalog["families"]["mgdb2026"]["namFounders"]),
      "->", ", ".join(catalog["families"]["mgdb2026"]["namFounders"]))
print("\nwrote:")
for p in ["data/accessions.tsv","data/projects.tsv","js/accessions.catalog.json","js/accessions.catalog.js"]:
    full = f"/home/claude/out/{p}"
    print(f"  {p:32s} {os.path.getsize(full):>9,d} bytes")
