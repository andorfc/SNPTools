#!/usr/bin/env python3
"""verify_ibs_ids.py — sanity-check the accession order for the precomputed
distance matrices used by ibsCompare.php.

The matrices (distance/maizegdb_allchr_final_similarity.csv, ..._missing_pct.csv)
are headerless and symmetric, so they only make sense alongside an ordered list
of accession IDs — distance/ids.txt — where line k is the accession for row/col k.

This script checks two things:
  1) ID match: every id in ids.txt exists in the SNPVersity accession catalog
     (so the browser can join Project / SRA / Name and the focal picker works).
  2) Order sanity: replicates of the same line (same 'founder') should be each
     other's most-similar partners. If ids.txt is in the wrong order, this fails.

Usage:
  python3 verify_ibs_ids.py \
      --ids distance/ids.txt \
      --sim distance/maizegdb_allchr_final_similarity.csv \
      --catalog js/accessions.catalog.json \
      [--family mgdb2026] [--topk 5] [--sample 300]

--sim is optional; omit it to run only the fast ID-match check.
"""
import argparse, json, sys, random

ap = argparse.ArgumentParser()
ap.add_argument('--ids', required=True)
ap.add_argument('--sim')
ap.add_argument('--catalog', required=True)
ap.add_argument('--family', default='mgdb2026')
ap.add_argument('--topk', type=int, default=5)
ap.add_argument('--sample', type=int, default=300)
a = ap.parse_args()

ids = [l.strip() for l in open(a.ids) if l.strip()]
print(f"ids.txt: {len(ids)} accessions")

cat = json.load(open(a.catalog))
fam = cat['families'].get(a.family)
if not fam:
    print(f"!! family {a.family} not in catalog"); sys.exit(1)
founder = {}
for p in fam['projects']:
    for g in p['groups']:
        for acc in g['accessions']:
            founder[acc['id']] = acc['founder']
print(f"catalog {a.family}: {len(founder)} accessions")

# ---- 1) ID match ----
missing = [i for i in ids if i not in founder]
print(f"\n[1] ID match: {len(ids)-len(missing)}/{len(ids)} ids found in catalog")
if missing:
    print(f"    !! {len(missing)} ids NOT in catalog (first 8): {missing[:8]}")
    print("    -> the matrix ID strings must match SNPVersity Final_IDs exactly, or the")
    print("       browser join (project/SRA/name) will be blank for those rows.")
else:
    print("    all ids.txt entries match catalog IDs — good.")

if not a.sim:
    print("\n(skip order check: no --sim given)"); sys.exit(0)

# ---- 2) order sanity via replicate similarity ----
import csv
# founders that have >=2 reps within ids.txt order
from collections import defaultdict
byf = defaultdict(list)
for k, i in enumerate(ids):
    byf[founder.get(i, '?')].append(k)
rep_founders = {f: idxs for f, idxs in byf.items() if len(idxs) >= 2 and f != '?'}
print(f"\n[2] order check: {len(rep_founders)} lines have >=2 replicates")
if not rep_founders:
    print("    no replicated lines to test order with."); sys.exit(0)

# read only the rows we need
want = set()
sample_founders = list(rep_founders)
random.seed(0); random.shuffle(sample_founders)
for f in sample_founders:
    for k in rep_founders[f]:
        want.add(k)
        if len(want) >= a.sample: break
    if len(want) >= a.sample: break

rows = {}
with open(a.sim) as fh:
    for k, line in enumerate(fh):
        if k in want:
            rows[k] = [float(x) for x in line.rstrip('\n').split(',')]
        if len(rows) == len(want):
            break

hit = 0; tot = 0
for f, idxs in rep_founders.items():
    idxset = set(idxs)
    for k in idxs:
        if k not in rows: continue
        r = rows[k]
        order = sorted(range(len(r)), key=lambda j: r[j], reverse=True)
        topk = [j for j in order if j != k][:a.topk]
        tot += 1
        if any(j in idxset for j in topk):
            hit += 1
if tot:
    rate = 100*hit/tot
    print(f"    of {tot} tested replicate accessions, {hit} ({rate:.1f}%) have a same-line")
    print(f"    replicate within their top-{a.topk} most-similar partners.")
    print("    -> expect a HIGH rate (usually >90%) if ids.txt order is correct.")
    print("    -> a low rate means ids.txt is likely mis-ordered relative to the matrix.")
