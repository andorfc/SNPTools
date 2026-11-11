#!/usr/bin/env python3
"""check_ibs_files.py <site_root> <family> -- consistency of the installed genome-wide IBS files
in <site_root>/distance/<family>/ (skipped, exit 0, when the directory is absent).

ids.txt: no duplicates, exactly the catalogue ids of the family (js/*.catalog.js), any order.
similarity(_snp).csv, missing_pct(_snp).csv, allele_distance.csv, co_called_sites.csv (whichever
exist): len(ids) x len(ids), numeric, symmetric; similarity diagonal 1 and values in 0..1;
missing fractions in 0..1 (not percent); allele_distance diagonal 0. tree_*.nwk: every tip is an
id of ids.txt and every id is a tip. Standard library only."""
import glob, json, os, re, sys

root, fam = sys.argv[1:3]
D = os.path.join(root, 'distance', fam)
if not os.path.isdir(D):
    print(f'{D}: not installed, skipped'); sys.exit(0)
bad = 0
ids = [l.strip() for l in open(f'{D}/ids.txt') if l.strip()]
cat = set()
for f in glob.glob(f'{root}/js/*catalog.js'):
    for line in open(f, encoding='utf-8'):
        m = re.match(r'window\.SNP_CATALOG(?:\.families\["%s"\])?\s*=\s*(\{.*\});?\s*$' % fam, line)
        if not m:
            continue
        o = json.loads(m.group(1))
        node = o['families'][fam] if 'families' in o and fam in o.get('families', {}) else (o if 'projects' in o else None)
        if node:
            cat |= {a['id'] for p in node['projects'] for g in p['groups'] for a in g['accessions']}
if len(set(ids)) != len(ids):
    bad += 1; print('duplicate ids')
if set(ids) != cat:
    bad += 1; print('ids vs catalogue: not in catalogue', sorted(set(ids) - cat)[:5], '; catalogue not in ids', sorted(cat - set(ids))[:5])
n = len(ids); out = {}
for name in ('similarity', 'similarity_snp', 'missing_pct', 'missing_pct_snp', 'allele_distance', 'co_called_sites'):
    p = f'{D}/{name}.csv'
    if not os.path.isfile(p):
        continue
    M = [[float(x) for x in l.split(',')] for l in open(p) if l.strip()]
    if len(M) != n or any(len(r) != n for r in M):
        bad += 1; print(name, 'shape', len(M)); continue
    asym = max(abs(M[i][j] - M[j][i]) for i in range(n) for j in range(i + 1, n))
    lo = min(min(r) for r in M); hi = max(max(r) for r in M)
    diag = [M[i][i] for i in range(n)]
    if asym > 1e-9: bad += 1; print(name, 'asymmetric', asym)
    if name.startswith('similarity') and (min(diag) != 1 or max(diag) != 1 or lo < 0 or hi > 1): bad += 1; print(name, 'range/diag', lo, hi)
    if name.startswith('missing') and (lo < 0 or hi > 1): bad += 1; print(name, 'not a fraction', lo, hi)
    if name == 'allele_distance' and max(abs(d) for d in diag) > 0: bad += 1; print(name, 'diag not 0')
    out[name] = f'{n}x{n} [{lo:.4g}, {hi:.4g}]'
for t in sorted(glob.glob(f'{D}/tree_*.nwk')):
    tips = re.findall(r'[(,]\s*([^():,;\s]+)\s*:', open(t).read())
    if sorted(tips) != sorted(ids):
        bad += 1; print(os.path.basename(t), 'tips differ from ids', len(tips))
    out[os.path.basename(t)] = f'{len(tips)} tips'
print(f'{fam} IBS files: {n} ids = catalogue; ' + '; '.join(f'{k} {v}' for k, v in out.items()) + f'; mismatches={bad}')
sys.exit(1 if bad else 0)
