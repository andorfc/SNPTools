#!/usr/bin/env python3
"""Independent check of PLINK --r2 values: recompute genotypic r2 (pairwise-complete allele-count correlation)
for 300 random pairs from the first 300,000 pair lines, from the PLINK .bed of the same chromosome."""
import sys, json, random
import numpy as np
O, C = sys.argv[1], sys.argv[2]
bim = [l.split('\t')[1] for l in open(f'{O}/g.bim')]
idx = {v: i for i, v in enumerate(bim)}
n = sum(1 for _ in open(f'{O}/g.fam'))
bpv = (n + 3) // 4
raw = np.fromfile(f'{O}/g.bed', dtype=np.uint8)[3:]
def geno(i):
    b = raw[i * bpv:(i + 1) * bpv]
    g = np.stack([(b >> (2 * k)) & 3 for k in range(4)], 1).reshape(-1)[:n]
    # PLINK .bed: 00 hom A1, 01 missing, 10 het, 11 hom A2  -> allele-count, -1 missing
    return np.select([g == 0, g == 2, g == 3], [2, 1, 0], -1)
lines = [l.split() for l in open(f'{O}/{C}.pairs_sample.txt')]
random.seed(3); pick = random.sample(lines, min(300, len(lines)))
diffs = []; nan = 0
for t in pick:
    x, y = geno(idx[t[2]]), geno(idx[t[5]]); ok = (x >= 0) & (y >= 0)
    x, y = x[ok].astype(float), y[ok].astype(float)
    if x.std() == 0 or y.std() == 0:
        nan += 1; continue
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    diffs.append(abs(r2 - float(t[6])))
d = np.array(diffs)
print(json.dumps(dict(pairs=len(pick), compared=int(len(d)), monomorphic_skipped=nan,
                      max_absdiff=float(d.max()) if len(d) else None, n_absdiff_gt_1e4=int((d > 1e-4).sum()))))

# --- part 2: recompute MAXR2 from the raw pair sample (plink output lines) for variants whose whole
# +/-5 kb neighbourhood lies inside the sample, and compare with the reduced chrN.maxr2.tsv.gz
import gzip, math
last_a = max(int(t[1]) for t in lines)
first_a = min(int(t[1]) for t in lines)
mx = {}
for t in lines:
    d = abs(int(t[4]) - int(t[1]))
    if d < 400 or t[6] in ('nan', '-nan'):
        continue
    v = float(t[6])
    for k in (t[2], t[5]):
        if k not in mx or v > mx[k]:
            mx[k] = v
got = {}; n_tab = 0
with gzip.open(f'{O}/{C}.maxr2.tsv.gz', 'rt') as fh:
    next(fh)
    for l in fh:
        c, p, r, a, m = l.rstrip('\n').split('\t')
        p = int(p)
        if p > last_a - 5000: break
        if p < first_a + 5000: continue
        got[f'{c}:{p}:{r}:{a}'] = m; n_tab += 1
agree = sum(1 for k, m in got.items() if (m == '' and k not in mx) or (m != '' and k in mx and abs(float(m) - mx[k]) < 1e-9))
print(json.dumps(dict(reducer_check_sites=n_tab, reducer_agree=agree, region=f'{C}:{first_a + 5000}-{last_a - 5000}')))
