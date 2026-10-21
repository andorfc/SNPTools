#!/usr/bin/env python3
"""build_score_cal.py -- the score calibration behind the composite AI score and its priority tiers.

    python3 tools/build_score_cal.py [--jobs N] [store.h5 ...]      writes data/qc/<family>.scorecal.json
                                                                    and js/<prefix>.scorecal.js
With no store named, the ten stores hdf5/grin2026/zmgrin2026_chr{1..10}_impute.h5 are used (read-only).

The composite score of a variant is the mean of the standardized values (z = (score - mean) / sd) of
the scores it has among PlantCAD1 (plantcad1_score), PlantCAD2 one-pass (plantcad2_onepass_score), ESM1
(ESM1_score) and ESM2 (ESM2_score). The reference distribution of each score, and of the composite, is
the release's missense sites that carry all four scores (every site of the stores, carriers or not;
the ESM scores of a site belong to its first missense consequence, as the app reads them). Standardizing
puts the four scores, whose scales differ (median missense PlantCAD1 about -0.9, PlantCAD2 one-pass -3.6,
ESM -4), on one scale before they are averaged, so no one model dominates.

Priority tiers are percentiles of the composite among those sites: TOP at or below the 5th percentile,
HIGH the 25th, MODERATE the 80th, LOW above (the proportions the former raw thresholds -7 / -4 / -1 on
mean(PlantCAD1, ESM1) gave). The cutoffs are stored as numbers, so the app and the Python mirrors
(tools/build_snpcurate.py, the harness) compare against the same values; the stored percentile table
(0..100) gives the percentile of any composite value for display.

Output (both hold the same object):
  data/qc/<family>.scorecal.json
  js/<prefix>.scorecal.js        window.SNP_SCORE_CAL[<family>] for the browser (index.html loads it)
Needs h5py + numpy.
"""
import argparse, datetime, json, os, re, sys, time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCORES = (('pc1', 'plantcad1_score'), ('pc2op', 'plantcad2_onepass_score'), ('esm1', 'ESM1_score'), ('esm2', 'ESM2_score'))
TIERS = (('top', 5), ('high', 25), ('moderate', 80))
PAT = {k: re.compile(r'(?:^|;)' + re.escape(key) + r'=([^;,]+)') for k, key in SCORES}
TYPE_PAT = re.compile(r'(?:^|;)TYPE=([^;]*)')


def scan_one(path):
    """Every missense site of one store that carries all four scores: an (n, 4) float array."""
    import h5py
    rows = []
    with h5py.File(path, 'r') as f:
        d = f['INFO']; n = d.shape[0]; step = 200000
        for a in range(0, n, step):
            for s in d[a:a + step]:
                s = s.decode() if isinstance(s, bytes) else str(s)
                if 'missense' not in s: continue
                m = TYPE_PAT.search(s)
                if not m or 'missense' not in m.group(1): continue
                v = []
                for k, _ in SCORES:
                    mm = PAT[k].search(s)
                    if not mm or mm.group(1) in ('.', ''): break
                    try: v.append(float(mm.group(1)))
                    except ValueError: break
                if len(v) == 4: rows.append(v)
    return os.path.basename(path), n, np.array(rows, dtype=np.float64).reshape(-1, 4)


def build(stores, jobs):
    t0 = time.time()
    parts = []
    if jobs <= 1 or len(stores) == 1:
        for p in stores: parts.append(scan_one(p))
    else:
        with ProcessPoolExecutor(max_workers=min(jobs, len(stores))) as ex:
            parts = list(ex.map(scan_one, stores))
    per_store = {name: {'sites': int(n), 'missense_scored': int(len(M))} for name, n, M in parts}
    M = np.concatenate([x[2] for x in parts]) if parts else np.zeros((0, 4))
    if len(M) < 100: raise SystemExit(f'only {len(M)} missense sites with all four scores: not enough for a calibration')
    mean, sd = M.mean(0), M.std(0)
    Z = (M - mean) / sd
    comp = Z.mean(1)
    pct = np.percentile(comp, np.arange(101))
    cal = {
        'format': 'snptools-scorecal-1',
        'built': datetime.date.today().isoformat(),
        'stores': per_store,
        'n_sites': int(sum(v['sites'] for v in per_store.values())),
        'n_reference': int(len(M)),
        'reference': 'missense sites of the stores carrying all four scores',
        'scores': {k: {'key': key, 'mean': round(float(mean[i]), 4), 'sd': round(float(sd[i]), 4),
                       'median': round(float(np.median(M[:, i])), 4)} for i, (k, key) in enumerate(SCORES)},
        'composite': {
            'of': [k for k, _ in SCORES],
            'definition': 'mean of (score - mean) / sd over the scores the variant has',
            'median': round(float(np.median(comp)), 4),
            'percentiles': [round(float(x), 4) for x in pct],
            'tiers': {name: round(float(np.percentile(comp, q)), 4) for name, q in TIERS},
            'tier_percentiles': {name: q for name, q in TIERS},
        },
        'seconds': round(time.time() - t0, 1),
    }
    return cal


def write(root, family, prefix, cal):
    jp = os.path.join(root, 'data', 'qc', f'{family}.scorecal.json')
    os.makedirs(os.path.dirname(jp), exist_ok=True)
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(cal, fh, indent=1); fh.write('\n')
    jsp = os.path.join(root, 'js', f'{prefix}.scorecal.js')
    with open(jsp, 'w', encoding='utf-8') as fh:
        fh.write(f'/* Generated by tools/build_score_cal.py on {cal["built"]} from the {len(cal["stores"])} stores of the release:\n'
                 f'   the reference distributions of PlantCAD1, PlantCAD2 one-pass, ESM1 and ESM2 over the release\'s\n'
                 f'   {cal["n_reference"]:,} missense sites carrying all four, the composite score\'s percentiles and the priority\n'
                 f'   tier cutoffs. Do not edit by hand: run `make score-cal`. */\n')
        fh.write('window.SNP_SCORE_CAL = window.SNP_SCORE_CAL || {};\n')
        fh.write(f'window.SNP_SCORE_CAL[{json.dumps(family)}] = ' + json.dumps(cal, separators=(',', ':')) + ';\n')
    return jp, jsp


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('stores', nargs='*', help='variant stores (default: the ten zmgrin2026 Imputed stores)')
    ap.add_argument('--root', default=ROOT)
    ap.add_argument('--family', default='zmgrin2026')
    ap.add_argument('--prefix', default='zmgrin')
    ap.add_argument('--jobs', type=int, default=1)
    a = ap.parse_args()
    stores = a.stores or [os.path.join(a.root, 'hdf5', 'grin2026', f'zmgrin2026_chr{c}_impute.h5') for c in range(1, 11)]
    for p in stores:
        if not os.path.isfile(p): raise SystemExit(f'no store: {p}')
    cal = build(stores, a.jobs)
    jp, jsp = write(a.root, a.family, a.prefix, cal)
    s = cal['scores']; t = cal['composite']['tiers']
    print(f"{cal['n_reference']:,} missense reference sites of {cal['n_sites']:,} in {len(stores)} stores ({cal['seconds']} s)")
    for k in s: print(f"  {k:<6} mean {s[k]['mean']:8.3f}  sd {s[k]['sd']:6.3f}  median {s[k]['median']:7.2f}")
    print(f"  composite median {cal['composite']['median']}; tiers TOP <= {t['top']}, HIGH <= {t['high']}, MODERATE <= {t['moderate']}")
    print('wrote', os.path.relpath(jp, a.root), os.path.relpath(jsp, a.root))


if __name__ == '__main__':
    main()
