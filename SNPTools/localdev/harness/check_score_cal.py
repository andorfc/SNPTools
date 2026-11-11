#!/usr/bin/env python3
"""check_score_cal.py -- the composite AI score's calibration (tools/build_score_cal.py).

    python3 check_score_cal.py <results.json> <root>

Checks that js/zmgrin.scorecal.js and data/qc/zmgrin2026.scorecal.json hold the same object, that the
tier cutoffs are the stated percentiles of the stored percentile table, that the composite's reference
is standardized (median near 0, every sd > 0), and, when the chr2 store is installed, that the builder's
own scan of it gives the count the file records and that the tiers cut chr2's missense composites at
about their stated shares. With REQUIRE_FULL=1 a missing store fails instead of being skipped.
Exit 1 on a failure."""
import json, os, re, sys

res, root = sys.argv[1:3]
sys.path.insert(0, os.path.join(root, 'tools'))
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

J = json.load(open(f'{root}/data/qc/zmgrin2026.scorecal.json', encoding='utf-8'))
js = open(f'{root}/js/zmgrin.scorecal.js', encoding='utf-8').read()
m = re.search(r'window\.SNP_SCORE_CAL\["zmgrin2026"\] = (\{.*\});\n', js, re.S)
check(m is not None and json.loads(m.group(1)) == J, 'js/zmgrin.scorecal.js holds the same object as data/qc/zmgrin2026.scorecal.json')
C = J['composite']; P = C['percentiles']; T = C['tiers']
check(J.get('format') == 'snptools-scorecal-1' and C['of'] == ['pc1', 'pc2op', 'esm1', 'esm2'], f"format {J.get('format')}, composite of {C['of']}")
check(len(P) == 101 and all(P[i] <= P[i + 1] for i in range(100)), 'percentile table: 101 non-decreasing values')
check(all(T[k] == P[q] for k, q in C['tier_percentiles'].items()) and C['tier_percentiles'] == {'top': 5, 'high': 25, 'moderate': 80},
      f"tiers are the 5th / 25th / 80th percentiles: TOP <= {T['top']}, HIGH <= {T['high']}, MODERATE <= {T['moderate']}")
check(abs(C['median']) < 0.2 and all(v['sd'] > 0 for v in J['scores'].values()) and J['n_reference'] >= 1000,
      f"reference standardized: composite median {C['median']}, {J['n_reference']:,} missense sites, sds {[v['sd'] for v in J['scores'].values()]}")
check(set(J['scores']) == set(C['of']) and all(J['scores'][k]['key'] == key for k, key in
      (('pc1', 'plantcad1_score'), ('pc2op', 'plantcad2_onepass_score'), ('esm1', 'ESM1_score'), ('esm2', 'ESM2_score'))), 'the four scores and their INFO keys')

store = f'{root}/hdf5/grin2026/zmgrin2026_chr2_impute.h5'
if os.path.isfile(store) and os.path.getsize(store) > 1e9:
    import numpy as np
    from build_score_cal import scan_one
    name, n, M = scan_one(store)
    rec = J['stores'].get(name, {})
    check((rec.get('sites'), rec.get('missense_scored')) == (n, len(M)), f'{name}: {n:,} sites, {len(M):,} scored missense sites, as the file records')
    mean = np.array([J['scores'][k]['mean'] for k in C['of']]); sd = np.array([J['scores'][k]['sd'] for k in C['of']])
    comp = ((M - mean) / sd).mean(1)
    shares = {k: float((comp <= T[k]).mean() * 100) for k in T}
    check(all(abs(shares[k] - q) <= 3 for k, q in C['tier_percentiles'].items()),
          'chr2 shares at or below the cutoffs: ' + ', '.join(f"{k} {shares[k]:.1f}% (release {q}%)" for k, q in C['tier_percentiles'].items()))
elif os.environ.get('REQUIRE_FULL') == '1':
    check(False, f'full store {store} is required (REQUIRE_FULL=1)')
else:
    print('skip the chr2 store is not installed: counts and tier shares not re-derived')

print(f'check_score_cal: {bad} failure(s)')
sys.exit(1 if bad else 0)
