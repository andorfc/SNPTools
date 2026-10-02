#!/usr/bin/env python3
"""check_global_compare.py <results.json> -- SNPCompare genome-wide scope per dataset family
and the SNPTree genome-wide Newick links, against the synthetic files the harness wrote
(make_synthetic_distance.py; the directory is recorded in results.json).

(a) zmgrin2026: the Genome-wide scope is enabled only when distance/zmgrin2026/ exists; the
    focal row returned for ZmG_B73 (all sites) and ZmG_CML103 (SNPs only) equals the matrix row
    of the files, missing as percent; every one of the 933 ids is present.
(b) mgdb2026: the unchanged layout (distance/maizegdb_allchr_final_*.csv) still serves, and the
    request carries no dataset parameter (same URL as before). Checked only when run_scenarios.js
    ran it, i.e. when MaizeGDB 2026 is offered (off in the GRIN-linked initial release).
(c) SNPTree offers exactly the trees present, and the downloads equal the files.
Standard library only."""
import json, sys

R = json.load(open(sys.argv[1]))['compare']
dd = R['distanceDir']
bad = 0


def load(d, sim, mis):
    ids = [l for l in open(f'{d}/ids.txt').read().split('\n') if l]
    S = [l.split(',') for l in open(f'{d}/{sim}').read().split('\n') if l]
    M = [l.split(',') for l in open(f'{d}/{mis}').read().split('\n') if l]
    return ids, S, M


def check(tag, res, d, sim, mis):
    global bad
    ids, S, M = load(d, sim, mis)
    if not res.get('gAvail') or res.get('scopeDisabled'):
        bad += 1; print(tag, 'genome-wide scope not enabled'); return
    i = ids.index(res['focal'])
    exp = {ids[j]: (float(S[i][j]), float(M[i][j]) * 100) for j in range(len(ids))}
    got = {r[0]: (r[1], r[2]) for r in res['rows']}
    if set(got) != set(exp):
        bad += 1; print(tag, 'id set differs', len(got), len(exp))
    d_ = [k for k in exp if k in got and (abs(got[k][0] - exp[k][0]) > 1e-9 or abs(got[k][1] - exp[k][1]) > 1e-6)]
    if d_:
        bad += 1; print(tag, 'values differ for', len(d_), d_[:3])
    return len(exp)


n1 = check('zmgrin2026/all', R['zmgrin_all'], f'{dd}/zmgrin2026', 'similarity.csv', 'missing_pct.csv')
n2 = check('zmgrin2026/snp', R['zmgrin_snp'], f'{dd}/zmgrin2026', 'similarity_snp.csv', 'missing_pct_snp.csv')
if not R['zmgrin_all'].get('sitesSelect'):
    bad += 1; print('SNP-only selector not offered')
n3 = None
if 'mgdb' in R:
    n3 = check('mgdb2026', R['mgdb'], dd, 'maizegdb_allchr_final_similarity.csv', 'maizegdb_allchr_final_missing_pct.csv')
    mq = [q for q in R['requests'] if 'focal=' in q and R['mgdb']['focal'] in q]
    if not mq or any('dataset=' in q for q in mq):
        bad += 1; print('mgdb2026 request URL changed', mq)
if not any('dataset=zmgrin2026' in q and 'sites=snp' in q for q in R['requests']):
    bad += 1; print('no SNP-only zmgrin2026 request', R['requests'])
if 'real' in R:   # installed matrices (./distance/zmgrin2026/)
    import os
    real = os.path.join(sys.argv[2] if len(sys.argv) > 2 else '../..', 'distance', 'zmgrin2026')
    nr = check('zmgrin2026/installed', R['real'], real, 'similarity.csv', 'missing_pct.csv')
    check('zmgrin2026/installed MO17', R['real_mo17'], real, 'similarity.csv', 'missing_pct.csv')
    check('zmgrin2026/installed snp', R['real_snp'], real, 'similarity_snp.csv', 'missing_pct_snp.csv')
    top = lambda r: [x[0] for x in sorted(r['rows'], key=lambda x: -x[1])[1:3]]
    print('installed matrices: ZmG_B73 nearest', top(R['real']), '; ZmG_MO17 nearest', top(R['real_mo17'])[:1],
          '; NJ tree', R['real_tree']['njBytes'], 'bytes')
nf = R['zmgrin_nofiles']
if nf.get('gAvail') or not nf.get('scopeDisabled') or 'No precomputed genome-wide IBS matrix' not in nf.get('note', ''):
    bad += 1; print('no-files case', nf)
T = R['tree']
if [l[1] for l in T['links']] != ['ibsCompare.php?tree=nj&dataset=zmgrin2026', 'ibsCompare.php?tree=upgma&dataset=zmgrin2026']:
    bad += 1; print('tree links', T['links'])
for t in ('nj', 'upgma'):
    if T[t] != open(f'{dd}/zmgrin2026/tree_{t}.nwk').read():
        bad += 1; print('tree download differs', t)
if R['tree_nofiles'].strip():
    bad += 1; print('tree card shown without files')
print(f'SNPCompare genome-wide: zmgrin2026 {n1} ids (all sites, focal {R["zmgrin_all"].get("focal")}), '
      f'{n2} ids (SNPs, focal {R["zmgrin_snp"].get("focal")}); mgdb2026 legacy layout {n3 if n3 is not None else "not offered"}; '
      f'disabled without files: {not nf.get("gAvail")}; SNPTree Newick downloads: {", ".join(l[0] for l in T["links"])}; mismatches={bad}')
sys.exit(1 if bad else 0)
