#!/usr/bin/env python3
"""check_performance.py -- the performance batch of 2026-10-22: scripts loaded on demand, SNPFold's
per-gene sites files, SNPGeo's paged table and cheap selection, SNPCompare's clustering, PCoA and
matrix gate (the performance scenario of run_scenarios.js, and the pergene scenario's SNPFold fetches).

    python3 check_performance.py <results.json> <root>

Exit 1 on a failure."""
import json, os, re, sys

res, root = sys.argv[1:3]
R = json.load(open(res)); P = R.get('performance')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not P:
    print('FAIL no performance scenario in results.json'); sys.exit(1)
html = open(f'{root}/index.html', encoding='utf-8').read()
check('data/structures/' not in html and 'js/pe/' not in html and 'd3js.org' not in html and 'paneffect-engine.js' not in html,
      'index.html loads no structure file, no PanEffect script and no d3 up front')
check(re.search(r"<script[^>]*\bdefer\b[^>]*issuecollector", html) is not None, 'the Jira collector script is deferred')
s = P['scripts']
check(s['d3Before'] == 'undefined' and s['d3AfterGeo'] == 'object' and s['peBefore'] == 'undefined' and s['peAfter'] == 'function',
      f"d3 appears when SNPGeo draws, the PanEffect engine when PanEffect opens ({s})")
f = ((R.get('siteqcui') or {}).get('on') or {}).get('fold') or {}   # the site-QC UI scenario's SNPFold pass (full stores)
sf = f.get('sitesFetch')
if sf is None:
    print('skip SNPFold sites fetches: the site-QC UI scenario did not run (full stores absent)')
else:
    check(any(u.endswith('data/sites/by_gene/MANIFEST.json') for u in sf) and any(u.endswith('data/sites/by_gene/Zm00001eb406050.tsv') for u in sf)
          and not any(u.endswith('results.sites.tsv') for u in sf), f"SNPFold fetched the gene's own sites file, not the 55 MB table: {sf}")
    man = f'{root}/data/sites/by_gene/MANIFEST.json'
    if os.path.isfile(man) and os.path.isfile(f'{root}/data/results.sites.tsv'):
        m = json.load(open(man)); n = 0
        big = {}
        for line in open(f'{root}/data/results.sites.tsv', encoding='utf-8'):
            g = line.split('\t', 1)[0].rsplit('_', 1)[0]
            if g in ('Zm00001eb406050', 'Zm00001eb062040', 'Zm00001eb174590'): big.setdefault(g, []).append(line)
        ok = True
        for g, lines in big.items():
            p = f'{root}/data/sites/by_gene/{g}.tsv'
            if not os.path.isfile(p) or open(p, encoding='utf-8').readlines() != lines: ok = False; print('   differs', g)
            n += len(lines)
        check(ok and m['genes'] > 10000 and m['format'] == 'snptools-sites-by-gene-1', f"per-gene files equal the table's rows for 3 genes ({n} rows); manifest {m['genes']} genes")
g = P['geo']
check(g['rowsOnPage'] == 50 and g['pager'].startswith('Variants 1–50 of') and g['pages'] >= 2,
      f"SNPGeo pages the variant table: {g['rowsOnPage']} rows on page 1, {g['pages']} pages, pager '{g['pager']}'")
check(g['sel60']['page'] == 1 and g['sel60']['active'] == 60 and g['sel61']['sameTbody'] is True and g['sel61']['active'] == 61 and g['sel61']['nActive'] == 1,
      f"selecting a variant on another page shows that page; the next selection moves the mark without rebuilding ({g['sel60']}, {g['sel61']})")
check(g['key']['active'] == 62 and g['key']['sameTbody'] is True, f"ArrowDown moves the mark in place ({g['key']})")
c = P['cluster']
check(c['heightsMatch'] and c['orderIsPermutation'], f"SNPCompare's UPGMA: merge heights equal a reference UPGMA on a random 12-line matrix; the order is a permutation ({c['nMerges']} merges)")
q = P['pcoa']
check(abs(q['line']['v1'] - 1) < 1e-6 and q['line']['v2'] < 1e-6 and q['line']['spaced'] and q['line']['nonEuclid'] < 1e-6,
      f"PCoA on four collinear lines: axis 1 carries {q['line']['v1']:.4f} of the variation, axis 2 {q['line']['v2']:.2e}, points evenly spaced")
check(q['square']['finite'] and q['square']['nonEuclid'] > 0.01 and q['square']['v1'] + q['square']['v2'] <= 1 + 1e-9,
      f"PCoA on a non-Euclidean configuration: finite points, non-Euclidean share {q['square']['nonEuclid']:.3f}, axes {q['square']['v1']:.3f} + {q['square']['v2']:.3f}")
m = P['matrixGate']
check(m['notice'] and 'Draw anyway' in m['noticeText'] and m['afterForce'], f"the matrix of {m['n']} lines waits behind a notice; Draw anyway draws it ({m['noticeText'][:60]!r})")
check(not P['errors'], f"no console errors ({len(P['errors'])})" + (': ' + '; '.join(e[:100] for e in P['errors'][:3]) if P['errors'] else ''))
print(f'check_performance: {bad} failure(s)')
sys.exit(1 if bad else 0)
