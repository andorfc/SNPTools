#!/usr/bin/env python3
"""check_correctness.py -- the correctness fixes of 2026-10-21 (the correctness scenario of run_scenarios.js,
and SNPGeo's country nodes in the snpgeo_gene scenario).

    python3 check_correctness.py <results.json> <root>

Exit 1 on a failure."""
import json, re, sys

res, root = sys.argv[1:3]
R = json.load(open(res)); K = R.get('correctness')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not K:
    print('FAIL no correctness scenario in results.json'); sys.exit(1)
m = K['metric']
check((m['AB'], m['AC'], m['BC']) == (0.5, 0, 0.5) and (m['treeAB'], m['treeAC'], m['treeBC']) == (0.5, 0, 0.5),
      f"SNPCompare similarity is allele sharing and equals 1 - SNPTree's distance: AB {m['AB']} AC {m['AC']} BC {m['BC']} (tree {m['treeAB']} {m['treeAC']} {m['treeBC']})")
if K.get('gAvail'):
    r = K['race']
    check(sorted(K['pendingKeys']) == ['ZmG_B73', 'ZmG_B97'] and r['ran'] and r['ranFocal'] == 'ZmG_B97' and r['b73'] == 0.99,
          f"both-scope race: focal B73's fetch answering after focal B97's leaves B97's rows ({r})")
else:
    print('skip both-scope race: no genome-wide matrix installed')
x = K['matrix']
check(x['first'] == ['ZmG_B73', 'ZmG_B97', 'ZmG_MO17'] and x['second'] == ['ZmG_B73', 'ZmG_B97', 'ZmG_OH43'] and x['firstAC'] == 1 and x['secondAB'] == 0,
      f"SNPMatrix rebuilds for a second set of the same size: labels {x['second']}, A-B distance {x['secondAB']} (first set A-C {x['firstAC']})")
t = K['tree']
check('ZmG_MO17' in t['text1'] and t['builtAfterNew']['built'] is None and t['builtAfterNew']['text'] == '' and 'ZmG_OH43' in t['text2'] and 'ZmG_MO17' not in t['text2'],
      f"SNPTree's text box: cleared on a new hand-off, then the new tree ({len(t['text1'])} -> {len(t['builtAfterNew']['text'])} -> {len(t['text2'])} chars)")
s = K['stale']
check(s['results'] is None and (s['chr'], s['start'], s['end']) == ('chr4', 1000, 2000),
      f"SNPVersity: a hand-off during a running query drops that query's answer ({s})")
e = K['effect']
check(e['fEffect'] == 'all' and e['shown'] == 1, f"the effect filter returns to all effects when the new result lacks the picked effect ({e})")
# SNPGeo: one node per country, none pooled under '???'; the 30 crosswalk lines join their countries
cat = open(f'{root}/js/zmgrin.catalog.js', encoding='utf-8').read()
fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', cat, re.S).group(1))
us = sum(1 for p in fam['projects'] for g in p['groups'] for a in g['accessions'] if a.get('country') == 'United States')
G = (R.get('snpgeo_gene') or {}).get('allStats')
if G:
    c = G[0]['c']
    check('???' not in c and not any(k.startswith('name:') and c[k]['total'] > 5 for k in c) and c.get('USA', {}).get('total') == us,
          f"SNPGeo: no '???' node; United States has one node of {c.get('USA', {}).get('total')} lines (catalogue {us}); keys {sorted(c)[:6]}...")
else:
    print('skip SNPGeo nodes: no snpgeo_gene scenario')
check(not K['errors'], f"no console errors ({len(K['errors'])})" + (': ' + '; '.join(x[:120] for x in K['errors'][:3]) if K['errors'] else ''))
print(f'check_correctness: {bad} failure(s)')
sys.exit(1 if bad else 0)
