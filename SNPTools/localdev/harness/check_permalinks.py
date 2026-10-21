#!/usr/bin/env python3
"""check_permalinks.py -- the hash router (core.js setRoute / applyRoute): the permalinks scenario of
run_scenarios.js (results.json "permalinks").

    python3 check_permalinks.py <results.json> <root>

Exit 1 on a failure."""
import json, sys

res, root = sys.argv[1:3]
P = json.load(open(res)).get('permalinks')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not P:
    print('FAIL no permalinks scenario in results.json'); sys.exit(1)
b = P['boot']
check(b['hash'] == '#snpversity' and b['tool'] == 'snpversity' and b['title'] == 'SNPMaize · SNPVersity',
      f"first load: {b['hash']}, title \"{b['title']}\"")
c = P['curate']
check(c['hash'] == '#snpcurate/SC0003' and c['title'].startswith('SNPMaize · SNPCurate · DGAT1-2'),
      f"an open record writes its hash: {c['hash']}, title \"{c['title']}\"")
f = P['func']
check(f['tool'] == 'snpfunction' and f['hash'] == '#snpfunction/Zm00001eb277490' and f['gene'] == 'Zm00001eb277490'
      and f['title'] == 'SNPMaize · SNPFunction · Zm00001eb277490',
      f"a pasted #snpfunction/<gene> opens the gene: tool {f['tool']}, gene box {f['gene']!r}, title \"{f['title']}\"")
k = P['back']
check(k['tool'] == 'snpcurate' and k['record'] is False and k['hash'] == '#snpcurate',
      f"#snpcurate (Back from a record) shows the table without a record: tool {k['tool']}, record {k['record']}")
v = P['versity']
q = v.get('q') or {}
check(v['tool'] == 'snpversity' and (q.get('dataset'), q.get('chr'), q.get('lo'), q.get('hi')) == ('zmgrin2026_imp', 'chr4', 46648300, 46648450)
      and v['hash'] == '#snpversity?ds=zmgrin2026_imp&chr=chr4&start=46648300&end=46648450' and (v.get('n') or 0) >= 1
      and v['title'].startswith('SNPMaize · SNPVersity · chr4:'),
      f"#snpversity?ds=&chr=&start=&end= runs the query: {q}, {v.get('n')} rows, hash {v['hash']}, title \"{v['title']}\"")
r = P['rerun']
check(r['requests'] == v['requests'] and r['tool'] == 'snpversity',
      f"returning to the same region does not query again ({r['requests']} processForm requests before and after)")
u = P['unknown']
check(u['tool'] == 'snpversity' and u['hash'] == '#nosuchtool', f"an unknown hash is ignored: tool {u['tool']}")
g = P['gwas']
check(g['tool'] == 'snpgwas' and g['active'] == 'tibbscortes2024_nam_asi_intcp' and g['hash'] == '#snpgwas/tibbscortes2024_nam_asi_intcp',
      f"#snpgwas/<id> opens that dataset: active {g['active']}, hash {g['hash']}")
e = P['errors']
check(not e, f"no console errors in the scenario ({len(e)})" + (': ' + '; '.join(x[:120] for x in e[:3]) if e else ''))
print(f'check_permalinks: {bad} failure(s)')
sys.exit(1 if bad else 0)
