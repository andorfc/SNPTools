#!/usr/bin/env python3
"""check_snptrait.py <results.json> <site_root> -- recompute SNPTrait facet counts and the
compound-filter result of run_scenarios.js from js/zmgrin.catalog.js. Stdlib only."""
import json, re, sys
res, root = sys.argv[1:3]
T = json.load(open(res))['snptrait']
cat = open(f'{root}/js/zmgrin.catalog.js', encoding='utf-8').read()
fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', cat, re.S).group(1))
sch = json.loads(re.search(r'SNPTRAIT_SCHEMA\["zmgrin2026"\]\s*=\s*(\{.*?\});', cat, re.S).group(1))
rows = [a for p in fam['projects'] for g in p['groups'] for a in g['accessions']]
val = lambda a, k: a.get(k) if a.get(k) not in (None, '') else 'Unknown'
bad = 0
for k, _ in sch['facets']:
    exp = {}
    for a in rows: exp[val(a, k)] = exp.get(val(a, k), 0) + 1
    if exp != T['facetCounts'][k]: bad += 1; print('facet differs:', k)
f = [a for a in rows if val(a, 'subpop') == 'SS' and val(a, 'inAmes282') == 'yes' and val(a, 'kernelType') == 'Dent']
if sorted(a['id'] for a in f) != T['filterSS_Ames_Dent']: bad += 1; print('compound filter differs')
hay = lambda a: ' '.join(str(a.get(k) or (a.get('label') or a['id'] if k == 'strain' else '')) for k in sch['search']).lower()
fi = [a for a in f if 'iowa' in hay(a)]
if sorted(a['id'] for a in fi) != T['filterPlusIowa']: bad += 1; print('search differs', len(fi), len(T['filterPlusIowa']))
tr = json.load(open(f'{root}/data/traits/zmgrin2026.traits.json'))
kw = sorted(a['id'] for a in rows if (tr['samples'].get(a['id'], {}).get('traits', {}).get('KERNEL-WEIGHT-1000') or [0, None])[1] is not None
            and 250 <= tr['samples'][a['id']]['traits']['KERNEL-WEIGHT-1000'][1] <= 300)
if kw != T['rangeKW']: bad += 1; print('range filter differs', len(kw), len(T['rangeKW']))
kwa = [i for i in kw if next(a for a in rows if a['id'] == i).get('panel') == 'Ames282']
if len(kwa) != T['rangeKW_Ames']: bad += 1; print('range x panel differs', len(kwa), T['rangeKW_Ames'])
# the Send dialog: replace pre-ticked; ticking add unticks replace; add keeps SNPVersity's
# selection and adds the SNPTrait lines (union), and says how many are new
D, B, H = T['dialog'], T['dialogBefore'], T['handoffAdd']
union = sorted(set(B['versity']) | set(B['trait']))
new = len(set(B['trait']) - set(B['versity']))
if not (D['open'] and D['replace'] and not D['add']): bad += 1; print('dialog defaults', D)
if D['afterAdd'] != {'replace': False, 'add': True, 'sendDisabled': False}: bad += 1; print('dialog after add', D['afterAdd'])
if f"replace the {len(B['versity'])} accessions" not in D['text'] or f"({new} new" not in D['text']: bad += 1; print('dialog text', D['text'])
if H['tool'] != 'snpversity' or not H['dialogClosed'] or H['selected'] != union: bad += 1; print('add hand-off', len(H['selected']), len(union))
if f"+{new} added" not in H['banner']: bad += 1; print('arrival banner', H['banner'])
if T['handoff']['selected'] != T['selectedAfterSelectVisible']: bad += 1; print('replace hand-off', T['handoff'])   # SNPTrait starts from SNPVersity's selection
print(f"Send dialog: replace {T['handoff']['selected']} lines; add {len(B['trait'])} to {len(B['versity'])} -> {len(H['selected'])} ({new} new)")
print(f'KW1000 in [250,300] g: {len(kw)} (x Ames282: {len(kwa)})')
print(f'rows={len(rows)} facets={len(sch["facets"])} SSxAmes282xDent={len(f)} +iowa={len(fi)} mismatches={bad}')
sys.exit(1 if bad else 0)
