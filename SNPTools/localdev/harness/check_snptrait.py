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
if T['selectedAfterSelectVisible'] != len(f): bad += 1; print('SNPTrait did not start empty', T['selectedAfterSelectVisible'], len(f))
if T['handoff']['selected'] != len(f): bad += 1; print('replace hand-off', T['handoff'])
print(f"Send dialog: replace {T['handoff']['selected']} lines; add {len(B['trait'])} to {len(B['versity'])} -> {len(H['selected'])} ({new} new)")
# SNPVersity step 3 and the SNPTrait drawer
Q = json.load(open(res))['step3']
nam = sorted(a['id'] for a in rows if a.get('namFounder'))
ames = sorted(a['id'] for a in rows if a.get('inAmes282') == 'yes')
if Q['before'] != nam: bad += 1; print('step 3 start (NAM founders)', len(Q['before']), len(nam))
if Q['initial']['oldPicker'] or Q['initial']['count'] != str(len(nam)): bad += 1; print('step 3 initial', Q['initial'])
picks = dict((t.rsplit(' ', 1)[0], on) for t, on in Q['initial']['picks'] if t != 'Clear')
if not picks.get('NAM founders + B73') or f"Ames 282 {len(ames)}" not in [t for t, _ in Q['initial']['picks']]:
    bad += 1; print('quick picks', Q['initial']['picks'])
if Q['pick']['selected'] != ames or not Q['pick']['on'] or 'Ames 282' not in Q['pick']['say'] or 'Undo' not in Q['pick']['say']:
    bad += 1; print('quick pick Ames 282', len(Q['pick']['selected']), len(ames), Q['pick']['say'])
if Q['undo'] != nam: bad += 1; print('undo after quick pick', len(Q['undo']))
D = Q['drawer']
if not D['open'] or D['mode'] != 'drawer' or D['draft'] != nam or D['facetBlocks'] < 10: bad += 1; print('drawer open', {k: v for k, v in D.items() if k != 'visible'})
if D['visible'] != sorted(a['id'] for a in f): bad += 1; print('drawer filter', len(D['visible']), len(f))
exp_applied = sorted(set(nam) | set(D['visible']))
added = len(set(D['visible']) - set(nam))
if f"+{added} vs the {len(nam)} selected in SNPVersity" not in D['foot']: bad += 1; print('drawer footer', D['foot'])
A = Q['applied']
if A['selected'] != exp_applied or not A['closed'] or A['mode'] != 'page' or f"+{added}" not in A['say']: bad += 1; print('drawer apply', len(A['selected']), len(exp_applied), A['say'])
if Q['reopenFilters'] != ['SS']: bad += 1; print('drawer filters not kept', Q['reopenFilters'])
if Q['afterCancel'] != {'selected': exp_applied, 'closed': True}: bad += 1; print('Escape changed the selection', len(Q['afterCancel']['selected']))
if Q['pageAfter'] != {'selected': Q['pageSel'], 'filtered': False}: bad += 1; print('SNPTrait page state disturbed', Q['pageAfter'])
mo17 = [a['id'] for a in rows if a.get('strain') == 'Mo17']
if Q['paste']['selected'] != sorted(['ZmG_B73'] + mo17) or Q['paste']['unmatched'] != ['not-a-line']: bad += 1; print('pasted list', Q['paste'])
print(f"step 3: quick pick Ames 282 = {len(ames)}, undo -> {len(Q['undo'])}; drawer +{added} -> {len(exp_applied)}; "
      f"pasted 'PI 550473', 'Mo17' -> {', '.join(Q['paste']['selected'])}")
print(f'KW1000 in [250,300] g: {len(kw)} (x Ames282: {len(kwa)})')
print(f'rows={len(rows)} facets={len(sch["facets"])} SSxAmes282xDent={len(f)} +iowa={len(fi)} mismatches={bad}')
sys.exit(1 if bad else 0)
