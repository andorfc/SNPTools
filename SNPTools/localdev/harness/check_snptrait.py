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
# with js/zmgrin.lineqc.js loaded the page's schema gains a "Sample heterozygosity" facet (SNPTrait
# adds it to a copy of the generated schema): count it like the others
FACETS = list(sch['facets'])
if 'sampleQC' in T['facetKeys']:
    _lq = json.loads(re.search(r'SNP_LINE_QC\["zmgrin2026"\]\s*=\s*(\{.*\});', open(f'{root}/js/zmgrin.lineqc.js', encoding='utf-8').read(), re.S).group(1))['lines']
    for a in rows: a['sampleQC'] = {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'}[_lq[a['id']][1]]
    FACETS.append(['sampleQC', 'Sample heterozygosity'])
if [k for k, _ in FACETS] != T['facetKeys']: bad += 1; print('facet keys', T['facetKeys'])
for k, _ in FACETS:
    exp = {}
    for a in rows: exp[val(a, k)] = exp.get(val(a, k), 0) + 1
    if exp != T['facetCounts'][k]: bad += 1; print('facet differs:', k)
f = [a for a in rows if val(a, 'subpop') == 'SS' and val(a, 'inAmes282') == 'yes' and val(a, 'kernelType') == 'Dent']
if sorted(a['id'] for a in f) != T['filterSS_Ames_Dent']: bad += 1; print('compound filter differs')
# section counts while filtered: each section under every filter but its own; ticked values listed
F = {'subpop': {'SS'}, 'inAmes282': {'yes'}, 'kernelType': {'Dent'}}
exp_counts = {}
for k, _ in FACETS:
    c = {}
    for a in rows:
        if all(val(a, kk) in vs for kk, vs in F.items() if kk != k):
            c[val(a, k)] = c.get(val(a, k), 0) + 1
    for v in F.get(k, ()): c.setdefault(v, 0)
    exp_counts[k] = c
if T['facetsWhileFiltered'] != exp_counts:
    bad += 1; print('section counts while filtered differ:', [k for k in exp_counts if T['facetsWhileFiltered'].get(k) != exp_counts[k]])
if sorted(T['subpopListed']) != sorted(exp_counts['subpop']): bad += 1; print('subpopulations listed', T['subpopListed'])
orf = [a for a in rows if val(a, 'subpop') in ('SS', 'NSS') and val(a, 'inAmes282') == 'yes' and val(a, 'kernelType') == 'Dent']
if sorted(a['id'] for a in orf) != T['filterSSorNSS_Ames_Dent']: bad += 1; print('SS or NSS x Ames282 x Dent differs', len(orf), len(T['filterSSorNSS_Ames_Dent']))
print(f"while SS x Ames282 x Dent is ticked, Subpopulation still lists {len(exp_counts['subpop'])} values; SS or NSS -> {len(orf)} lines")
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
# the top bar's selection chip
C = json.load(open(res))['chip']
other = [a['id'] for a in rows if a.get('panel') == 'Other GRIN']
grown = len(set(nam) | set(other))
if C['onGeo'] != str(len(nam)): bad += 1; print('chip count on SNPGeo', C['onGeo'])
if not C['open']['drawer'] or C['open']['draft'] != len(nam) or 'SNPVersity' not in C['open']['title']: bad += 1; print('chip drawer', C['open'])
if C['applied']['n'] != grown or C['applied']['chip'] != str(grown) or C['applied']['tool'] != 'snpgeo' \
        or f"+{grown - len(nam)}" not in C['applied']['toast'] or 'Undo' not in C['applied']['toast']: bad += 1; print('chip apply', C['applied'])
if C['undone'] != {'n': len(nam), 'chip': str(len(nam))}: bad += 1; print('chip undo', C['undone'])
if C['onTrait'] != {'grids': 1, 'pageSetAside': True} or not C['traitBack']: bad += 1; print('chip on the SNPTrait page', C['onTrait'], C['traitBack'])
print(f"selection chip: {C['onGeo']} on SNPGeo; apply + Other GRIN -> {C['applied']['n']}; undo -> {C['undone']['n']}; SNPTrait page set aside and restored")
print(f'KW1000 in [250,300] g: {len(kw)} (x Ames282: {len(kwa)})')
# line QC (js/zmgrin.lineqc.js): the Sample heterozygosity facet, its exports, the hand-off of the
# heterozygous samples, the "het" marker; and the same pages without the file
LQ = json.load(open(res)).get('lineqc')
if LQ:
    lq_js = open(f'{root}/js/zmgrin.lineqc.js', encoding='utf-8').read()
    lq = json.loads(re.search(r'SNP_LINE_QC\["zmgrin2026"\]\s*=\s*(\{.*\});', lq_js, re.S).group(1))['lines']
    names = {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'}
    W, WO = LQ['with'], LQ['without']
    exp_counts = {}
    for a in rows: exp_counts[names[lq[a['id']][1]]] = exp_counts.get(names[lq[a['id']][1]], 0) + 1
    H = sorted(a['id'] for a in rows if lq[a['id']][1] == 'H')
    if W['counts'] != exp_counts or W['facetKeys'][-1] != 'sampleQC' or not W['facetBlock']: bad += 1; print('line QC facet', W['counts'], exp_counts)
    if W['selectedTrait'] != H or W['versity'] != {'tool': 'snpversity', 'selected': H}: bad += 1; print('heterozygous samples hand-off', W['versity'])
    hdr = W['exports']['csv'][0].split(',')
    csv_ok = 'sampleQC' in hdr and 'hetShare' in hdr and len(W['exports']['csv']) == len(H) + 1 and all(
        dict(zip(hdr, l.split(',')))['sampleQC'] == 'Heterozygous sample' and abs(float(dict(zip(hdr, l.split(',')))['hetShare']) - lq[l.split(',')[0]][0]) < 1e-9
        for l in W['exports']['csv'][1:])
    json_ok = len(W['exports']['json']) == len(H) and all(o['sampleQC'] == 'Heterozygous sample' and o['hetShare'] == lq[o['id']][0] for o in W['exports']['json'])
    if not (csv_ok and json_ok): bad += 1; print('line QC exports', hdr)
    ch9 = f"This sample is heterozygous at {lq['ZmG_CH9'][0] * 100:.1f}% of clean sites. Inbred lines are near 0.8%."
    if W['header']['marked'] != H or W['header']['n'] != len(H) or W['header']['ch9'] != {'text': 'het', 'tt': ch9}:
        bad += 1; print('header marker', W['header'])
    if sorted(c['run'] for c in W['chips'] if c['het']) != sorted(a['run'] for a in rows if a['id'] in set(H)) or not all(c['het'] for c in W['chips']):
        bad += 1; print('selected-line chips', W['chips'][:3])
    C = W['carriers']
    exp_marked = sorted(i for i in C['hom'] + C['het'] if lq[i][1] == 'H')
    if C['marked'] != exp_marked or C['shown'] != len(C['hom']) + len(C['het']) or not exp_marked: bad += 1; print('carrier chips', C['allele'], C['marked'], exp_marked)
    if WO['hasLineQc'] or 'sampleQC' in WO['facetKeys'] or WO['facetBlock'] or WO['header']['marked'] or any(c['het'] for c in WO['chips']) \
            or WO['carriers']['marked'] or WO['consoleErrors'] or W['consoleErrors'] or 'hetShare' in WO['exports']['csv'][0].split(','):
        bad += 1; print('without js/zmgrin.lineqc.js', {k: WO[k] for k in ('facetKeys', 'facetBlock', 'consoleErrors')}, W['consoleErrors'])
    print(f"line QC: {exp_counts}; Heterozygous sample -> SNPVersity selects {len(W['versity']['selected'])} lines; ZmG_CH9 marked het ({lq['ZmG_CH9'][0]:.4f}); "
          f"{len(C['marked'])} marked carriers of {C['allele']}; without the file: no facet, no marker, no error")
print(f'rows={len(rows)} facets={len(sch["facets"])} SSxAmes282xDent={len(f)} +iowa={len(fi)} mismatches={bad}')
sys.exit(1 if bad else 0)
