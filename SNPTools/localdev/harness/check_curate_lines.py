#!/usr/bin/env python3
"""check_curate_lines.py -- SNPCurate's Your lines panel (the allele passport): the curateLines scenario of
run_scenarios.js against js/snpcurate.data.js.

    python3 check_curate_lines.py <results.json> <root>

Exit 1 on a failure."""
import json, re, sys

res, root = sys.argv[1:3]
L = json.load(open(res)).get('curateLines')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not L:
    print('skip no curateLines scenario (the full stores are not installed)'); sys.exit(0)
js = open(f'{root}/js/snpcurate.data.js', encoding='utf-8').read()
D = json.loads(js[js.index('{'):js.rstrip().rstrip(';').rindex('}') + 1])
MO = {'gold': 0, 'silver': 1, 'white': 2}
sites = sorted([e for e in D['entries'] if e.get('site')], key=lambda e: (MO[e['mark']], e['symbol'].lower(), e['id']))
def gt(e, i):
    c = e['site']['carriers']
    return 'alt' if i in c['hom'] else 'het' if i in c['het'] else 'missing' if i in c.get('missing', []) else 'ref'
GLYPH = {'alt': 'H', 'het': 'h', 'ref': '·', 'missing': '?'}
def fav(e, g):
    f = e.get('favorable')
    if f == 'alt': return 'hom' if g == 'alt' else 'het' if g == 'het' else ''
    if f == 'ref': return 'hom' if g == 'ref' else 'het' if g == 'het' else ''
    return ''

st = L['state']
check(st['ids'] == ['ZmG_B73', 'ZmG_IL14H', 'ZmG_P39'] and st['report']['unmatched'] == ['notaline'] and st['report']['matched'] == 3,
      f"pasted 'B73', 'Il14H', 'PI 550473' (B73 again), 'P39', 'notaline' resolve to {st['ids']}, not found {st['report']['unmatched']}")
cols = [e['id'] for e in sites]
check(L['cols'] == cols, f"{len(L['cols'])} columns: the genotyped alleles in the table's order (gold, silver, white; then gene)")
ok = True
for row in L['rows']:
    i = row['id']
    want = [GLYPH[gt(e, i)] for e in sites]
    want_fav = sum(fav(e, gt(e, i)) == 'hom' for e in sites); want_alt = sum(gt(e, i) == 'alt' for e in sites)
    favn = sum(e.get('favorable') in ('alt', 'ref') for e in sites)
    if row['cells'] != want or row['fav'] != f"{want_fav} of {favn}" or row['alt'] != str(want_alt):
        ok = False; print('   row', i, row['cells'], want, row['fav'], want_fav, row['alt'], want_alt)
    if i == 'ZmG_IL14H':
        check(gt(next(e for e in sites if e['id'] == 'SC0007'), i) == 'alt' and gt(next(e for e in sites if e['id'] == 'SC0003'), i) == 'alt',
              'Il14H is homozygous alternate at su1 F163L (SC0007) and the DGAT1-2 F469 insertion (SC0003) in the data')
check(ok, f"every cell, favorable count and alternate count of the {len(L['rows'])} rows match the registry's carrier lists")
b73 = next(r for r in L['rows'] if r['id'] == 'ZmG_B73')
check(all(c == '·' for c in b73['cells']) and b73['alt'] == '0', f"B73 is reference at every site ({b73['alt']} alternate)")
il = next(r for r in L['rows'] if r['id'] == 'ZmG_IL14H')
check(il['favClasses']['SC0003'] == 'hom' and il['favClasses']['SC0007'] == '' and L['tsv'].splitlines()[0].startswith('accession\tline\tsubpopulation\tSC0003 DGAT1-2'),
      f"the F469 cell of Il14H is coloured as the favorable allele, the su1 cell (depends) is not; TSV header {L['tsv'].splitlines()[0][:60]!r}")
tsv = [l.split('\t') for l in L['tsv'].strip().split('\n')]
check(len(tsv) == 4 and tsv[2][0] == 'ZmG_IL14H' and tsv[2][3:3 + len(sites)] == [gt(e, 'ZmG_IL14H') for e in sites] and tsv[2][-2:] == [il['fav'].split(' ')[0], il['alt']],
      f"TSV: header + 3 rows, Il14H row carries the genotype words and the counts")
check(L['hash'] == '#snpcurate?lines=ZmG_B73,ZmG_IL14H,ZmG_P39', f"the hash carries the lines: {L['hash']}")
sel = L['selectCarriers']
e3 = next(e for e in sites if e['id'] == 'SC0003')
check(sorted(sel['selected']) == sorted(e3['site']['carriers']['hom'] + e3['site']['carriers']['het']) and sel['n'] == e3['site']['nHet'] + e3['site']['nHom'],
      f"clicking the DGAT1-2 F469 column selects its {sel['n']} carriers for SNPVersity")
check(L['route']['ids'] == ['ZmG_B73', 'ZmG_IL14H', 'ZmG_P39'] and L['route']['rows'] == 3, f"a pasted #snpcurate?lines= link rebuilds the panel: {L['route']}")
check(L['founders'] == 26, f"NAM founders: {L['founders']} lines")
check(not L['errors'], f"no console errors ({len(L['errors'])})")
print(f'check_curate_lines: {bad} failure(s)')
sys.exit(1 if bad else 0)
