#!/usr/bin/env python3
"""check_gwas_handoff.py <results.json> <site_root> -- check the GWAS Explorer -> SNPVersity
hand-off of run_scenarios.js against the catalogues, independently of the JS:
(a) MaizeGDB 2026 HQ/HC: every run of B73 + the 25 NAM founders in the B73-reference and NAM
    projects (B73 included); lines absent from the whole catalogue reported as not available;
(b) GRIN-linked set: the sample_ids whose manifest nam_founder is set; CML103 not available.
Standard library only."""
import csv, json, re, sys
res, root = sys.argv[1:3]
W = json.load(open(res))['gwas']
src = open(f'{root}/js/accessions.catalog.js', encoding='utf-8').read()
cat = json.loads(src[src.index('{'):src.rindex('}') + 1])
nam = [a for p in cat['families']['mgdb2026']['projects']
       if re.search(r'nested association mapping|Zm-B73-REFERENCE-NAM', p['title'], re.I)
       for g in p['groups'] for a in g['accessions']]
NAM = ['B73', 'B97', 'CML52', 'CML69', 'CML103', 'CML228', 'CML247', 'CML277', 'CML322', 'CML333', 'Hp301',
       'Il14H', 'Ki3', 'Ki11', 'Ky21', 'M37W', 'M162W', 'Mo18W', 'MS71', 'NC350', 'NC358', 'Oh7B', 'Oh43', 'P39',
       'Tx303', 'Tzi8']
norm = lambda v: re.sub(r'[^A-Z0-9]', '', str(v or '').upper())
allm = [a for p in cat['families']['mgdb2026']['projects'] for g in p['groups'] for a in g['accessions']]
exp_ids = sorted(a['id'] for a in nam if norm(a.get('founder')) in {norm(n) for n in NAM})
present = {norm(a.get('founder')) for a in allm} | {norm(a.get('namFounder')) for a in allm}
exp_missing = [n for n in NAM if norm(n) not in present]
S = list(csv.DictReader(open(f'{root}/data/zmgrin2026_samples.tsv'), delimiter='\t'))
grin = sorted(r['sample_id'] for r in S if r['nam_founder'])
bad = 0
if W['defaultSend']['dataset'] != 'mgdb2026_hq' or W['defaultSend']['ids'] != exp_ids: bad += 1; print('MaizeGDB 2026 ids differ')
b73 = [i for i in W['defaultSend']['ids'] if i.startswith('B73_')]
if not b73: bad += 1; print('B73 not sent')
for ds in ('mgdb2026_hq', 'mgdb2026_hc'):
    m = W['perSet'][ds]
    if f"{len(NAM) - len(exp_missing)} of {len(NAM)}" not in m['map'] or m['missing'] != 'Not available in this set: ' + ', '.join(exp_missing):
        bad += 1; print('popup text', ds, m)
g = W['perSet']['zmgrin2026_imp']
if '25 of 26' not in g['map'] or g['missing'] != 'Not available in this set: CML103': bad += 1; print('popup text GRIN', g)
if not W['replacePreticked'] or not all(v['replaceChecked'] for v in W['perSet'].values()): bad += 1; print('replace not pre-ticked')
if W['grinSend']['dataset'] != 'zmgrin2026_imp': bad += 1; print('dataset not switched')
if W['grinSend']['selected'] != grin: bad += 1; print('GRIN ids differ', set(W['grinSend']['selected']) ^ set(grin))
if 'CML103' not in W['grinMissing'] or W['grinMissing'].count(',') != 0: bad += 1; print('missing list:', W['grinMissing'])
if (W['grinSend']['chr'], W['grinSend']['start'], W['grinSend']['end']) != ('chr2', 4491424, 4499434): bad += 1; print('region differs')
print(f"MaizeGDB 2026 HQ/HC: {len(NAM) - len(exp_missing)} of {len(NAM)} lines, {len(exp_ids)} runs (B73: {len(b73)}), "
      f"not available: {', '.join(exp_missing)}; GRIN: {len(grin)} samples "
      f"({', '.join(grin[:3])}, ...); not available: {W['grinMissing'].split(':')[-1].strip()}; mismatches={bad}")
sys.exit(1 if bad else 0)
