#!/usr/bin/env python3
"""check_gwas_handoff.py <results.json> <site_root> -- check the GWAS Explorer -> SNPVersity
hand-off of run_scenarios.js: (a) default VCF set unchanged (every run of the MaizeGDB 2026 NAM
projects); (b) GRIN-linked set = the sample_ids whose manifest nam_founder is set, with CML103
reported as not available. Standard library only."""
import csv, json, re, sys
res, root = sys.argv[1:3]
W = json.load(open(res))['gwas']
src = open(f'{root}/js/accessions.catalog.js', encoding='utf-8').read()
cat = json.loads(src[src.index('{'):src.rindex('}') + 1])
nam = [a for p in cat['families']['mgdb2026']['projects']
       if re.search(r'nested association mapping|Zm-B73-REFERENCE-NAM', p['title'], re.I)
       for g in p['groups'] for a in g['accessions']]
S = list(csv.DictReader(open(f'{root}/data/zmgrin2026_samples.tsv'), delimiter='\t'))
grin = sorted(r['sample_id'] for r in S if r['nam_founder'])
bad = 0
if W['defaultSend']['dataset'] != 'mgdb2026_hq' or W['defaultSend']['selected'] != len(nam): bad += 1; print('default hand-off changed')
if W['grinSend']['dataset'] != 'zmgrin2026_imp': bad += 1; print('dataset not switched')
if W['grinSend']['selected'] != grin: bad += 1; print('GRIN ids differ', set(W['grinSend']['selected']) ^ set(grin))
if 'CML103' not in W['grinMissing'] or W['grinMissing'].count(',') != 0: bad += 1; print('missing list:', W['grinMissing'])
if (W['grinSend']['chr'], W['grinSend']['start'], W['grinSend']['end']) != ('chr2', 4491424, 4499434): bad += 1; print('region differs')
print(f"default: {W['defaultSend']['selected']} MaizeGDB 2026 NAM runs; GRIN: {len(grin)} samples "
      f"({', '.join(grin[:3])}, ...); not available: {W['grinMissing'].split(':')[-1].strip()}; mismatches={bad}")
sys.exit(1 if bad else 0)
