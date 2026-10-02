#!/usr/bin/env python3
"""check_snpgeo_counts.py <results.json> <fixture.vcf.gz> <site_root>

Independent recomputation of SNPGeo's per-site x per-country and per-state
statistics (total, called, carriers, het, altAlleles, calledAlleles) for the
gene-search scenario of run_scenarios.js, straight from the fixture VCF text,
js/snpgeo.regions.js and the dataset's sample list in js/zmgrin.catalog.js.
Exits non-zero on any mismatch. Standard library only."""
import gzip, json, re, sys
res, vcf, root = sys.argv[1:4]
R = json.load(open(res))['snpgeo_gene']
reg = open(f'{root}/js/snpgeo.regions.js', encoding='utf-8').read()
REG = json.loads(re.search(r'window\.SNPGEO_REGIONS\s*=\s*(\{.*?\});', reg, re.S).group(1))
cat = open(f'{root}/js/zmgrin.catalog.js', encoding='utf-8').read()
fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', cat, re.S).group(1))
ds = {a['id'] for p in fam['projects'] for g in p['groups'] for a in g['accessions']}
lo, hi = 12838008, 12843999            # Zm00001eb374090 as returned by lookupGeneModel.php
sites = {}
for line in gzip.open(vcf, 'rt'):
    if line.startswith('##'): continue
    t = line.rstrip('\n').split('\t')
    if line.startswith('#CHROM'): cols = t[9:]; continue
    p = int(t[1])
    if lo <= p <= hi: sites[p] = dict(zip(cols, t[9:]))
def exp_site(gts):
    out = {}
    for sid, geo in REG.items():
        if sid not in ds: continue
        iso = geo.get('iso3') or '???'
        c = out.setdefault(iso, dict(total=0, called=0, count=0, het=0, altAlleles=0, calledAlleles=0, states={}))
        st = c['states'].setdefault(geo.get('state') or geo.get('county') or 'Unknown', dict(total=0, called=0, count=0, het=0, altAlleles=0, calledAlleles=0))
        c['total'] += 1; st['total'] += 1
        g = gts.get(sid, './.').split(':')[0]
        al = re.split(r'[/|]', g)
        if '.' in al or '' in al: continue
        a = sum(x != '0' for x in al)
        for o in (c, st):
            o['called'] += 1; o['altAlleles'] += a; o['calledAlleles'] += len(al)
            if a: o['count'] += 1
            if 0 < a < len(al): o['het'] += 1
    return out
bad = 0; ncmp = 0
assert len(R['allStats']) == len(sites), (len(R['allStats']), len(sites))
for s in R['allStats']:
    e = exp_site(sites[s['pos']])
    if set(e) != set(s['c']): bad += 1; print('country set differs at', s['pos']); continue
    for iso, x in e.items():
        y = s['c'][iso]
        for k in ('total', 'called', 'count', 'het', 'altAlleles', 'calledAlleles'):
            ncmp += 1
            if x[k] != y[k]: bad += 1; print('mismatch', s['pos'], iso, k, x[k], y[k])
        for stn, xs in x['states'].items():
            ys = y['states'].get(stn)
            if ys is None: bad += 1; print('missing state', s['pos'], iso, stn); continue
            for k in ('total', 'called', 'count', 'het', 'altAlleles', 'calledAlleles'):
                ncmp += 1
                if xs[k] != ys[k]: bad += 1; print('state mismatch', s['pos'], iso, stn, k, xs[k], ys[k])
# the default view and the score columns (DNA: PlantCAD1, PlantCAD2, Evo2; protein: ESM1-3, ESM-C)
if R.get('defaultMode') != 'refalt' or not any(o and o[0] == 'refalt' for o in R.get('modeSelects', [])):
    bad += 1; print('default colour mode', R.get('defaultMode'), R.get('modeSelects'))
SCORES = ['PlantCAD1', 'PlantCAD2', 'Evo2', 'ESM1', 'ESM2', 'ESM3', 'ESM-C']
if [h for h in R['tableHead'] if h in SCORES] != SCORES: bad += 1; print('score columns', R['tableHead'])
print(f'sites={len(sites)} comparisons={ncmp} default view {R.get("defaultMode")}; score columns {" | ".join(SCORES)}; mismatches={bad}')
sys.exit(1 if bad else 0)
