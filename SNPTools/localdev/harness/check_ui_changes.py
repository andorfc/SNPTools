#!/usr/bin/env python3
"""check_ui_changes.py <results.json> <site_root> -- the change list of 3 October 2026.

run_scenarios.js (uiList) drives the pages; recomputed here from the VCF, the catalogue and the
SNPCurate data file:
(a) the name SNPMaize and the tagline in the left panel, the masthead, the breadcrumb and the title;
(b) SNPVersity (tga1 interval, 26 NAM lines): maxR² shown to hundredths (from the VCF's MAXR2), the
    PlantCAD1 / PlantCAD2 headers broken in two lines (the other score headers not), the curated row
    on gold (its site and annotation cells, not the score or genotype cells), and each gene model of
    "Gene models in this region" linked to its MaizeGDB page;
(c) SNPTree, SNPImpact, SNPCompare: no "Load data from SNPVersity" before SNPVersity has a result;
    afterwards the offer names its region, lines and variants, and the button loads it;
(d) SNPCompare's PI number column equals each line's GRIN accession in the catalogue;
(e) SNPFunction's "Reference lines" hand-off sends exactly the lines homozygous for the reference
    allele (the panel minus the carriers and the missing calls).
(Widths at 375 px are checked in a browser; jsdom does no layout.) Standard library only."""
import gzip, json, os, re, sys

res, root = sys.argv[1:3]
L = json.load(open(res)).get('uiList')
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond: fails.append(msg)


if not L:
    print('change list of 3 October: skipped (the full chr2 / chr4 stores are not installed)'); sys.exit(0)
TAG = 'A SNP toolkit to explore variant diversity across maize germplasm'
b = L['brand']
check(b == {'title': 'SNPMaize — SNPVersity 2.1', 'rail': 'SNPMaize', 'railTag': TAG, 'mast': 'SNPMaize', 'mastTag': TAG, 'crumb': 'SNPMaize'},
      f"name and tagline: {b['rail']} / {b['mast']} / crumb {b['crumb']} / title {b['title']}")

# (b) SNPVersity
v = L['versity']
p = os.path.join(root, v['vcf'][2:] if v['vcf'].startswith('./') else v['vcf'])
r2 = {}
for line in gzip.open(p, 'rt'):
    if line.startswith('#'): continue
    t = line.split('\t', 8)
    I = dict(kv.split('=', 1) for kv in t[7].split(';') if '=' in kv)
    r2[f"{int(t[1]):,}"] = I.get('MAXR2')
want = [[pos, 'NA' if r2.get(pos) in (None, '', '.') else f"{float(r2[pos]):.2f}"] for pos, _ in v['r2']]
check(v['r2'] == want and len(v['r2']) == len(r2), f"maxR² to hundredths on all {len(v['r2'])} rows (e.g. {v['r2'][1]})")
check(v['pcHead'] == ['Plant<br>CAD1', 'Plant<br>CAD2'] and all('<br>' not in h for h in v['otherHead']),
      f"PlantCAD headers in two lines {v['pcHead']}, the others in one {v['otherHead']}")
cur = json.loads(re.search(r'window\.SNP_CURATE = (\{.*\});', open(f'{root}/js/snpcurate.data.js', encoding='utf-8').read(), re.S).group(1))
want_cur = [e for e in cur['entries'] if e['status'] == 'site' and e['chr'] == 'chr4' and 46647932 <= e['pos'] <= 46652896]
cr = v['curRows']
check(len(cr) == len(want_cur) == 1 and cr[0]['id'] == want_cur[0]['id'] and cr[0]['cls'] == f"cur-row cur-row-{want_cur[0]['mark']}"
      and cr[0]['pos'] == f"{want_cur[0]['pos']:,}" and cr[0]['score'] == 7 and cr[0]['gt'] == 26 and cr[0]['plain'] == 4 + 7,
      f"curated row on gold: {[(c['id'], c['pos'], c['cls']) for c in cr]}, {cr[0]['plain'] if cr else 0} site and annotation cells, scores and genotypes not")
g = v['genes']
check(g and all(x['href'] == f"https://www.maizegdb.org/gene_center/gene/{x['gene']}" and x['maizegdb'] == [x['href']] for x in g),
      f"gene models in this region linked to MaizeGDB: {[x['gene'] for x in g]}")

# (c) Load data from SNPVersity
check(not any(L['emptyBefore'].values()), f"no offer before SNPVersity has a result {L['emptyBefore']}")
for t, x in L['load'].items():
    check(x['offer'] and 'chr4:46,647,932–46,652,896 · 26 lines' in x['offer'] and x['tool'] == t and x['loaded'] and not x['block'],
          f"{t}: offers \"{(x['offer'] or '')[:60]}...\" and the button loads it")

# (d) SNPCompare PI number
cat_js = open(f'{root}/js/zmgrin.catalog.js', encoding='utf-8').read()
fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', cat_js, re.S).group(1))
grin = {a['id']: a.get('grin') or '' for p_ in fam['projects'] for g_ in p_['groups'] for a in g_['accessions']}
c = L['compare']
check('PI number' in c['head'] and len(c['rows']) == len(grin) and all(pi == grin[i] for i, pi in c['rows']),
      f"SNPCompare PI number column: {len(c['rows'])} rows equal the catalogue's GRIN accession (e.g. {c['rows'][0]})")

# (e) reference lines
rf = L['reference']
want_ref = sorted(set(grin) - set(rf['het']) - set(rf['hom']) - set(rf['miss']))
check(rf['tool'] == 'snpversity' and rf['selected'] == want_ref and rf['nRef'] == len(want_ref) and f"Reference lines ({len(want_ref)}) →" in rf['buttons'],
      f"SNPFunction ZmWAK T102A: Reference lines ({len(want_ref)}) sends exactly the lines homozygous for the reference allele")
check(not L.get('consoleErrors'), f"no script errors {L.get('consoleErrors')}")
print(f"change list of 3 October: {len(fails)} failed; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
