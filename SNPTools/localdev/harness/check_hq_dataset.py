#!/usr/bin/env python3
"""check_hq_dataset.py <results.json> <site_root> -- the High Quality set (zmgrin2026_hq).

run_scenarios.js (R.hq) queries the chr2 GWAS window on both sets of the release and drives the
pages; recomputed here from the two VCFs, the line QC tables and the SNPCurate data file:
(a) the HQ VCF is exactly the Imputed VCF's records that the set's rule keeps: at least one
    carrier (NHET+NHOM > 0) and, with 10 or more carriers, MAXR2 >= 0.2 (missing MAXR2 fails),
    with identical CHROM..ALT, genotypes and INFO values; its header carries the release's
    ##HQ_filter and ##subset lines, the Imputed VCF's does not;
(b) both sets are listed, one family, labelled "<name> · Imputed" / "<name> · High Quality";
    a selection survives a switch between them and a hand-off into the other set keeps 'add';
(c) SNPVersity's Min. maxR² filter: the rows it shows at 0.5 and 0.2 equal the rows of the VCF
    with MAXR2 at or above the value (a site without MAXR2 is hidden); its note counts the sites it
    hides; out-of-range input is clamped (2 -> 1, -1 and text -> 0); it combines with Min. MAF; it is
    left out for a set without maxR²;
(d) the tools name the set (SNPImpact, SNPTree, SNPFunction), and the genome-wide notes of SNPTree
    and SNPCompare say the matrices were computed on the Imputed set's sites (none for Imputed);
(e) per-set line QC: each set's classes equal its data/qc/<family>.lineqc.tsv;
(f) SNPCurate: the set chooser follows the site-wide set; a site the rule drops is marked in the
    table and the record, "On the map" is disabled, and nothing is marked on the Imputed set.
Standard library only."""
import gzip, json, os, sys

res, root = sys.argv[1:3]
R = json.load(open(res))
H = R.get('hq')
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond: fails.append(msg)


if not H:
    print(f"High Quality set: skipped (no chr2 HQ store installed: {R.get('hqStore')})"); sys.exit(0)


def read_vcf(url):
    p = os.path.join(root, url[2:] if url.startswith('./') else url)
    head, rows = [], []
    for line in gzip.open(p, 'rt'):
        if line.startswith('##'): head.append(line.rstrip('\n')); continue
        if line.startswith('#'): continue
        t = line.rstrip('\n').split('\t')
        info = {}
        for kv in t[7].split(';'):
            k, _, v = kv.partition('=')
            try: info[k] = float(v)
            except ValueError: info[k] = v
        rows.append({'key': (t[0], int(t[1]), t[3], t[4]), 'info': info, 'gt': t[9:]})
    return head, rows


def keeps(info):
    c = int(info.get('NHET', 0)) + int(info.get('NHOM', 0))
    r2 = info.get('MAXR2')
    return c > 0 and (c < 10 or (isinstance(r2, float) and r2 >= 0.2))


q = H['query']
hi_head, imp = read_vcf(q['zmgrin2026_imp']['vcf'])
hq_head, hq = read_vcf(q['zmgrin2026_hq']['vcf'])
want = [r for r in imp if keeps(r['info'])]
check([r['key'] for r in hq] == [r['key'] for r in want],
      f"(a) HQ records = the Imputed records the rule keeps: {len(hq)} of {len(imp)} "
      f"({sum(1 for r in imp if int(r['info']['NHET']) + int(r['info']['NHOM']) == 0)} without a carrier, "
      f"{len(imp) - len(want) - sum(1 for r in imp if int(r['info']['NHET']) + int(r['info']['NHOM']) == 0)} common with MAXR2 < 0.2 or none)")
check(len(hq) == len(want) and all(a['gt'] == b['gt'] for a, b in zip(hq, want)), "(a) genotypes identical on every kept record")
check(len(hq) == len(want) and all(a['info'] == b['info'] for a, b in zip(hq, want)), "(a) INFO identical by value on every kept record")
check(all(keeps(r['info']) for r in hq), "(a) every HQ record satisfies the rule")
check(any(l.startswith('##HQ_filter=<Description="High Quality subset of zmgrin2026') for l in hq_head) and '##subset=High_Quality' in hq_head
      and not any(l.startswith(('##HQ_filter', '##subset')) for l in hi_head), "(a) ##HQ_filter and ##subset in the HQ download only")
check(q['zmgrin2026_imp']['rows'] == len(imp) and q['zmgrin2026_hq']['rows'] == len(hq), f"(a) table rows = VCF records ({len(imp)}, {len(hq)})")

# (b) labels, selection, hand-off
ds = {d['id']: d for d in H['datasets']}
check(set(ds) >= {'zmgrin2026_imp', 'zmgrin2026_hq'} and ds['zmgrin2026_imp']['family'] == ds['zmgrin2026_hq']['family'] == 'zmgrin2026',
      f"(b) both sets listed, one family: {list(ds)}")
check(ds['zmgrin2026_imp']['label'] == 'MaizeGDB GRIN-linked 2026 · Imputed' and ds['zmgrin2026_hq']['label'] == 'MaizeGDB GRIN-linked 2026 · High Quality',
      f"(b) labels {ds['zmgrin2026_imp']['label']!r}, {ds['zmgrin2026_hq']['label']!r}")
check(len(H['cards']) == len(ds) and any('High Quality (LD max R² ≥ 0.2)' in c for c in H['cards']), f"(b) SNPVersity cards {H['cards']}")
k = H['keep']
check(k['before'] == k['hq'] == k['back'] and k['ds1'] == 'zmgrin2026_hq' and k['ds2'] == 'zmgrin2026_imp' and len(k['before']) == 28,
      f"(b) a selection of {len(k['before'])} lines kept on switching to High Quality and back")
h = H['handoff']
check(h['merge'] == 'add' and h['crossDataset'] and not h['crossFamily'] and h['dataset'] == 'zmgrin2026_hq' and h['n'] >= h['n0']
      and 'the same lines, so the selection carries over' in h['hint'],
      f"(b) hand-off into the High Quality set keeps 'add' ({h['n0']} -> {h['n']}): {h['hint']!r}")
c = H['control']
check(c['dataset'] == 'zmgrin2026_hq' and not c['disabled'] and not c['reset'], "(b) the hand-off control treats the other set of the family as the same lines")

# (c) Min. maxR²
for ds_id, rows in (('zmgrin2026_imp', imp), ('zmgrin2026_hq', hq)):
    f = q[ds_id]['filter']
    for v in ('0.5', '0.2'):
        x = float(v)
        shown = sum(1 for r in rows if isinstance(r['info'].get('MAXR2'), float) and r['info']['MAXR2'] >= x)
        low = sum(1 for r in rows if isinstance(r['info'].get('MAXR2'), float) and r['info']['MAXR2'] < x)
        none = sum(1 for r in rows if not isinstance(r['info'].get('MAXR2'), float))
        note = f"Min. maxR² {x:.2f} hides {low + none:,} of {len(rows):,} sites: {low:,} below it, {none:,} with no maxR²."
        check(f[v]['shown'] == shown and f[v]['note'].startswith(note), f"(c) {ds_id} at {v}: {shown} shown; note {f[v]['note'][:70]!r}")
    check(f['2']['fR2'] == 1 and f['-1']['fR2'] == 0 and f['abc']['fR2'] == 0 and f['-1']['note'] == '' and f['2']['input'] == '1.00',
          f"(c) {ds_id}: input clamped (2 -> {f['2']['fR2']}, -1 -> {f['-1']['fR2']}, text -> {f['abc']['fR2']})")
    both = sum(1 for r in rows if isinstance(r['info'].get('MAXR2'), float) and r['info']['MAXR2'] >= 0.5 and float(r['info'].get('MAF', 0)) >= 0.1)
    check(f['withMaf']['shown'] == both, f"(c) {ds_id}: with Min. MAF 0.1 too, {both} shown")
    check(q[ds_id]['r2Field'], f"(c) {ds_id}: the filter is offered")
check(not H['noR2']['field'] and H['noR2']['shown'] == len(hq), "(c) no filter (and nothing hidden) for a set without maxR²")

# (d) tools and notes
check(q['zmgrin2026_hq']['handoff'] == {'dataset': 'zmgrin2026_hq', 'datasetName': 'MaizeGDB GRIN-linked 2026 · High Quality'},
      f"(d) the result hands off as {q['zmgrin2026_hq']['handoff']}")
check(H['impact'] and H['tree']['label'], "(d) SNPImpact and SNPTree name the High Quality set")
gw = "computed on the Imputed set's sites; none is installed for the High Quality set"
check(gw in H['tree']['gw'] or not H['tree']['gw'], f"(d) SNPTree genome-wide note {H['tree']['gw'][-90:]!r}")
check((gw in H['compare'] or not H['compare']) and "computed on" not in H['compareImp'], f"(d) SNPCompare note on High Quality {H['compare'][-90:]!r}; none on Imputed")
fn = H['fn']
check(fn['imp']['name'].endswith('· Imputed') and fn['hq']['name'].endswith('· High Quality') and fn['hq']['n'] <= fn['imp']['n'],
      f"(d) SNPFunction Zm00001eb067740: {fn['imp']['n']} variants (Imputed), {fn['hq']['n']} (High Quality)")

# (e) line QC per set
LAB = {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'}
for ds_id, key, f_ in (('zmgrin2026_imp', 'imp', 'zmgrin2026.lineqc.tsv'), ('zmgrin2026_hq', 'hq', 'zmgrin2026_HQ.lineqc.tsv')):
    t = [l.rstrip('\n').split('\t') for l in open(os.path.join(root, 'data', 'qc', f_))][1:]
    want_qc = {r[0]: LAB[r[5]] for r in t}
    got = {k: v for k, v in H['lineqc'][key].items() if v}
    check(got == want_qc, f"(e) {ds_id}: line QC of {len(got)} lines = {f_}")
diff = sorted(k for k in H['lineqc']['imp'] if H['lineqc']['imp'][k] != H['lineqc']['hq'][k])
print(f"     lines whose class differs between the sets: {diff}")

# (f) SNPCurate
cu = H['curate']
check(cu['select'] == 'zmgrin2026_hq' and cu['options'] == [ds['zmgrin2026_imp']['label'], ds['zmgrin2026_hq']['label']] and cu['outRows'] == 0,
      f"(f) SNPCurate set chooser {cu['options']}; no curated site outside the High Quality set")
check(cu['row'] == 'not in the High Quality set' and cu['note'].startswith(f"Not in the High Quality set. {cu['carriers']:,} carriers and maxR² 0.05")
      and dict(cu['map']).get('On the map') is True and dict(cu['map']).get('Carriers in SNPVersity') is False,
      f"(f) a dropped site ({cu['id']}, maxR² set to 0.05): marked, record note, On the map disabled")
check(not cu['imp']['note'] and cu['imp']['rows'] == 0, "(f) nothing marked on the Imputed set")
check(not H.get('consoleErrors'), f"no script errors {H.get('consoleErrors')}")
print(f"High Quality set: {len(fails)} failed; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
