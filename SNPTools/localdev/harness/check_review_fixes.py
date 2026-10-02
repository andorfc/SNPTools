#!/usr/bin/env python3
"""check_review_fixes.py <results.json> <siteRoot> -- the code-review fixes of patch 0038.

From the jsdom run (results.json -> fixes):
  SNPTrait "Random %" with no line visible does nothing (it threw); randomSample() draws k
  distinct ids, flat over catalogue order (the sort-based shuffle picked the first lines ~5x
  as often as the last); a SNPVersity result keeps the region it was queried for after the
  form moves on (CHR column, JBrowse link, heading, hand-off); the newer of two Builds, two
  SNPFunction genes, two SNPGeo searches and two SNPFold genes stands when the older answers
  last; a SNPGeo search answered after the user left does not draw over SNPVersity; SNPFold
  makes one variant query and one gene lookup per gene.

On the server, in a temporary site root (so the real vcf/ folder is not pruned):
  processForm.php refuses a request over SNPTOOLS_MAX_CELLS with a message and writes no
  file, keeps string ids once each and fills a reserved name (INFO) with ./., never reuses an
  existing file name, ignores a path in outName, refuses a non-integer interval, and removes
  its own VCFs older than the TTL; ibsCompare.php rebuilds an unreadable offset index.
"""
import gzip, json, math, os, shutil, subprocess, sys, tempfile, time

res = json.load(open(sys.argv[1]))
ROOT = os.path.abspath(sys.argv[2])
HERE = os.path.dirname(os.path.abspath(__file__))
PHP = os.environ.get('PHP_BIN', 'php')
bad = []


def check(ok, what):
    print(('ok   ' if ok else 'FAIL ') + what)
    if not ok:
        bad.append(what)


F = res['fixes']
e = F['traitRandomEmpty']
check(e['visible'] == 0 and e['threw'] is None and e['selected'] == ['ZmG_B73'],
      f"SNPTrait Random 5% with nothing visible: no error, selection kept ({e})")
t = F['traitRandom10']
check(t['selected'] == max(1, round(t['visible'] * .10)), f"SNPTrait Random 10%: {t['selected']} of {t['visible']}")
s = F['sample']
flat = abs(s['first100'] / s['expected'] - 1) < .10 and abs(s['last100'] / s['expected'] - 1) < .10
check(s['distinct'] and flat and s['zero'] == 0 and s['over'] == 3,
      f"randomSample: {s['k']} distinct of {s['N']}; first 100 lines {s['first100']:.1f}, last 100 {s['last100']:.1f} picks "
      f"(expected {s['expected']:.1f} over {s['T']} draws)")
q = F['quickPick5']
check(q['selected'] == math.ceil(q['of'] * .05), f"SNPVersity quick pick Random 5%: {q['selected']} of {q['of']}")

qd = F['queried']
check([qd['chr'], qd['lo'], qd['hi']] == ['chr9', 13118306, 13124164], f"result records its query {qd}")
a = F['afterFormChange']
check(a['form'][0] == 'chr2' and a['chrCell'] == '9' and 'loc=chr9:' in a['jbrowse'] and 'highlight=chr9:' in a['jbrowse']
      and 'chr9:13,118,306' in a['heading'],
      f"form moved to {a['form']}: table still reads chr9 (CHR {a['chrCell']}, {a['jbrowse'][:60]}..., '{a['heading']}')")
ti = F['treeInput']
check([ti['chr'], ti['start'], ti['end']] == ['chr9', 13118306, 13124164] and ti['rows'] > 0,
      f"Send to SNPTree hands off the queried region {ti}")
check('chr9:13,118,306' in F['backHeading'], f"back in SNPVersity: '{F['backHeading']}'")

r = F['runRace']
check(r['calls'] == 2 and r['resultChr'] == 'chr10', f"two Builds, the first answering last: result is the second's ({r})")
fr = F['functionRace']
check(fr['showsB'] and not fr['showsA'], f"SNPFunction: gene B's answer stands over A's late one ({fr})")
g = F['geoLeft']
check(g['tool'] == 'snpversity' and not g['geoMap'] and g['runbar'] and not g['geoInput'],
      f"SNPGeo search answered after leaving SNPGeo: SNPVersity untouched ({g})")
g2 = F['geoTwo']
check(g2['gene'] == 'Zm00001eb374090' and g2['rows'], f"SNPGeo: the newer search's gene stands ({g2})")
fq = F['foldRequests']
check(fq['processForm'] == 1 and fq['lookupGeneModel'] == 1 and fq['shown'] == 'Zm00001eb406050',
      f"SNPFold loads a gene with one variant query and one gene lookup ({fq})")
fo = F['foldRace']
check(fo['calls'] == ['Zm00001eb378140', 'Zm00001eb406050'] and fo['shown'] == 'Zm00001eb406050',
      f"SNPFold: both genes waiting on variants, A answering last: B stands ({fo})")
errs = [x for x in F.get('consoleErrors', []) if 'TypeError' in x or 'is not a function' in x or 'ReferenceError' in x]
check(not errs, f"no script errors in the scenarios ({errs[:2]})")


# ---------------- server side, in a temporary root ----------------
def php(root, script, method='GET', get=None, post=None, env=None):
    req = {'method': method, 'get': get or {}, 'post': post or {}, 'root': root, 'script': script}
    out = subprocess.run([PHP, os.path.join(HERE, 'php_shim.php')], capture_output=True, text=True,
                         env=dict(os.environ, SNPT_REQ=json.dumps(req), **(env or {})))
    try:
        return json.loads(out.stdout)
    except ValueError:
        return {'raw': out.stdout[:300], 'stderr': out.stderr[:300]}


tmp = tempfile.mkdtemp(prefix='snpt_fixes_')
try:
    for f in ('processForm.php', 'h5_to_vcf.py', 'ibsCompare.php'):
        shutil.copy(os.path.join(ROOT, f), tmp)
    os.makedirs(os.path.join(tmp, 'hdf5'))
    os.symlink(os.path.join(ROOT, 'hdf5', 'version3'), os.path.join(tmp, 'hdf5', 'version3'))
    vdir = os.path.join(tmp, 'vcf'); os.makedirs(vdir)
    old = time.time() - 3 * 86400
    for name in ('snpv_1_old_1_2.vcf.gz', 'snpv_1_old_1_2.vcf.gz.acc.json', 'notes.txt'):
        open(os.path.join(vdir, name), 'w').write('x'); os.utime(os.path.join(vdir, name), (old, old))
    open(os.path.join(vdir, 'snpv_2_recent_1_2.vcf.gz'), 'w').write('x')

    win = {'chr': 'chr9', 'start': '13118306', 'end': '13124164', 'dataSet': 'zmgrin2026_imp'}
    ids = json.dumps(['ZmG_B73', 'ZmG_B73', 5, 'INFO', 'ZmG_CML103'])
    j = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=ids, outName='vcf/snpv_7_t_13118306_13124164.vcf.gz'))
    left = sorted(os.listdir(vdir))
    check(j.get('status') == 'success' and 'snpv_1_old_1_2.vcf.gz' not in left and 'snpv_1_old_1_2.vcf.gz.acc.json' not in left
          and 'notes.txt' in left and 'snpv_2_recent_1_2.vcf.gz' in left,
          f"old VCFs pruned, other files kept: {left}")
    first = os.path.join(vdir, 'snpv_7_t_13118306_13124164.vcf.gz')
    with gzip.open(first, 'rt') as fh:
        lines = [l.rstrip('\n').split('\t') for l in fh if not l.startswith('##')]
    head, body = lines[0], lines[1:]
    check(head[9:] == ['ZmG_B73', 'INFO', 'ZmG_CML103'] and {r[10] for r in body} == {'./.'} and len(body) == j.get('variants'),
          f"ids once each, numbers dropped, INFO filled with ./.: {head[9:]}, {len(body)} rows")
    before = open(first, 'rb').read()
    j2 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73']), outName='vcf/snpv_7_t_13118306_13124164.vcf.gz'))
    check(j2.get('status') == 'success' and os.path.basename(j2.get('outFile', '')) != 'snpv_7_t_13118306_13124164.vcf.gz'
          and open(first, 'rb').read() == before, f"an existing name is not reused: {j2.get('outFile')}")
    j3 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73']), outName='../index.html'))
    check(j3.get('status') == 'success' and j3.get('outFile', '').startswith('./vcf/snpv_') and not os.path.exists(os.path.join(tmp, 'index.html')),
          f"a path in outName is ignored: {j3.get('outFile')}")
    j4 = php(tmp, 'processForm.php', 'POST', post=dict(win, start='1e5', genotypes='[]', outName=''))
    check(j4.get('status') == 'error' and 'whole numbers' in j4.get('message', ''), f"start=1e5 refused: {j4.get('message')}")
    n_before = len(os.listdir(vdir))
    j5 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73', 'ZmG_B97']), outName=''),
             env={'SNPTOOLS_MAX_CELLS': '100'})
    check(j5.get('status') == 'error' and j5.get('tooLarge') and 'too large' in j5.get('message', '')
          and len(os.listdir(vdir)) == n_before, f"over SNPTOOLS_MAX_CELLS: refused, nothing written: {j5.get('message', '')[:110]}")

    # ibsCompare.php: an unreadable offset index is rebuilt, written whole
    dist = os.path.join(tmp, 'distance', 'zz'); os.makedirs(dist)
    open(os.path.join(dist, 'ids.txt'), 'w').write('a\nb\nc\n')
    open(os.path.join(dist, 'similarity.csv'), 'w').write('1,0.5,0.25\n0.5,1,0.75\n0.25,0.75,1\n')
    open(os.path.join(dist, 'missing_pct.csv'), 'w').write('0,0.1,0.2\n0.1,0,0.3\n0.2,0.3,0\n')
    for f in ('similarity.csv', 'missing_pct.csv'):
        idx = os.path.join(dist, f + '.offidx'); open(idx, 'w').write('a:3:{i:0;i:0;i:1')   # a torn write
        t2 = time.time() + 5; os.utime(idx, (t2, t2))
    jb = php(tmp, 'ibsCompare.php', get={'focal': 'b', 'dataset': 'zz'}, env={'SNPTOOLS_DISTANCE_DIR': os.path.join(tmp, 'distance')})
    rows = jb.get('rows') or []
    leftovers = [f for f in os.listdir(dist) if f.endswith('.tmp')]
    check([(x['id'], x['similarity'], round(x['missing'], 6)) for x in rows] == [('a', 0.5, 10.0), ('b', 1.0, 0.0), ('c', 0.75, 30.0)]
          and not leftovers, f"ibsCompare.php with a torn .offidx: {rows}, temp files left {leftovers}")
    jb2 = php(tmp, 'ibsCompare.php', get={'focal': 'c', 'dataset': 'zz'}, env={'SNPTOOLS_DISTANCE_DIR': os.path.join(tmp, 'distance')})
    check([x['similarity'] for x in jb2.get('rows') or []] == [0.25, 0.75, 1.0], "the rebuilt index reads back")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"review fixes: {len(bad)} failed; mismatches={len(bad)}")
sys.exit(1 if bad else 0)
