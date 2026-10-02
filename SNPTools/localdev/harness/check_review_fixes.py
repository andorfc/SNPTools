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

Gene lookups (patch 0039): every gene of gff/genes_data.serialized reads back identically from
gff/genes_index.txt, and misses stay misses; in a temporary root, lookupGeneModel.php answers
from a current index without touching it, rebuilds a stale, missing or truncated one, and reads
the store itself when the index is stale and gff/ is read-only.

Patch 0040: no tooltip listeners are added to [data-tt] elements however often pages are shown,
and the delegated tooltip shows the innermost element's text, hides off it and when its element
is re-rendered away; SNPFold's pLDDT legend names the model shown (AlphaFold2, ESMFold);
PlantCAD1/2, Evo2 and ESM-C are available columns of the GRIN-linked set; SNPFunction's default
gene has variants; processForm.php and ibsCompare.php replies name no server path (the details
come back only with SNPTOOLS_DEBUG=1), and ibsCompare.php without a dataset serves zmgrin2026.
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
check(F['ttListeners'] == 0, f"tooltip listeners added to [data-tt] elements over 3 visits of 3 tools: {F['ttListeners']}")
tb = F['ttBehaviour']
check(tb == ['outer A', 'inner B', False, False], f"delegated tooltip: outer, inner, off, element removed -> {tb}")
fl = F['foldLegend']
check(fl.get('alphafold', {}).get('legend') == 'AlphaFold2 confidence score (pLDDT):'
      and fl.get('esmfold', {}).get('legend') == 'ESMFold confidence score (pLDDT):',
      f"SNPFold pLDDT legend names the model shown: {fl}")
fs = F['fieldStatus']
check(all(fs.get(k) == 'ok' for k in ('pc1', 'pc2', 'evo2', 'esmc')) and fs.get('mq') == 'hidden' and fs.get('comp') == 'hidden',
      f"GRIN-linked score columns available: {fs}")
fd = F['functionDefault']
check(fd['gene'] == 'Zm00001eb406050' and fd['n'] > 0 and not fd['error'],
      f"SNPFunction default gene has variants: {fd}")
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

    # ---- gene lookups: gff/genes_index.txt (gene_index_lib.php) ----
    probe = subprocess.run([PHP, '-r', r"""
        require $argv[1] . '/gene_index_lib.php';
        list($s, $i) = gene_index_paths($argv[1] . '/gff');
        $d = unserialize(file_get_contents($s)); $ix = gene_index_open($s, $i); $bad = 0;
        if (!$ix) { echo json_encode(['open' => false]); exit; }
        foreach ($d as $id => $g) if (json_encode(gene_index_find($ix, $id)) !== json_encode($g)) $bad++;
        foreach (['', 'zm00001eb067740', 'Zm00001eb06774', 'Zm00001eb0677400', 'Zm00001eb067740 ', 'AAAA', 'zzzz'] as $m)
            if (gene_index_find($ix, $m) !== null) $bad++;
        echo json_encode(['open' => true, 'n' => count($d), 'bad' => $bad]);""", ROOT], capture_output=True, text=True)
    gi = json.loads(probe.stdout or '{}')
    check(gi.get('open') and gi.get('n') == 39756 and gi.get('bad') == 0,
          f"gff/genes_index.txt is current and answers every gene of the store like the store does ({gi})")

    gdir = os.path.join(tmp, 'gff'); os.makedirs(gdir)
    for f in ('lookupGeneModel.php', 'gene_index_lib.php'):
        shutil.copy(os.path.join(ROOT, f), tmp)
    src, idx = os.path.join(gdir, 'genes_data.serialized'), os.path.join(gdir, 'genes_index.txt')
    shutil.copy(os.path.join(ROOT, 'gff', 'genes_data.serialized'), src)
    shutil.copy(os.path.join(ROOT, 'gff', 'genes_index.txt'), idx)
    look = lambda g: php(tmp, 'lookupGeneModel.php', get={'geneModelId': g})
    want = {'chromosome': 'chr2', 'start': '4493424', 'end': '4497434', 'ID': 'Zm00001eb067740'}
    m0 = os.stat(idx).st_mtime_ns
    check(look('Zm00001eb067740') == want and os.stat(idx).st_mtime_ns == m0
          and look('Zm00001eb06774') == {'chromosome': 'chr1', 'start': '0', 'end': '0', 'id': 'empty'},
          "lookupGeneModel.php answers from the current index and leaves it alone")

    def add_gene(gid, chrom):      # change the store (its byte size changes, so the index is stale)
        subprocess.run([PHP, '-r', '$d = unserialize(file_get_contents($argv[1])); '
                        '$d[$argv[2]] = ["chromosome" => $argv[3], "start" => "100", "end" => "200", "ID" => $argv[2]]; '
                        'file_put_contents($argv[1], serialize($d));', src, gid, chrom], check=True)
    add_gene('Zm00001eb999990', 'chr10')
    r1 = look('Zm00001eb999990')
    head = open(idx).readline()
    check(r1.get('chromosome') == 'chr10' and f"source_bytes={os.path.getsize(src)} records=39757" in head,
          f"store changed: the index is rebuilt ({head.strip()[24:]})")
    add_gene('Zm00001eb999991', 'chr3')
    m1 = os.stat(idx).st_mtime_ns
    os.chmod(gdir, 0o555)
    try:
        r2 = look('Zm00001eb999991'); r2b = look('Zm00001eb067740')
    finally:
        os.chmod(gdir, 0o755)
    check(r2.get('chromosome') == 'chr3' and r2b == want and os.stat(idx).st_mtime_ns == m1,
          "store changed, gff/ read-only: the store itself answers, the stale index is left alone")
    os.remove(idx)
    check(look('Zm00001eb999991').get('chromosome') == 'chr3' and os.path.exists(idx), "missing index: rebuilt")
    with open(idx, 'r+b') as fh:
        fh.truncate(os.path.getsize(idx) - 7)
    check(look('Zm00001eb067740') == want and open(idx).readline() == head.replace('39757', '39758').replace(
          f"source_bytes={head.split('source_bytes=')[1].split()[0]}", f"source_bytes={os.path.getsize(src)}"),
          "truncated index: rebuilt")
    check(not [f for f in os.listdir(gdir) if f.endswith('.tmp')], "no temporary index files left")

    # ---- replies name no server path (patch 0040) ----
    def pathless(j):
        txt = json.dumps(j)
        return not any(k in j for k in ('command', 'output', 'store', 'vcfDir')) and 'hdf5/' not in txt \
            and tmp not in txt and 'h5_to_vcf' not in txt and 'distance/' not in txt
    p1 = php(tmp, 'processForm.php', 'POST', post=dict(win, chr='chr11', genotypes='[]', outName=''))
    p2 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73']), outName=''),
             env={'PYTHON_PATH': '/nonexistent/python'})
    p3 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73']), outName=''))
    check(p1.get('status') == 'error' and 'No variant store for zmgrin2026_imp chr11' in p1.get('message', '') and pathless(p1)
          and p2.get('status') == 'error' and 'extraction failed' in p2.get('message', '') and pathless(p2)
          and p3.get('status') == 'success' and pathless(p3),
          f"processForm.php replies name no path: {p1.get('message')} | {p2.get('message')} | {sorted(p3)}")
    p4 = php(tmp, 'processForm.php', 'POST', post=dict(win, chr='chr11', genotypes='[]', outName=''), env={'SNPTOOLS_DEBUG': '1'})
    p5 = php(tmp, 'processForm.php', 'POST', post=dict(win, genotypes=json.dumps(['ZmG_B73']), outName=''),
             env={'PYTHON_PATH': '/nonexistent/python', 'SNPTOOLS_DEBUG': '1'})
    check('hdf5/' in p4.get('store', '') and 'h5_to_vcf' in p5.get('command', ''), "SNPTOOLS_DEBUG=1 returns the details")
    denv = {'SNPTOOLS_DISTANCE_DIR': os.path.join(tmp, 'distance')}
    i1 = php(tmp, 'ibsCompare.php', get={'focal': 'a', 'dataset': 'nosuchset'}, env=denv)
    i2 = php(tmp, 'ibsCompare.php', get={'focal': 'zz9', 'dataset': 'zz'}, env=denv)
    check(i1.get('error') == 'Genome-wide matrices for nosuchset are not installed on this server.' and pathless(i1)
          and i2.get('error') == 'Focal id not found in the zz matrices.' and pathless(i2),
          f"ibsCompare.php errors name no path: {i1.get('error')} | {i2.get('error')}")
    i3 = php(tmp, 'ibsCompare.php', get={'probe': '1'}, env=denv)
    check(i3.get('dataset') == 'zmgrin2026', f"ibsCompare.php without a dataset serves zmgrin2026: {i3}")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"review fixes: {len(bad)} failed; mismatches={len(bad)}")
sys.exit(1 if bad else 0)
