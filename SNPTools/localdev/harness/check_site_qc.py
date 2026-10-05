#!/usr/bin/env python3
"""check_site_qc.py <results.json> <site_root> -- the site QC counts every VCF carries.

Every row's INFO ends with NHET, NHOM and SITEQC: a store built from a release VCF that has
them (release v1.4.2) carries them itself, and for a store that does not h5_to_vcf.py appends
them from the store's sidecar (<store>.siteqc.h5, tools/build_site_qc.py). Recomputed here,
independently:
1. The page's queries (run_scenarios.js, Data.queryVariants -> processForm.php) built VCFs.
2. In the all-lines VCFs every row's NHET / NHOM equal the heterozygous / homozygous-alternate
   cells of its 933 genotypes, and SITEQC follows the rule (written again below).
3. The five-line VCF carries the same INFO as the all-lines VCF of the same interval: the
   counts are panel-wide, not selection-wide.
4. With SNPTOOLS_SITEQC=0 the VCF equals the normal one minus the three fields and their three
   ##INFO lines, whether the fields came from the store's own INFO or from the sidecar.
5. In a temporary root that symlinks the store: with no sidecar, or one whose store_bytes is
   wrong, a Note: line says so, and a valid copy of the sidecar gives the normal VCF. Without a
   usable sidecar a store that carries the fields (read from the store itself) still writes the
   normal VCF, its rows as stored; one that does not writes the SNPTOOLS_SITEQC=0 one.
   SNPTOOLS_SITEQC=0 beside the stale sidecar gives the VCF without the fields and no Note.
5b. The merge itself, on a store without the three fields whatever stores are installed: one
   built in a temporary root from the chr9 release fixture. No sidecar: a Note, every INFO as in
   the fixture. Sidecar built: NHET / NHOM / SITEQC on every row equal the fixture's 933
   genotypes, the rest unchanged. SNPTOOLS_SITEQC=0, and a sidecar with a wrong store_bytes
   (with its Note): the VCF without them again.
6. The Python implementations of the rule (build_site_qc.py, h5_to_vcf.py, annotate_release_info.py)
   agree with the rule on every (het, hom) pair from 0 to 60.
7. tools/annotate_release_info.py, run on each fixture VCF with the inputs in
   fixtures/annotation/, writes NHET / NHOM / SITEQC equal to the fixture's own genotype counts and
   to the sidecar values at the same sites, and leaves every other INFO key of the fixture as it was
   (ESMC_score, which the fixtures predate, is the one other key it adds). A store built from its
   chr9 output carries the three keys: h5_to_vcf.py writes its rows as stored, with a sidecar and
   without one, and with SNPTOOLS_SITEQC=0 every row without the three, its other keys untouched.
8. SNPTOOLS_SITEQC=0 takes the three fields out wherever they stand in a stored INFO (last, first,
   alone, apart, in another order) and leaves keys with similar names and rows without them as
   stored: a store of written-out layouts, each with the INFO expected.
Standard library only; h5py work (sidecar values, the stale copy, a store's own INFO) runs under
PYTHON_PATH."""
import glob, gzip, json, os, re, shutil, subprocess, sys, tempfile

res, root = sys.argv[1:3]
root = os.path.abspath(root)
R = json.load(open(res))
Q = R['siteqc']
HERE = os.path.dirname(os.path.abspath(__file__))
FX = os.path.join(root, 'localdev', 'fixtures')
PY = os.environ.get('PYTHON_PATH') or 'python3'
PHP = os.environ.get('PHP_BIN') or 'php'
NAMES = ['NO_CARRIER', 'HET_ONLY', 'HET_EXCESS', 'HET_ELEVATED', 'PASS']
QC_IDS = ('NHET', 'NHOM', 'SITEQC')
SUFFIX = re.compile(r'(?:^|;)NHET=(\d+);NHOM=(\d+);SITEQC=([A-Z_]+)$')
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond:
        fails.append(msg)


def rule(het, hom):        # the plan's rule, first match wins
    if het + hom == 0: return 'NO_CARRIER'
    if hom == 0: return 'HET_ONLY'
    if het > hom: return 'HET_EXCESS'
    if 4 * het >= het + hom: return 'HET_ELEVATED'
    return 'PASS'


def dose(cell):            # vcf_to_h5.py gt_code: None = missing
    a = cell.split(':')[0].replace('|', '/').split('/')
    if len(a) < 2 or '.' in a or '' in a: return None
    return (a[0] != '0') + (a[1] != '0')


def read_vcf(path):
    op = gzip.open if path.endswith('.gz') else open
    meta, head, rows = [], None, []
    with op(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('##'):
                meta.append(line)
            elif line.startswith('#'):
                head = line[1:].split('\t')
            else:
                rows.append(line.split('\t'))
    return meta, head, rows


def info_ids(meta):
    return [m[len('##INFO=<ID='):].split(',', 1)[0] for m in meta if m.startswith('##INFO=<ID=')]


def strip_qc(meta, head, rows):
    """The VCF as it would be without the merge: the three fields and their ##INFO lines removed."""
    meta = [m for m in meta if not m.startswith('##fileDate') and m[len('##INFO=<ID='):].split(',', 1)[0] not in QC_IDS]
    out = []
    for r in rows:
        r = list(r)
        r[7] = SUFFIX.sub('', r[7]) or '.'
        out.append(r)
    return meta, head, out


def no_date(meta):
    return [m for m in meta if not m.startswith('##fileDate')]


def vpath(v):
    return os.path.join(root, v[2:] if v.startswith('./') else v)


# ---------------- 1-4: the page's queries ----------------
V = {}
for name, q in sorted(Q['queries'].items()):
    ok = q.get('vcf') and os.path.isfile(vpath(q['vcf'])) and q.get('rows', 0) > 0 and not q.get('error')
    check(bool(ok), f"{name}: {q.get('chr')}:{q.get('start')}-{q.get('end')} -> {q.get('rows')} rows x {q.get('accs')} lines")
    if ok:
        V[name] = read_vcf(vpath(q['vcf']))
check(not Q.get('consoleErrors'), f"no script errors ({(Q.get('consoleErrors') or [])[:2]})")

classes = {}
for name in ('five_all', 'gene_all', 'region_all'):
    if name not in V: continue
    meta, head, rows = V[name]
    ids = info_ids(meta)
    check(all(ids.count(k) == 1 for k in QC_IDS) and len(head) - 9 == 933,
          f"{name}: one ##INFO line each for NHET, NHOM, SITEQC; {len(head) - 9} sample columns")
    bad, cnt = [], dict.fromkeys(NAMES, 0)
    for r in rows:
        m = SUFFIX.search(r[7])
        d = [dose(c) for c in r[9:]]
        het, hom = d.count(1), d.count(2)
        if not m or (int(m.group(1)), int(m.group(2)), m.group(3)) != (het, hom, rule(het, hom)):
            bad.append((r[1], r[7][-60:], het, hom))
        else:
            cnt[m.group(3)] += 1
    classes[name] = cnt
    check(not bad, f"{name}: NHET / NHOM / SITEQC equal the 933 genotypes on all {len(rows)} rows {cnt} {bad[:2]}")

if 'five' in V and 'five_all' in V:
    a, b = V['five'][2], V['five_all'][2]
    same = len(a) == len(b) and all(x[:8] == y[:8] for x, y in zip(a, b))
    check(same and len(V['five'][1]) - 9 == 5, f"five lines: INFO equals the all-lines query of the interval on all {len(a)} rows (panel-wide counts)")

for name in ('five', 'gene_all', 'region_all'):
    if name in V and name + '_off' in V:
        on, off = V[name], V[name + '_off']
        exp = strip_qc(*on)
        got = (no_date(off[0]), off[1], off[2])
        check(exp == got and not any(k in info_ids(off[0]) for k in QC_IDS) and all(SUFFIX.search(r[7]) for r in on[2]),
              f"{name}: SNPTOOLS_SITEQC=0 output = normal output minus the three fields and ##INFO lines ({len(off[2])} rows)")

# ---------------- 5: stale, missing and valid sidecars in a temporary root ----------------
def php(rt, post, env=None):
    req = {'method': 'POST', 'get': {}, 'post': post, 'root': rt, 'script': 'processForm.php'}
    out = subprocess.run([PHP, os.path.join(HERE, 'php_shim.php')], capture_output=True, text=True,
                         env=dict(os.environ, SNPT_REQ=json.dumps(req), SNPTOOLS_DEBUG='1', **(env or {})))
    try:
        return json.loads(out.stdout)
    except ValueError:
        return {'raw': out.stdout[:300], 'stderr': out.stderr[:300]}


def temp_root():
    """A site root of its own: the two scripts, hdf5/version3/ and vcf/. Returns it and its store folder."""
    tmp = tempfile.mkdtemp(prefix='snpt_siteqc_')
    for f in ('processForm.php', 'h5_to_vcf.py'):
        shutil.copy(os.path.join(root, f), tmp)
    hd = os.path.join(tmp, 'hdf5', 'version3'); os.makedirs(hd); os.makedirs(os.path.join(tmp, 'vcf'))
    return tmp, hd


def request(tmp, post, env=None):
    """One processForm.php request in a temporary root: (its VCF or None, the extractor's Note: lines, the reply)."""
    j = php(tmp, post, env)
    p = os.path.join(tmp, j.get('outFile', '')[2:]) if j.get('outFile') else ''
    if j.get('status') != 'success' or not os.path.isfile(p):
        return None, [], j
    return read_vcf(p), [l for l in (j.get('output') or '').splitlines() if l.startswith('Note:')], j


def run(tmp, post, label, want, note, env=None):
    """The request's VCF equals `want` (the file date apart) and a Note: line holds `note`, or there is none."""
    got, noted, j = request(tmp, post, env)
    if got is None:
        check(False, f"{label}: no VCF ({j})"); return
    ok = (no_date(got[0]), got[1], got[2]) == (no_date(want[0]), want[1], want[2])
    check(ok and (any(note in l for l in noted) if note else not noted),
          f"{label}: VCF {'equals' if ok else 'DIFFERS from'} the expected one; {noted[0][:110] if noted else 'no Note: line'}")


def wrong_store_bytes(side):
    bump = subprocess.run([PY, '-c', 'import h5py,sys\nf=h5py.File(sys.argv[1],"r+")\nf.attrs["store_bytes"]=int(f.attrs["store_bytes"])+1\nf.close()', side],
                          capture_output=True, text=True)
    check(bump.returncode == 0, f"sidecar copy given a wrong store_bytes {bump.stderr[:200]}")


# how many stored rows of an interval carry SITEQC in the store's own INFO (read from the store, not a VCF)
stored_snippet = r'''
import h5py, sys
import numpy as np
f = h5py.File(sys.argv[1], 'r'); pos = f['POS'][:].astype(np.int64)
a, b = int(np.searchsorted(pos, int(sys.argv[2]), 'left')), int(np.searchsorted(pos, int(sys.argv[3]), 'right'))
print(sum(b'SITEQC=' in (v if isinstance(v, bytes) else str(v).encode()) for v in f['INFO'][a:b]), b - a)
'''

q5 = Q['queries'].get('five')
if q5 and 'five' in V and 'five_off' in V:
    store = f"zmgrin2026_{q5['chr']}_impute"
    real_store = os.path.join(root, 'hdf5', 'version3', store + '.h5')
    real_side = os.path.join(root, 'hdf5', 'version3', store + '.siteqc.h5')
    sn = subprocess.run([PY, '-c', stored_snippet, real_store, str(q5['start']), str(q5['end'])], capture_output=True, text=True)
    try:
        n_qc, n_rows = [int(x) for x in sn.stdout.split()]
    except ValueError:
        n_qc = n_rows = -1
    native = n_qc == n_rows > 0
    check(n_rows == len(V['five'][2]) and n_qc in (0, n_rows),
          f"{store}: {n_qc} of {n_rows} stored rows of the interval carry SITEQC "
          f"({'the store carries the three fields itself' if native else 'the fields come from the sidecar'}) {sn.stderr[-200:]}")
    # without a usable sidecar: a store that carries the fields writes its rows as stored (the normal
    # VCF), one that does not writes the VCF of SNPTOOLS_SITEQC=0
    bare, kind = (V['five'], ' (fields in the store: rows as stored)') if native else (V['five_off'], '')
    tmp, hd = temp_root()
    try:
        os.symlink(real_store, os.path.join(hd, store + '.h5'))
        side = os.path.join(hd, store + '.siteqc.h5')
        post = {'chr': q5['chr'], 'start': str(q5['start']), 'end': str(q5['end']), 'dataSet': 'zmgrin2026_imp',
                'genotypes': json.dumps(V['five'][1][9:]), 'outName': ''}
        run(tmp, post, 'no sidecar' + kind, bare, 'no site QC sidecar')
        shutil.copy(real_side, side)
        run(tmp, post, 'valid sidecar copy', V['five'], None)
        wrong_store_bytes(side)
        run(tmp, post, 'sidecar with a wrong store_bytes' + kind, bare, 'does not match its store')
        run(tmp, post, 'SNPTOOLS_SITEQC=0 beside the stale sidecar', V['five_off'], None, {'SNPTOOLS_SITEQC': '0'})
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

# ---------------- 5b: the merge, on a store without the three fields ----------------
# The release v1.4.2 stores carry NHET / NHOM / SITEQC in their own INFO, so the sidecar adds nothing
# to their rows. A store built here from the chr9 release fixture, whose INFO has none of the three,
# keeps the merge and the three sidecar states under test whatever stores are installed.
fx9 = os.path.join(FX, 'zmgrin2026_v1.4_chr9_testregions.vcf.gz')
m9, h9, r9 = read_vcf(fx9)
tmp, hd = temp_root()
try:
    st = os.path.join(hd, 'zmgrin2026_chr9_impute.h5')
    b9 = subprocess.run([PY, os.path.join(root, 'tools', 'vcf_to_h5.py'), fx9, st], capture_output=True, text=True)
    check(b9.returncode == 0 and os.path.isfile(st) and not any(k + '=' in r[7] for r in r9 for k in QC_IDS),
          f"plain store: built from the chr9 release fixture, whose INFO has none of the three fields ({len(r9)} sites) {b9.stderr[-200:]}")
    post = {'chr': 'chr9', 'start': r9[0][1], 'end': r9[-1][1], 'dataSet': 'zmgrin2026_imp', 'genotypes': json.dumps(h9[9:14]), 'outName': ''}
    bare, noted, j = request(tmp, post)
    if bare is None:
        check(False, f"plain store, no sidecar: no VCF ({j})")
    else:
        site = lambda rows: [(r[1], r[3], r[4], r[7]) for r in rows]
        check(site(bare[2]) == site(r9) and not any(k in info_ids(bare[0]) for k in QC_IDS) and any('no site QC sidecar' in l for l in noted),
              f"plain store, no sidecar: every INFO as in the fixture, no ##INFO line for the three ({len(bare[2])} rows); "
              f"{noted[0][:110] if noted else 'no Note: line'}")
        s9 = subprocess.run([PY, os.path.join(root, 'tools', 'build_site_qc.py'), st], capture_output=True, text=True)
        merged, noted, j = request(tmp, post)
        if s9.returncode or merged is None:
            check(False, f"plain store, sidecar built: no VCF ({s9.stderr[-200:]} {j})")
        else:
            bad9 = []
            for r, f in zip(merged[2], r9):
                d = [dose(c) for c in f[9:]]
                het, hom = d.count(1), d.count(2)
                mm = SUFFIX.search(r[7])
                if not mm or (int(mm.group(1)), int(mm.group(2)), mm.group(3)) != (het, hom, rule(het, hom)):
                    bad9.append(r[1])
            check(len(merged[2]) == len(r9) and not bad9 and all(info_ids(merged[0]).count(k) == 1 for k in QC_IDS) and not noted
                  and strip_qc(*merged) == (no_date(bare[0]), bare[1], bare[2]),
                  f"plain store, sidecar built: NHET / NHOM / SITEQC merged on all {len(merged[2])} rows equal the fixture's 933 genotypes, "
                  f"the rest as without it; no Note: line {bad9[:3]}")
            run(tmp, post, 'plain store, SNPTOOLS_SITEQC=0', bare, None, {'SNPTOOLS_SITEQC': '0'})
            wrong_store_bytes(st[:-3] + '.siteqc.h5')
            run(tmp, post, 'plain store, sidecar with a wrong store_bytes', bare, 'does not match its store')
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# ---------------- 6: the rule, in every Python copy ----------------
snippet = r'''
import ast, json, os, sys
import numpy as np
root = sys.argv[1]
sys.path.insert(0, os.path.join(root, 'tools'))
import build_site_qc, annotate_release_info
src = open(os.path.join(root, 'h5_to_vcf.py')).read()
fn = [n for n in ast.parse(src).body if isinstance(n, ast.FunctionDef) and n.name == 'site_qc_codes'][0]
ns = {'np': np}; exec(compile(ast.Module(body=[fn], type_ignores=[]), 'h5_to_vcf.py', 'exec'), ns)
het = np.repeat(np.arange(61), 61); hom = np.tile(np.arange(61), 61)
names = build_site_qc.SITEQC_NAMES
out = {'build_site_qc': [names[c] for c in build_site_qc.site_qc_codes(het, hom)],
       'h5_to_vcf': [names[c] for c in ns['site_qc_codes'](het.astype(np.int16), hom.astype(np.int16))],
       'annotate_release_info': [annotate_release_info.site_qc(int(a), int(b)) for a, b in zip(het, hom)],
       'pairs': [[int(a), int(b)] for a, b in zip(het, hom)], 'names': names}
print(json.dumps(out))
'''
p6 = subprocess.run([PY, '-c', snippet, root], capture_output=True, text=True)
try:
    J = json.loads(p6.stdout)
    want = [rule(a, b) for a, b in J['pairs']]
    check(len(want) == 3721 and J['names'] == NAMES, f"rule table: {len(want)} pairs (0-60 x 0-60), codes {J['names']}")
    for impl in ('build_site_qc', 'h5_to_vcf', 'annotate_release_info'):
        diff = [J['pairs'][i] for i, (x, y) in enumerate(zip(J[impl], want)) if x != y]
        check(not diff and len(J[impl]) == 3721, f"{impl}: agrees with the rule on all 3,721 pairs {diff[:3]}")
except (ValueError, KeyError) as e:
    check(False, f"rule implementations could not be loaded: {e} {p6.stderr[-300:]}")

# ---------------- 7: annotate_release_info.py on the fixtures ----------------
AN = os.path.join(FX, 'annotation')
side_snippet = r'''
import h5py, json, sys
import numpy as np
store, want = sys.argv[1], json.load(sys.stdin)
f = h5py.File(store, 'r'); s = h5py.File(store[:-3] + '.siteqc.h5', 'r')
pos = s['POS'][:]; het = s['NHET']; hom = s['NHOM']
out = {}
for key in want:
    p, ref, alt = key.split(':')
    i = int(np.searchsorted(pos, int(p)))
    while i < len(pos) and pos[i] == int(p):
        if f['REF'][i].decode() == ref and f['ALT'][i].decode() == alt:
            out[key] = [int(het[i]), int(hom[i])]; break
        i += 1
print(json.dumps(out))
'''
tmp = tempfile.mkdtemp(prefix='snpt_annot_')
try:
    for fx in sorted(glob.glob(os.path.join(FX, 'zmgrin2026_v1.4_chr*_testregions.vcf.gz'))):
        chrom = re.search(r'_(chr\d+)_', os.path.basename(fx)).group(1)
        out = os.path.join(tmp, chrom + '.vcf.gz')
        p = subprocess.run([PY, os.path.join(root, 'tools', 'annotate_release_info.py'), '--vcf', fx,
                            '--snpeff', os.path.join(AN, 'schnable_scored_testregions.sites.vcf.gz'),
                            '--esm', os.path.join(AN, 'grz2023_missense_esm_testregions.tsv.gz'), '--out', out],
                           capture_output=True, text=True)
        if p.returncode != 0:
            check(False, f"annotate_release_info.py on {chrom}: {p.stderr[-300:]}"); continue
        m0, h0, r0 = read_vcf(fx)
        m1, h1, r1 = read_vcf(out)
        own, kept = [], True
        for a, b in zip(r0, r1):
            d = [dose(c) for c in a[9:]]
            het, hom = d.count(1), d.count(2)
            mm = SUFFIX.search(b[7])
            if a[:7] != b[:7] or not mm or (int(mm.group(1)), int(mm.group(2)), mm.group(3)) != (het, hom, rule(het, hom)):
                own.append(a[1])
            # every key the fixture has keeps its value; the only keys added are the three, and
            # ESMC_score, which the script writes since the fixtures were cut (2026-09-29)
            ia = dict(kv.split('=', 1) for kv in a[7].split(';') if '=' in kv)
            ib = dict(kv.split('=', 1) for kv in b[7].split(';') if '=' in kv)
            if any(ib.get(k) != v for k, v in ia.items()) or set(ib) - set(ia) - set(QC_IDS) - {'ESMC_score'}:
                kept = False
        check(len(r0) == len(r1) and not own and all(info_ids(m1).count(k) == 1 for k in QC_IDS),
              f"annotate_release_info.py {chrom}: NHET / NHOM / SITEQC equal the fixture's genotypes on {len(r1)} rows {own[:3]}")
        check(kept, f"annotate_release_info.py {chrom}: every INFO key of the fixture unchanged")
        store = os.path.join(root, 'hdf5', 'version3', f'zmgrin2026_{chrom}_impute.h5')
        if not os.path.isfile(store[:-3] + '.siteqc.h5'):
            check(False, f"{chrom}: no sidecar beside {os.path.basename(store)} to compare with"); continue
        keys = [f'{r[1]}:{r[3]}:{r[4]}' for r in r1]
        s = subprocess.run([PY, '-c', side_snippet, store], input=json.dumps(keys), capture_output=True, text=True)
        try:
            S = json.loads(s.stdout)
        except ValueError:
            check(False, f"{chrom}: sidecar values unreadable {s.stderr[-300:]}"); continue
        diff = [k for k, r in zip(keys, r1) if S.get(k) != [int(x) for x in SUFFIX.search(r[7]).groups()[:2]]]
        check(len(S) == len(keys) and not diff, f"annotate_release_info.py {chrom}: counts equal the sidecar at all {len(S)} sites {diff[:3]}")
        if chrom != 'chr9':
            continue
        # a store built from the annotated VCF carries the three keys itself: h5_to_vcf.py writes
        # its rows as stored (no second NHET/NHOM/SITEQC), with the sidecar and without one
        st = os.path.join(tmp, 'zmgrin2026_chr9_impute.h5')
        b1 = subprocess.run([PY, os.path.join(root, 'tools', 'vcf_to_h5.py'), out, st], capture_output=True, text=True)
        ids = os.path.join(tmp, 'ids.json'); json.dump(h1[9:], open(ids, 'w'))
        lo, hi = r1[0][1], r1[-1][1]
        for with_side in (False, True):
            if with_side:
                subprocess.run([PY, os.path.join(root, 'tools', 'build_site_qc.py'), st], capture_output=True, text=True)
            ov = os.path.join(tmp, f'native{int(with_side)}.vcf')
            h = subprocess.run([PY, os.path.join(root, 'h5_to_vcf.py'), st, ov, lo, hi, ids], capture_output=True, text=True, cwd=tmp)
            if b1.returncode or h.returncode or not os.path.isfile(ov):
                check(False, f"store built from the annotated chr9 fixture: {b1.stderr[-200:]} {h.stdout[-200:]}"); break
            m2, h2, r2 = read_vcf(ov)
            same = [x[7] for x in r2] == [x[7] for x in r1] and all(x[7].count('SITEQC=') == 1 for x in r2)
            check(same and all(info_ids(m2).count(k) == 1 for k in QC_IDS) and h2[9:] == h1[9:],
                  f"store with the keys in its INFO, {'with' if with_side else 'without'} a sidecar: rows as stored, "
                  f"one NHET/NHOM/SITEQC each, one ##INFO line each ({len(r2)} rows)")
            # the switch on such a store: the three taken out of every row, nothing else changed
            oz = os.path.join(tmp, f'native{int(with_side)}_off.vcf')
            z = subprocess.run([PY, os.path.join(root, 'h5_to_vcf.py'), st, oz, lo, hi, ids], capture_output=True, text=True, cwd=tmp,
                               env=dict(os.environ, SNPTOOLS_SITEQC='0'))
            if z.returncode or not os.path.isfile(oz):
                check(False, f"store with the keys in its INFO, SNPTOOLS_SITEQC=0: {z.stdout[-200:]} {z.stderr[-200:]}"); break
            m3, h3, r3 = read_vcf(oz)
            check(strip_qc(m2, h2, r2) == (no_date(m3), h3, r3) and not any(k in info_ids(m3) for k in QC_IDS)
                  and not [l for l in z.stdout.splitlines() if l.startswith('Note:')],
                  f"store with the keys in its INFO, {'with' if with_side else 'without'} a sidecar, SNPTOOLS_SITEQC=0: the output "
                  f"minus the three fields and ##INFO lines, no Note: line ({len(r3)} rows)")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

# ---------------- 8: SNPTOOLS_SITEQC=0 wherever the three fields stand in a stored INFO ----------------
LAYOUTS = [  # (INFO as stored, INFO written with SNPTOOLS_SITEQC=0)
    ('TYPE=x;MAF=0.0011;NHET=2;NHOM=13;SITEQC=PASS', 'TYPE=x;MAF=0.0011'),        # the builders' layout: the three last
    ('NHET=0;NHOM=0;SITEQC=NO_CARRIER', '.'),                                      # nothing else in the INFO
    ('NHET=2;NHOM=1;SITEQC=HET_EXCESS;TYPE=x', 'TYPE=x'),                          # first
    ('TYPE=x;NHET=1;NHOM=2;SITEQC=HET_ONLY;MAF=0.0011', 'TYPE=x;MAF=0.0011'),      # in the middle
    ('A=1;NHET=2;B=x,y;NHOM=3;FLAG;SITEQC=PASS;D=5', 'A=1;B=x,y;FLAG;D=5'),       # apart
    ('A=1;SITEQC=PASS;NHOM=3;NHET=0', 'A=1'),                                      # another order
    ('SITEQC=HET_ELEVATED;A=1;NHET=1', 'A=1'),
    ('NHOM=5', '.'),                                                               # one of them alone
    ('XNHET=9;NHETX=1;SITEQCS=no', 'XNHET=9;NHETX=1;SITEQCS=no'),                  # other keys with similar names
    ('XNHET=9;NHET=1;NHOM=2;SITEQC=HET_ONLY', 'XNHET=9'),
    ('TYPE=x;NHET;MAF=0.5', 'TYPE=x;NHET;MAF=0.5'),                                # a flag of that name is not the field
    ('TYPE=x;MAF=0.5', 'TYPE=x;MAF=0.5'),                                          # none of the three
    ('.', '.'),
]
tmp = tempfile.mkdtemp(prefix='snpt_layout_')
try:
    src = os.path.join(tmp, 'layouts.vcf')
    with open(src, 'w') as fh:
        fh.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n')
        for i, (stored, _) in enumerate(LAYOUTS):
            fh.write(f'chrT\t{100 + i}\t.\tA\tG\t.\t.\t{stored}\tGT\t0/0\t1/1\n')
    st = os.path.join(tmp, 'zmgrin2026_chrT_impute.h5')
    b8 = subprocess.run([PY, os.path.join(root, 'tools', 'vcf_to_h5.py'), src, st], capture_output=True, text=True)
    ids = os.path.join(tmp, 'ids.json'); json.dump(['S1', 'S2'], open(ids, 'w'))
    outs = {}
    for name, env in (('as stored', {}), ('SNPTOOLS_SITEQC=0', {'SNPTOOLS_SITEQC': '0'})):
        ov = os.path.join(tmp, f'layouts{len(outs)}.vcf')
        p8 = subprocess.run([PY, os.path.join(root, 'h5_to_vcf.py'), os.path.basename(st), ov, '1', '1000', ids], capture_output=True, text=True, cwd=tmp,
                            env=dict(os.environ, **env))
        outs[name] = read_vcf(ov) if p8.returncode == 0 and os.path.isfile(ov) else None
    if b8.returncode or None in outs.values():
        check(False, f"layout store: {b8.stderr[-200:]} {p8.stdout[-200:]} {p8.stderr[-200:]}")
    else:
        a8, z8 = outs['as stored'], outs['SNPTOOLS_SITEQC=0']
        check([r[7] for r in a8[2]] == [s for s, _ in LAYOUTS] and all(info_ids(a8[0]).count(k) == 1 for k in QC_IDS),
              f"layout store, no switch: all {len(a8[2])} INFO values as stored")
        wrong = [(s, r[7], w) for (s, w), r in zip(LAYOUTS, z8[2]) if r[7] != w]
        check(len(z8[2]) == len(LAYOUTS) and not wrong and not any(k in info_ids(z8[0]) for k in QC_IDS)
              and [r[:7] + r[8:] for r in z8[2]] == [r[:7] + r[8:] for r in a8[2]],
              f"layout store, SNPTOOLS_SITEQC=0: the three fields taken out of {len(LAYOUTS)} written-out layouts (last, first, alone, apart, "
              f"another order), similar key names and rows without them as stored {wrong[:2]}")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"site QC: {len(fails)} failed; classes {json.dumps(classes)}; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
