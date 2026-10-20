"""store_windows.py -- the fixture windows, read from an installed full-chromosome store.

The test-region fixtures (fixtures/zmgrin2026_v1.4_<chr>_testregions.vcf.gz) were cut from the
release before the score fields existed: their INFO has TYPE/EFFECT/GENEMODEL/SUB, MAF and
ESM1-3, but no PlantCAD1/2, Evo2, ESM-C or MAXR2. With a full store installed the app reads
those too, so a check that recomputes what the app shows needs them.

window_rows(root, fixture, store, spans) exports the fixture's windows (and any extra [lo, hi]
spans, e.g. a tested gene the windows only partly cover) from the store with h5_to_vcf.py (the
fixture's own samples, in its order) and returns the store's rows
[(chrom, pos, ref, alt, info dict, genotype cells)] sorted by position, and how many lie outside
the fixture -- after checking every fixture row against the store:
the same sites with the same REF/ALT, the same genotypes, and the same value for every INFO key
the fixture carries (numbers compared as numbers). Any difference raises StoreMismatch, so the
score fields taken from the store sit on records the independent release cut confirms.
"""
import gzip, json, os, subprocess, sys, tempfile


class StoreMismatch(Exception):
    pass


def read_vcf(path):
    samples, rows = [], []
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            t = line.rstrip('\n').split('\t')
            if line.startswith('#'):
                samples = t[9:]
                continue
            info = dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in t[7].split(';') if kv)
            rows.append((t[0], int(t[1]), t[3], t[4], info, t[9:]))
    return samples, rows


def dosage(g):
    a = g.split(':')[0].replace('|', '/').split('/')
    if len(a) < 2 or '.' in a or '' in a:
        return None
    return (a[0] != '0') + (a[1] != '0')


def same_value(a, b):
    try:
        return abs(float(a.split(',')[0]) - float(b.split(',')[0])) < 1e-9 and a.count(',') == b.count(',')
    except ValueError:
        return a == b


def windows(positions, gap=50_000):
    """[lo, hi] spans covering the positions, split where they are more than gap apart."""
    out = []
    for p in sorted(positions):
        if out and p - out[-1][1] <= gap:
            out[-1][1] = p
        else:
            out.append([p, p])
    return out


def window_rows(root, fixture, store, spans=()):
    samples, fx = read_vcf(fixture)
    python = os.environ.get('PYTHON_PATH') or sys.executable
    got = {}
    with tempfile.TemporaryDirectory(prefix='snpt_win_') as tmp:
        ids = os.path.join(tmp, 'ids.json')
        json.dump(samples, open(ids, 'w'))
        for i, (lo, hi) in enumerate(windows(r[1] for r in fx) + [list(x) for x in spans]):
            out = os.path.join(tmp, f'w{i}.vcf')
            subprocess.run([python, os.path.join(root, 'h5_to_vcf.py'), store, out, str(lo), str(hi), ids],
                           check=True, capture_output=True, env=dict(os.environ, SNPTOOLS_MAX_CELLS='1e12'))
            if not os.path.exists(out):
                continue
            s2, rows = read_vcf(out)
            if s2 != samples:
                raise StoreMismatch(f'{store}: sample columns differ from {fixture}')
            for r in rows:
                got[(r[1], r[2], r[3])] = r
    for chrom, pos, ref, alt, info, gts in fx:
        r = got.get((pos, ref, alt))
        if r is None:
            raise StoreMismatch(f'{store}: no {chrom}:{pos} {ref}>{alt} (in {os.path.basename(fixture)})')
        if [dosage(g) for g in gts] != [dosage(g) for g in r[5]]:
            raise StoreMismatch(f'{store}: genotypes differ at {chrom}:{pos}')
        for k, v in info.items():
            if not same_value(v, r[4].get(k, '')):
                raise StoreMismatch(f'{store}: INFO {k} differs at {chrom}:{pos}: fixture {v!r}, store {r[4].get(k)!r}')
    extra = len(got) - len(fx)
    return sorted(got.values(), key=lambda r: (r[1], r[2], r[3])), extra
