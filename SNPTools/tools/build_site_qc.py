#!/usr/bin/env python3
"""build_site_qc.py -- per-site genotype counts beside each variant store (site QC), and the
release's per-line heterozygosity table (line QC).

    python3 tools/build_site_qc.py [--jobs N] [store.h5 ...]             build the sidecars
    python3 tools/build_site_qc.py --check [store.h5 ...]                verify them (exit 1 on a mismatch)
    python3 tools/build_site_qc.py --summary [--jobs N] [store.h5 ...]   summary + line QC files

With no store named, the ten stores hdf5/grin2026/zmgrin2026_chr{1..10}_impute.h5 are used.

Sidecar. For <name>.h5 the sidecar is <name>.siteqc.h5 in the same folder. The store is opened
read-only and never changed. The sidecar holds, per site of the store and in the store's order,
chunked (65,536,) with gzip level 4 like the store:

  POS     int64   a copy of the store's POS, for alignment checks
  NHET    int16   heterozygous carriers among the store's samples (genotype code 1)
  NHOM    int16   homozygous-alternate carriers (code 2)
  NMISS   int16   missing calls (code 3)

and the attributes format ("snptools-siteqc-1"), store (base name), store_bytes, n_sites,
n_samples and built (ISO date). A sidecar is valid while its n_sites, store_bytes and POS match
its store (--check). The file is written under a temporary name and renamed.

Only counts are stored. The class of a site is derived from NHET and NHOM by one rule, first
match wins (site_qc_codes below):

  NO_CARRIER    NHET + NHOM == 0
  HET_ONLY      NHOM == 0
  HET_EXCESS    NHET > NHOM
  HET_ELEVATED  4 * NHET >= NHET + NHOM
  PASS          everything else

Usable = PASS or HET_ELEVATED. Flagged = HET_ONLY or HET_EXCESS.

Summary (--summary) reads the stores and their sidecars again and writes, under --root:

  data/qc/<family>.siteqc.summary.json   class counts per chromosome and in total, clean sites
  data/qc/<family>.lineqc.tsv            per sample: het_calls, hom_calls, missing_calls,
                                         het_share_clean, class
  js/<prefix>.lineqc.js                  window.SNP_LINE_QC[<family>] for the browser

A clean site is a PASS site with NHOM >= 3. A line's share is its number of heterozygous calls at
clean sites divided by the number of clean sites. Class H (heterozygous sample) when the share is
above 0.05, E (elevated) above 0.02, I (inbred) otherwise, on unrounded values.

Memory: one sample column at a time plus a few per-site vectors (about 150 MB for chr1).
Requires h5py + numpy.
"""
import argparse, datetime, json, os, re, sys, time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import h5py

FIXED = ('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO')
FORMAT = 'snptools-siteqc-1'
CHUNK = 65536
SITEQC_NAMES = ['NO_CARRIER', 'HET_ONLY', 'HET_EXCESS', 'HET_ELEVATED', 'PASS']
RULE = ('First match wins. NO_CARRIER: NHET+NHOM = 0; HET_ONLY: NHOM = 0; HET_EXCESS: NHET > NHOM; '
        'HET_ELEVATED: 4*NHET >= NHET+NHOM; PASS: everything else. Usable = PASS or HET_ELEVATED; '
        'flagged = HET_ONLY or HET_EXCESS.')
CLEAN_RULE = 'PASS and NHOM >= 3'
LINE_THRESHOLDS = {'elevated': 0.02, 'heterozygous': 0.05}
LINE_LABELS = {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'}
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def site_qc_codes(het, hom):            # numpy integer arrays -> int8 codes 0..4
    het = np.asarray(het, dtype=np.int32); hom = np.asarray(hom, dtype=np.int32)
    c = het + hom
    f = np.full(len(het), 4, dtype=np.int8)
    f[het * 4 >= c] = 3
    f[(het > hom) & (hom >= 1)] = 2
    f[(hom == 0) & (het >= 1)] = 1
    f[c == 0] = 0
    return f


def line_class(share):
    if share > LINE_THRESHOLDS['heterozygous']:
        return 'H'
    if share > LINE_THRESHOLDS['elevated']:
        return 'E'
    return 'I'


def sidecar_path(store):
    return (store[:-3] if store.endswith('.h5') else store) + '.siteqc.h5'


def default_stores(root):
    return [os.path.join(root, 'hdf5', 'grin2026', f'zmgrin2026_chr{c}_impute.h5') for c in range(1, 11)]


def samples_of(f, store):
    """The store's sample columns: every top-level dataset but the six fixed columns, each int8
    with the length of POS, as many as the store's n_samples attribute says."""
    n = f['POS'].shape[0]
    ids = [k for k in f.keys() if k not in FIXED]
    for k in ids:
        d = f[k]
        if not isinstance(d, h5py.Dataset) or d.dtype != np.int8 or d.shape != (n,):
            raise SystemExit(f'{store}: column {k} is not an int8 dataset of length {n}')
    if 'n_samples' in f.attrs and int(f.attrs['n_samples']) != len(ids):
        raise SystemExit(f'{store}: {len(ids)} sample columns, but the store says n_samples = {int(f.attrs["n_samples"])}')
    return ids


def chrom_of(f, store):
    if f['POS'].shape[0]:
        c = f['CHROM'][0]
        return c.decode() if isinstance(c, bytes) else str(c)
    m = re.search(r'_(chr[0-9A-Za-z]+)_', os.path.basename(store))
    return m.group(1) if m else os.path.basename(store)


def check_sidecar(store, read_pos=True):
    """(ok, reason). A sidecar is valid for its store while n_sites, store_bytes and POS agree."""
    sc = sidecar_path(store)
    if not os.path.isfile(sc):
        return False, f'no sidecar {os.path.basename(sc)}'
    try:
        with h5py.File(sc, 'r') as s, h5py.File(store, 'r') as f:
            a = s.attrs
            if a.get('format') != FORMAT:
                return False, f'format {a.get("format")!r}, expected {FORMAT!r}'
            n = f['POS'].shape[0]
            if int(a['n_sites']) != n:
                return False, f'n_sites {int(a["n_sites"])} != {n} sites in the store'
            if int(a['store_bytes']) != os.path.getsize(store):
                return False, f'store_bytes {int(a["store_bytes"])} != store size {os.path.getsize(store)}'
            for k in ('NHET', 'NHOM', 'NMISS', 'POS'):
                if s[k].shape != (n,):
                    return False, f'{k} has shape {s[k].shape}, expected ({n},)'
            if read_pos and not np.array_equal(s['POS'][:], f['POS'][:].astype(np.int64)):
                return False, 'POS differs from the store'
    except (OSError, KeyError) as e:
        return False, f'unreadable: {e}'
    return True, 'ok'


def build_one(store):
    """Count het / hom-alt / missing per site, one sample column at a time, and write the sidecar."""
    t0 = time.time()
    store_bytes = os.path.getsize(store)
    with h5py.File(store, 'r') as f:
        pos = f['POS'][:].astype(np.int64)
        n = len(pos)
        ids = samples_of(f, store)
        het = np.zeros(n, dtype=np.int32); hom = np.zeros(n, dtype=np.int32); miss = np.zeros(n, dtype=np.int32)
        for k in ids:
            col = f[k][:]
            het += col == 1; hom += col == 2; miss += col == 3
        chrom = chrom_of(f, store)
    if len(ids) > np.iinfo(np.int16).max:
        raise SystemExit(f'{store}: {len(ids)} samples do not fit the int16 counts')
    sc = sidecar_path(store)
    tmp = f'{sc}.tmp-{os.getpid()}'
    cz = dict(chunks=(min(CHUNK, n),), compression='gzip', compression_opts=4) if n else {}
    try:
        with h5py.File(tmp, 'w') as s:
            s.create_dataset('POS', data=pos, dtype='i8', **cz)
            s.create_dataset('NHET', data=het.astype(np.int16), **cz)
            s.create_dataset('NHOM', data=hom.astype(np.int16), **cz)
            s.create_dataset('NMISS', data=miss.astype(np.int16), **cz)
            s.attrs['format'] = FORMAT
            s.attrs['store'] = os.path.basename(store)
            s.attrs['store_bytes'] = store_bytes
            s.attrs['n_sites'] = n
            s.attrs['n_samples'] = len(ids)
            s.attrs['built'] = datetime.date.today().isoformat()
        os.replace(tmp, sc)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    cnt = np.bincount(site_qc_codes(het, hom), minlength=5)
    return dict(store=os.path.basename(store), chrom=chrom, sites=n, samples=len(ids),
                classes={SITEQC_NAMES[i]: int(cnt[i]) for i in range(5)},
                sidecar_bytes=os.path.getsize(sc), seconds=round(time.time() - t0, 1))


def summarize_one(store):
    """Class counts from the sidecar, and per sample: het / hom / missing calls over all sites and
    het calls at clean sites (PASS with NHOM >= 3)."""
    t0 = time.time()
    ok, why = check_sidecar(store, read_pos=False)
    if not ok:
        raise SystemExit(f'{store}: {why} (run tools/build_site_qc.py first)')
    sc = sidecar_path(store)
    with h5py.File(sc, 'r') as s:
        nhet = s['NHET'][:].astype(np.int32); nhom = s['NHOM'][:].astype(np.int32)
    codes = site_qc_codes(nhet, nhom)
    clean = (codes == 4) & (nhom >= 3)
    cnt = np.bincount(codes, minlength=5)
    per = {}
    with h5py.File(store, 'r') as f:
        ids = samples_of(f, store)
        chrom = chrom_of(f, store)
        for k in ids:
            col = f[k][:]
            h = col == 1
            per[k] = (int(np.count_nonzero(h)), int(np.count_nonzero(col == 2)), int(np.count_nonzero(col == 3)),
                      int(np.count_nonzero(h & clean)))
    return dict(store=os.path.basename(store), store_bytes=os.path.getsize(store), sidecar=os.path.basename(sc),
                sidecar_bytes=os.path.getsize(sc), chrom=chrom, sites=int(len(codes)),
                classes={SITEQC_NAMES[i]: int(cnt[i]) for i in range(5)}, clean=int(clean.sum()),
                per=per, seconds=round(time.time() - t0, 1))


def run_all(fn, stores, jobs):
    if jobs <= 1 or len(stores) == 1:
        for st in stores:
            yield fn(st)
        return
    with ProcessPoolExecutor(max_workers=min(jobs, len(stores))) as ex:
        yield from ex.map(fn, stores)


def chrom_key(c):
    m = re.match(r'chr(\d+)$', c)
    return (0, int(m.group(1)), c) if m else (1, 0, c)


def write_summary(results, root, family, prefix, release, stores):
    results = sorted(results, key=lambda r: chrom_key(r['chrom']))
    ids = list(results[0]['per'])
    for r in results[1:]:
        if set(r['per']) != set(ids):
            raise SystemExit(f'{r["store"]}: its samples differ from {results[0]["store"]}')
    n_clean = sum(r['clean'] for r in results)
    if not n_clean:
        raise SystemExit('no clean site (PASS with NHOM >= 3): cannot compute line shares')
    total = {k: sum(r['classes'][k] for r in results) for k in SITEQC_NAMES}
    lines = []
    for k in ids:
        het = sum(r['per'][k][0] for r in results); hom = sum(r['per'][k][1] for r in results)
        miss = sum(r['per'][k][2] for r in results); hc = sum(r['per'][k][3] for r in results)
        share = hc / n_clean
        lines.append((k, het, hom, miss, share, line_class(share)))
    today = datetime.date.today().isoformat()
    os.makedirs(os.path.join(root, 'data', 'qc'), exist_ok=True)
    counts = {c: sum(1 for x in lines if x[5] == c) for c in ('H', 'E', 'I')}
    summary = {
        'format': 'snptools-siteqc-summary-1', 'family': family, 'release': release, 'built': today,
        'rule': RULE, 'codes': SITEQC_NAMES, 'usable': ['PASS', 'HET_ELEVATED'], 'flagged': ['HET_ONLY', 'HET_EXCESS'],
        'counts_over': f'all {len(ids)} samples of the release',
        'by_chromosome': {r['chrom']: dict(sites=r['sites'], **r['classes'], clean=r['clean']) for r in results},
        'total': dict(sites=sum(r['sites'] for r in results), **total),
        'clean_site': CLEAN_RULE, 'clean_sites': n_clean,
        'line_qc': {'thresholds': LINE_THRESHOLDS, 'labels': LINE_LABELS, 'counts': counts,
                    'share': 'heterozygous calls at clean sites / clean sites'},
        'stores': {r['chrom']: {'store': r['store'], 'store_bytes': r['store_bytes'], 'sidecar': r['sidecar'],
                                'sidecar_bytes': r['sidecar_bytes']} for r in results},
    }
    sp = os.path.join(root, 'data', 'qc', f'{family}.siteqc.summary.json')
    with open(sp, 'w') as fh:
        json.dump(summary, fh, indent=1)
        fh.write('\n')
    tp = os.path.join(root, 'data', 'qc', f'{family}.lineqc.tsv')
    with open(tp, 'w') as fh:
        fh.write('sample\thet_calls\thom_calls\tmissing_calls\thet_share_clean\tclass\n')
        for k, het, hom, miss, share, cl in lines:
            fh.write(f'{k}\t{het}\t{hom}\t{miss}\t{share:.6f}\t{cl}\n')
    jp = os.path.join(root, 'js', f'{prefix}.lineqc.js')
    obj = {'release': release, 'cleanSites': n_clean, 'thresholds': LINE_THRESHOLDS,
           'lines': {k: [round(share, 5), cl] for k, _, _, _, share, cl in lines}}
    with open(jp, 'w') as fh:
        fh.write(f'/* Generated by tools/build_site_qc.py --summary on {today} from {len(results)} stores and their\n'
                 f'   .siteqc.h5 sidecars ({release}). Do not edit by hand: run `make site-qc`.\n'
                 f'   Per sample: [share of clean sites ({CLEAN_RULE}) where the sample is heterozygous,\n'
                 f'   class I = inbred, E = share > {LINE_THRESHOLDS["elevated"]}, H = share > {LINE_THRESHOLDS["heterozygous"]}]. */\n')
        fh.write('window.SNP_LINE_QC = window.SNP_LINE_QC || {};\n')
        fh.write(f'window.SNP_LINE_QC[{json.dumps(family)}] = ' + json.dumps(obj, separators=(',', ':')) + ';\n')
    return summary, (sp, tp, jp)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('stores', nargs='*', help='variant stores (default: the ten zmgrin2026 stores)')
    ap.add_argument('--check', action='store_true', help='verify existing sidecars against their stores')
    ap.add_argument('--summary', action='store_true', help='write the summary and the line QC files')
    ap.add_argument('--jobs', type=int, default=1, help='stores processed in parallel (default 1)')
    ap.add_argument('--root', default=ROOT, help='SNPTools web root for --summary outputs (default: this tree)')
    ap.add_argument('--family', default='zmgrin2026')
    ap.add_argument('--prefix', default='zmgrin', help='js/<prefix>.lineqc.js')
    ap.add_argument('--release', default='zmgrin2026 v1.4')
    a = ap.parse_args()
    stores = a.stores or default_stores(a.root)
    missing = [s for s in stores if not os.path.isfile(s)]
    if missing:
        raise SystemExit('no such store: ' + ', '.join(missing))
    t0 = time.time()
    if a.check:
        bad = 0
        for st in stores:
            ok, why = check_sidecar(st)
            bad += not ok
            print(f'{"ok  " if ok else "FAIL"} {os.path.basename(sidecar_path(st))}: {why}', flush=True)
        print(f'{len(stores) - bad} of {len(stores)} sidecars valid ({time.time() - t0:.0f} s)')
        sys.exit(1 if bad else 0)
    if a.summary:
        res = []
        for r in run_all(summarize_one, stores, a.jobs):
            print(f'{r["chrom"]}: {r["sites"]:,} sites, {r["clean"]:,} clean, {r["seconds"]} s', flush=True)
            res.append(r)
        summary, paths = write_summary(res, a.root, a.family, a.prefix, a.release, stores)
        print('total', summary['total'], 'clean', summary['clean_sites'], 'lines', summary['line_qc']['counts'])
        for p in paths:
            print('wrote', os.path.relpath(p, a.root))
        print(f'{time.time() - t0:.0f} s')
        return
    tot = dict.fromkeys(SITEQC_NAMES, 0); size = 0
    for r in run_all(build_one, stores, a.jobs):
        print(f'{r["chrom"]}: {r["sites"]:,} sites x {r["samples"]} samples -> {sidecar_path(r["store"])} '
              f'({r["sidecar_bytes"] / 1e6:.1f} MB, {r["seconds"]} s) {r["classes"]}', flush=True)
        for k in SITEQC_NAMES:
            tot[k] += r['classes'][k]
        size += r['sidecar_bytes']
    print(f'total {tot}; sidecars {size / 1e6:.1f} MB; {time.time() - t0:.0f} s')


if __name__ == '__main__':
    main()
