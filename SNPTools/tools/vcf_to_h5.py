#!/usr/bin/env python3
"""vcf_to_h5.py -- build a SNPTools HDF5 genotype store from a (b)gzipped VCF.

    python3 tools/vcf_to_h5.py <in.vcf[.gz]> <out.h5> [--samples ids.txt] [--chunk 200000] [--h5-chunk 65536]

Writes the layout h5_to_vcf.py / processForm.php read (one file per chromosome;
name it hdf5/version3/<family>_<chr>_<quality>.h5, e.g. zmgrin2026_chr10_impute.h5):

  POS                   int64 [n]            sorted positions
  CHROM, REF, ALT, QUAL vlen bytes [n]
  INFO                  vlen bytes [n]       verbatim VCF INFO (GENEMODEL, TYPE, EFFECT, SUB,
                                             MAXR2, MAF, plantcad1/2_score, ESM*_score, ...)
  <sample_id>           int8 [n]             0 = 0/0, 1 = het, 2 = 1/1, 3 = missing

Multi-allelic records are not expected (release VCFs are biallelic, norm -d all);
any allele index > 0 is treated as ALT. Only one CHROM per input is allowed.
--samples restricts/reorders the sample columns (one id per line; ids absent
from the VCF are skipped with a warning). Rows are buffered --chunk at a time in an
int8 matrix (memory ~ chunk x samples bytes, 186 MB at 200,000 x 933), and datasets are
stored in HDF5 chunks of --h5-chunk positions (default 65,536) with gzip level 4, so a
gene-sized query decompresses one small chunk per sample column.
Requires h5py + numpy.
"""
import argparse, gzip, sys
import numpy as np
import h5py

CODE = {b'0/0': 0, b'0|0': 0, b'0/1': 1, b'1/0': 1, b'0|1': 1, b'1|0': 1,
        b'1/1': 2, b'1|1': 2, b'./.': 3, b'.|.': 3, b'.': 3}


def gt_code(cell):
    c = CODE.get(cell)
    if c is not None:
        return c
    g = cell.split(b':', 1)[0]
    c = CODE.get(g)
    if c is not None:
        return c
    a = g.replace(b'|', b'/').split(b'/')
    if len(a) < 2 or b'.' in a or b'' in a:
        return 3
    return (a[0] != b'0') + (a[1] != b'0')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('vcf'); ap.add_argument('h5')
    ap.add_argument('--samples', help='file with one sample id per line (subset / order)')
    ap.add_argument('--chunk', type=int, default=200_000, help='rows buffered per write')
    ap.add_argument('--h5-chunk', type=int, default=65_536, help='HDF5 chunk length (positions)')
    a = ap.parse_args()
    op = gzip.open if a.vcf.endswith('.gz') else open
    vlen = h5py.special_dtype(vlen=bytes)
    with op(a.vcf, 'rb') as fh, h5py.File(a.h5, 'w') as h5:
        for line in fh:
            if line.startswith(b'#CHROM'):
                cols = line.rstrip(b'\n').split(b'\t'); break
        names = [c.decode() for c in cols[9:]]
        idx = {n: i + 9 for i, n in enumerate(names)}
        want = names
        if a.samples:
            want = [l.strip() for l in open(a.samples) if l.strip()]
            miss = [w for w in want if w not in idx]
            if miss:
                print(f'warning: {len(miss)} requested sample(s) not in VCF, skipped: {miss[:8]}', file=sys.stderr)
            want = [w for w in want if w in idx]
        ci = [idx[w] for w in want]
        cz = dict(chunks=(a.h5_chunk,), compression='gzip', compression_opts=4)
        ds = {k: h5.create_dataset(k, (0,), maxshape=(None,), dtype=vlen, **cz) for k in ('CHROM', 'REF', 'ALT', 'QUAL', 'INFO')}
        ds['POS'] = h5.create_dataset('POS', (0,), maxshape=(None,), dtype='i8', **cz)
        gds = [h5.create_dataset(w, (0,), maxshape=(None,), dtype='i1', **cz) for w in want]
        n, chrom0, lastpos = 0, None, -1
        buf = {k: [] for k in ('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO')}
        G = np.empty((a.chunk, len(gds)), dtype=np.int8)

        def flush():
            nonlocal n
            m = len(buf['POS'])
            if not m:
                return
            for k in ('CHROM', 'REF', 'ALT', 'QUAL', 'INFO'):
                ds[k].resize((n + m,)); ds[k][n:] = np.array(buf[k], dtype=object)
            ds['POS'].resize((n + m,)); ds['POS'][n:] = np.array(buf['POS'], dtype=np.int64)
            for j, d in enumerate(gds):
                d.resize((n + m,)); d[n:] = G[:m, j]
            n += m
            for v in buf.values(): v.clear()

        for line in fh:
            t = line.rstrip(b'\n').split(b'\t')
            if len(t) < 10:
                continue
            if chrom0 is None:
                chrom0 = t[0]
            elif t[0] != chrom0:
                sys.exit(f'more than one CHROM in input ({chrom0!r}, {t[0]!r}); split per chromosome first')
            pos = int(t[1])
            if pos < lastpos:
                sys.exit(f'input not position-sorted at {t[0].decode()}:{pos}')
            lastpos = pos
            buf['CHROM'].append(t[0]); buf['POS'].append(pos); buf['REF'].append(t[3]); buf['ALT'].append(t[4])
            buf['QUAL'].append(t[5]); buf['INFO'].append(t[7])
            G[len(buf['POS']) - 1, :] = [gt_code(t[k]) for k in ci]
            if len(buf['POS']) >= a.chunk:
                flush()
        flush()
        h5.attrs['source_vcf'] = a.vcf
        h5.attrs['n_samples'] = len(gds)
    print(f'{a.h5}: {n} sites x {len(want)} samples')


if __name__ == '__main__':
    main()
