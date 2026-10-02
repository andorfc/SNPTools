#!/usr/bin/env python3
"""check_annotation_columns.py <results.json> <site_root> -- SNPVersity annotation columns.

(a) The 13 annotation columns (Gene model ... ESM3) are rendered, in the same order, for
    zmgrin2026_imp, mgdb2026_hq and mgdb2026_hc; none disappears when a set lacks a field.
(b) Per set and column, the number of filled cells equals what the INFO predicts
    (zmgrin2026_imp: fixtures/zmgrin2026_v1.4_chr2_testregions.vcf.gz; MaizeGDB 2026:
    fixtures/mgdb2026_{hq,hc}_chr2_4491424_4499434.sites.vcf.gz), and columns a set lacks
    are empty with status na/pending.
(c) zmgrin2026_imp INFO: MAF recomputed from the fixture genotypes; ESM1/2/3 equal the
    rounded llr_esm1b / llr_esm2 / llr_esm3 of fixtures/annotation/grz2023_missense_esm_testregions.tsv.gz
    (llr_esm2 = the store ESM-2 layer, i.e. esm2_store_score); TYPE/EFFECT/GENEMODEL/SUB equal
    the Schnable scored VCF rows. Standard library only."""
import csv, gzip, io, json, sys

res, root = sys.argv[1:3]
A = json.load(open(res))['annot']
FX = f'{root}/localdev/fixtures'
LABELS = ['Gene model', 'Effect', 'SNPEff Impact', 'Domain', 'MQ', 'COMP', 'maxR²', 'MAF',
          'PlantCAD1', 'PlantCAD2', 'ESM1', 'ESM2', 'ESM3']
KEYS = ['gene', 'effect', 'impact', 'domain', 'mq', 'comp', 'r2', 'maf', 'pc1', 'pc2', 'esm1', 'esm2', 'esm3']
bad = 0


def info(s):
    return dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in s.split(';')) if s not in ('.', '') else {}


def present(v):
    return v is not None and v.split(',')[0] not in ('', '.')


def records(path, lo=4491424, hi=4499434):
    for line in gzip.open(path, 'rt'):
        if line[0] == '#':
            continue
        t = line.rstrip('\n').split('\t')
        if lo <= int(t[1]) <= hi:
            yield t


def expected(path, status):
    n = 0; f = dict.fromkeys(KEYS, 0)
    for t in records(path):
        I = info(t[7]); n += 1
        f['gene'] += present(I.get('GENEMODEL')); f['effect'] += present(I.get('TYPE'))
        f['impact'] += 1                                  # pill always shown (MODIFIER default)
        f['mq'] += present(I.get('MQ')); f['comp'] += present(I.get('CVP')); f['r2'] += present(I.get('MAXR2'))
        f['maf'] += present(I.get('MAF'))
        f['pc1'] += present(I.get('plantcad1_score')) or present(I.get('DNA_SCORE'))
        f['pc2'] += present(I.get('plantcad2_score'))
        f['esm1'] += present(I.get('ESM1_score')) or present(I.get('AA_SCORE'))
        f['esm2'] += present(I.get('ESM2_score')); f['esm3'] += present(I.get('ESM3_score'))
    for k in KEYS:
        if status.get(k, 'ok') != 'ok':
            f[k] = 0                                      # unavailable / pending: empty cells
    return n, f


SETS = {'zmgrin2026_imp': f'{FX}/zmgrin2026_v1.4_chr2_testregions.vcf.gz',
        'mgdb2026_hq': f'{FX}/mgdb2026_hq_chr2_4491424_4499434.sites.vcf.gz',
        'mgdb2026_hc': f'{FX}/mgdb2026_hc_chr2_4491424_4499434.sites.vcf.gz'}
summary = {}
for ds, path in SETS.items():
    a = A.get(ds)
    if not a:
        bad += 1; print('no render for', ds); continue
    hdr = a['headers']
    if hdr[:4] != ['CHR', 'POS', 'REF', 'ALT'] or hdr[4:17] != LABELS:
        bad += 1; print('columns differ', ds, hdr)
    status = {c['key']: c['status'] for c in a['cols']}
    n, f = expected(path, status)
    if a['rows'] != n:
        bad += 1; print('rows', ds, a['rows'], n)
    for k in KEYS:
        if a['filled'][k] != f[k]:
            bad += 1; print('filled', ds, k, a['filled'][k], 'expected', f[k])
    if any(v != 'ok' for v in status.values()) != bool(a['note']):
        bad += 1; print('availability note', ds, repr(a['note'][:80]))
    summary[ds] = {k: (status[k] if status[k] != 'ok' else f"{a['filled'][k]}/{n}") for k in KEYS}

# (c) the GRIN-linked INFO against its sources, all 3,495 fixture sites
snp = {}
for t in records(f'{FX}/annotation/schnable_scored_testregions.sites.vcf.gz', 0, 10**10):
    snp[(t[0], t[1], t[3], t[4])] = info(t[7])
esm = {}
for r in csv.DictReader(io.TextIOWrapper(gzip.open(f'{FX}/annotation/grz2023_missense_esm_testregions.tsv.gz')), delimiter='\t'):
    esm.setdefault((r['chrom'], r['pos'], r['ref'], r['alt']), []).append(r)
nsite = nesm = 0
for c in ('chr1', 'chr2', 'chr5', 'chr7', 'chr8', 'chr9', 'chr10'):
    for t in records(f'{FX}/zmgrin2026_v1.4_{c}_testregions.vcf.gz', 0, 10**10):
        nsite += 1; I = info(t[7]); k = (t[0], t[1], t[3], t[4]); S = snp[k]
        nz = lambda v: None if v in (None, '', '.') else v      # a bare '.' is not copied
        if any(I.get(x) != nz(S.get(x)) for x in ('TYPE', 'EFFECT', 'GENEMODEL', 'SUB')):
            bad += 1; print('snpeff', k)
        if any(x in I for x in ('MQ', 'CVP', 'MAXR2', 'DNA_SCORE', 'AA_SCORE', 'plantcad1_score', 'plantcad2_score')):
            bad += 1; print('unexpected key', k)
        al = [a for g in t[9:] for a in g.split(':')[0].replace('|', '/').split('/') if a != '.']
        p = sum(a != '0' for a in al) / len(al)
        if f'{min(p, 1 - p):.4f}' != I.get('MAF'):
            bad += 1; print('maf', k, I.get('MAF'), min(p, 1 - p))
        if k in esm:
            nesm += 1
            types, genes, subs = S['TYPE'].split(','), S['GENEMODEL'].split(','), S['SUB'].split(',')
            j = next(i for i, x in enumerate(types) if 'missense' in x)
            r = next((r for r in esm[k] if r['vcf_protein'] == genes[j] and r['variant'] == subs[j]), esm[k][0])
            for key, col in (('ESM1_score', 'llr_esm1b'), ('ESM2_score', 'llr_esm2'), ('ESM3_score', 'llr_esm3')):
                if I.get(key) != f'{round(float(r[col]), 1):.1f}':
                    bad += 1; print('esm', k, key, I.get(key), r[col])
        elif any(x in I for x in ('ESM1_score', 'ESM2_score', 'ESM3_score')):
            bad += 1; print('esm on non-missense', k)
print('annotation columns:', json.dumps(summary, ensure_ascii=False))
print(f'GRIN-linked INFO: {nsite} sites (SnpEff fields, MAF), {nesm} missense sites with ESM1/2/3; mismatches={bad}')
sys.exit(1 if bad else 0)
