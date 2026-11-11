#!/usr/bin/env python3
"""check_annotation_columns.py <results.json> <site_root> -- SNPVersity annotation columns.

(a) The annotation columns (Gene model ... ESM3) are rendered, in the same order, for
    zmgrin2026_imp (and mgdb2026_hq / mgdb2026_hc, when those sets are offered): a field a set
    lacks stays as an empty column ('na'), except one the call set never records, which is
    left out ('hidden': MQ and COMP for zmgrin2026_imp); every row has as many cells as the
    header. Site QC ('auto') is shown exactly when the result carries NHET/NHOM/SITEQC: on the
    rows whose INFO has them in the release itself (v1.4.2 on: the window cut from the release VCF
    gives the count), and on every row when the queried store has its site-QC sidecar
    (h5_to_vcf.py then adds them where a row has none); the classes themselves are checked by
    check_site_qc_ui.py.
(b) Per set and column, the number of filled cells equals what the INFO predicts
    (zmgrin2026_imp: fixtures/zmgrin2026_v1.4_chr2_testregions.vcf.gz; MaizeGDB 2026:
    fixtures/mgdb2026_{hq,hc}_chr2_4491424_4499434.sites.vcf.gz), and columns a set lacks
    are empty with status na/pending. Domain: a site is filled when data/domains/ (by_chr/<chr>.json,
    else domains.by_chr.json) has a block covering it, the lookup Data.domainAt() does; 0 when
    neither file is installed.
(c) zmgrin2026_imp INFO: MAF recomputed from the fixture genotypes; ESM1/2/3 equal the
    rounded llr_esm1b / llr_esm2 / llr_esm3 of fixtures/annotation/grz2023_missense_esm_testregions.tsv.gz
    (llr_esm2 = the store ESM-2 layer, i.e. esm2_store_score); TYPE/EFFECT/GENEMODEL/SUB equal
    the Schnable scored VCF rows.
(d) The PlantCAD2 column shows the one-pass score (INFO plantcad2_onepass_score, release v1.4.3), not
    plantcad2_score (512 bp): its filled cells are the rows with that key, and in the window cut from
    the release VCF every SNP has it and no indel does. Standard library only."""
import csv, gzip, io, json, os, sys

res, root = sys.argv[1:3]
A = json.load(open(res))['annot']
FX = f'{root}/localdev/fixtures'
LABELS = ['Gene model', 'Effect', 'SNPEff Impact', 'Domain', 'MQ', 'COMP', 'maxR²', 'MAF', 'Site QC',
          'PlantCAD1', 'PlantCAD2', 'Evo2', 'ESM1', 'ESM2', 'ESM3', 'ESM-C']
KEYS = ['gene', 'effect', 'impact', 'domain', 'mq', 'comp', 'r2', 'maf', 'qc', 'pc1', 'pc2op', 'evo2', 'esm1', 'esm2', 'esm3', 'esmc']
# the chr2 store's site-QC sidecar: with it every row of the result carries a class (without it,
# the rows whose INFO has SITEQC in the release itself)
SIDECAR = os.path.exists(f'{root}/hdf5/grin2026/zmgrin2026_chr2_impute.siteqc.h5')
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


_DOM = {}
def domain_blocks(chrom):
    """The Pfam blocks Data.ensureDomains() loads: [g_start, g_end, name, pfam, type] rows."""
    if chrom not in _DOM:
        import os
        one, every = f'{root}/data/domains/by_chr/{chrom}.json', f'{root}/data/domains/domains.by_chr.json'
        _DOM[chrom] = json.load(open(one)) if os.path.exists(one) else \
                      (json.load(open(every)).get(chrom, []) if os.path.exists(every) else [])
    return _DOM[chrom]


def covered(chrom, pos):
    # Data.domainAt(): any block covering pos, looking back at most 100 kb (blocks are exon-sized)
    return any(b[0] <= pos <= b[1] and b[0] >= pos - 100000 for b in domain_blocks(chrom))


def expected(path, status):
    n = 0; f = dict.fromkeys(KEYS, 0)
    for t in records(path):
        I = info(t[7]); n += 1
        f['domain'] += covered(t[0], int(t[1]))
        f['gene'] += present(I.get('GENEMODEL')); f['effect'] += present(I.get('TYPE'))
        f['impact'] += 1                                  # pill always shown (MODIFIER default)
        f['mq'] += present(I.get('MQ')); f['comp'] += present(I.get('CVP')); f['r2'] += present(I.get('MAXR2'))
        f['maf'] += present(I.get('MAF'))
        f['qc'] += present(I.get('SITEQC')) or SIDECAR    # a class on the row (pass shown faintly)
        f['pc1'] += present(I.get('plantcad1_score')) or present(I.get('DNA_SCORE'))
        f['pc2op'] += present(I.get('plantcad2_onepass_score'))   # one-pass, not plantcad2_score
        f['evo2'] += present(I.get('evo2_score'))
        f['esmc'] += present(I.get('ESMC_score'))
        f['esm1'] += present(I.get('ESM1_score')) or present(I.get('AA_SCORE'))
        f['esm2'] += present(I.get('ESM2_score')); f['esm3'] += present(I.get('ESM3_score'))
    for k in KEYS:
        if status.get(k, 'ok') != 'ok':
            f[k] = 0                                      # unavailable / pending: empty cells
    return n, f


# chr2 query source: the demo store (fixture) or a full chr2 store (its window cut from the
# annotated chr2 release VCF on Ceres, release v1.4.3), whichever run_scenarios.js found installed
FULL_CHR2 = (json.load(open(res)).get('store') or {}).get('chr2', 0) > 100000
SETS = {'zmgrin2026_imp': f'{FX}/chr2_store/zmgrin2026_v1.4.3_chr2_4491424_4499434.annotated.vcf.gz' if FULL_CHR2
                          else f'{FX}/zmgrin2026_v1.4_chr2_testregions.vcf.gz',
        # MaizeGDB 2026, when offered again (run_scenarios.js renders them under the same switch):
        # 'mgdb2026_hq': f'{FX}/mgdb2026_hq_chr2_4491424_4499434.sites.vcf.gz',
        # 'mgdb2026_hc': f'{FX}/mgdb2026_hc_chr2_4491424_4499434.sites.vcf.gz',
        }
summary = {}
PREC = {}
for ds, path in SETS.items():
    a = A.get(ds)
    if not a:
        bad += 1; print('no render for', ds); continue
    hdr = a['headers']
    # a 'hidden' column (MQ and COMP for the GRIN-linked call set) is left out of the table
    hidden = {x['key'] for x in a.get('fields', []) if x['status'] == 'hidden'}
    # an 'auto' column (Site QC) is shown only when the result carries its values: the sidecar
    # gives them to every row, the release's own INFO to the rows that have SITEQC there
    auto = {x['key'] for x in a.get('fields', []) if x['status'] == 'auto'}
    if ds == 'zmgrin2026_imp' and auto != {'qc'}:
        bad += 1; print('auto columns', ds, sorted(auto))
    with_qc = SIDECAR or any(present(info(t[7]).get('SITEQC')) for t in records(path))
    shown = [k for k in KEYS if k not in hidden and (k not in auto or with_qc)]
    if hdr[:4] != ['CHR', 'POS', 'REF', 'ALT'] or hdr[4:] != [LABELS[KEYS.index(k)] for k in shown]:
        bad += 1; print('columns differ', ds, hdr)
    if a.get('cellsPerRow') != a.get('headerCells'):
        bad += 1; print('row cells != header cells', ds, a.get('cellsPerRow'), a.get('headerCells'))
    if ds == 'zmgrin2026_imp' and hidden != {'mq', 'comp'}:
        bad += 1; print('hidden columns', ds, sorted(hidden))
    status = {c['key']: c['status'] for c in a['cols']}
    n, f = expected(path, status)
    if a['rows'] != n:
        bad += 1; print('rows', ds, a['rows'], n)
    for k in shown:
        if a['filled'][k] != f[k]:
            bad += 1; print('filled', ds, k, a['filled'][k], 'expected', f[k])
    # numeric precision as rendered never exceeds the INFO's (PlantCAD1/2 and Evo2: 0.1 in the rebuilt stores)
    for k, info_key in (('pc1', 'plantcad1_score'), ('pc2op', 'plantcad2_onepass_score'), ('evo2', 'evo2_score'), ('esmc', 'ESMC_score'), ('r2', 'MAXR2')):
        dec = max([len(t.split('.')[1]) if '.' in t else 0 for t in (info(x[7]).get(info_key) for x in records(path)) if t not in (None, '', '.')] or [0])
        if a.get('decimals', {}).get(k, 0) > dec:
            bad += 1; print('precision', ds, k, a['decimals'][k], '>', dec)
        if ds == 'zmgrin2026_imp':
            PREC[k] = {'info_max_decimals': dec, 'rendered_max_decimals': a.get('decimals', {}).get(k), 'rendered_examples': a.get('sample', {}).get(k)}
    if any(v != 'ok' for v in status.values()) != bool(a['note']):
        bad += 1; print('availability note', ds, repr(a['note'][:80]))
    if ds == 'zmgrin2026_imp' and FULL_CHR2:
        # (d) the release window: the one-pass score on every SNP, on no indel
        ns = ni = s_with = i_with = 0
        for t in records(path):
            w = present(info(t[7]).get('plantcad2_onepass_score'))
            if len(t[3]) == 1 and len(t[4]) == 1: ns += 1; s_with += w
            else: ni += 1; i_with += w
        if not (ns and s_with == ns and i_with == 0 and a['filled']['pc2op'] == ns):
            bad += 1; print('PlantCAD2 one-pass', ds, f'SNPs {s_with}/{ns}, indels {i_with}/{ni}, cells {a["filled"]["pc2op"]}')
        else:
            print(f'PlantCAD2 column = plantcad2_onepass_score: {s_with} of {ns} SNPs, {i_with} of {ni} indels, {a["filled"]["pc2op"]} cells')
    summary[ds] = {k: ('hidden' if k in hidden else 'absent' if k not in status else status[k] if status[k] != 'ok' else f"{a['filled'][k]}/{n}") for k in KEYS}

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
print('chr2 store:', 'full chromosome' if FULL_CHR2 else 'demo fixture')
print('annotation columns:', json.dumps(summary, ensure_ascii=False))
print('zmgrin2026 precision:', json.dumps(PREC, ensure_ascii=False))
print(f'GRIN-linked INFO: {nsite} sites (SnpEff fields, MAF), {nesm} missense sites with ESM1/2/3; mismatches={bad}')
sys.exit(1 if bad else 0)
