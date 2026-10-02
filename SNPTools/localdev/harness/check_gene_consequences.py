#!/usr/bin/env python3
"""check_gene_consequences.py <results.json> <site_root> -- one gene's consequences at sites
SnpEff annotates in several genes (GRIN-linked 2026), recomputed from the fixture VCFs the
local stores are built from, independently of the JS.

SnpEff writes GENEMODEL, TYPE, EFFECT and SUB as parallel comma lists, one entry per gene
(most severe first). A tool showing one gene must use that gene's entry (its most severe,
the first on a tie), not the row's first entry:
(a) SNPFunction (Data.geneFunction): every site with an entry for the gene, its consequence
    class, residue substitution, ESM1/2/3, and homozygous / heterozygous carriers among the
    933 lines; the burden by class and the damaging-allele count.
(b) SNPFold (Data.queryFoldVariants): the coding entries with a residue.
(c) ESM: a site carries one ESM score, computed for the FIRST missense entry of TYPE
    (tools/annotate_release_info.py, pick_esm); it may appear only on that entry, and every
    score shown for gene G and substitution S must equal llr_esm1b of the ESM table row for
    (G, S) -- fixtures/annotation/grz2023_missense_esm_testregions.tsv.gz.
(d) SNPGeo gene search: each site in the gene interval reads the gene's entry where it has
    one, else the row's first entry.
(e) SNPVersity: "+N" counts the other distinct entries of a site; the region's gene list
    names every single gene model with an entry, not only first-listed ones.
Standard library only."""
import csv, glob, gzip, io, json, re, sys

res, root = sys.argv[1:3]
P = json.load(open(res))['pergene']
FX = f'{root}/localdev/fixtures'
SEV = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
bad = 0


def clean(tok):            # Data's cleanTok on one list entry
    t = (tok or '').split(';')[0].replace('_', ' ').strip()
    return '' if t == '.' else t


def klass(effect):         # Data.impactClass, by effect text
    e = re.sub(r'\s+', '_', (effect or '').lower())
    if re.search(r'missense|protein_altering|non[_-]?synonymous', e): return 'missense', False
    if re.search(r'stop_gained|nonsense|frameshift|start_lost|initiator|stop_lost', e): return 'lof', True
    if re.search(r'splice_(acceptor|donor)', e): return 'splice', True
    if 'splice' in e: return 'splice', False
    if re.search(r'inframe_(insertion|deletion)', e): return 'indel', False
    if re.search(r'synonymous|stop_retained', e): return 'syn', False
    return 'other', False


def fold_coding(effect):   # Data.classifyConsequence is not null
    e = re.sub(r'\s+', '_', (effect or '').lower())
    return bool(re.search(r'missense|protein_altering|non[_-]?synonymous|stop_gained|nonsense|frameshift|'
                          r'start_lost|initiator_codon|stop_lost|inframe_insertion|inframe_deletion', e))


def dose(cell):            # Data's gtCode: 0/1/2, None = missing
    g = cell.split(':')[0].replace('|', '/')
    a = g.split('/')
    if len(a) < 2 or '.' in a or '' in a: return None
    return (a[0] != '0') + (a[1] != '0')


def num(v):
    return None if v in (None, '', '.') else float(v.split(',')[0])


ROWS = {}                  # chrom -> [(pos, ref, alt, info dict, genotype cells)]
for f in glob.glob(f'{FX}/zmgrin2026_v1.4_chr*_testregions.vcf.gz'):
    for line in gzip.open(f, 'rt'):
        if line[0] == '#': continue
        t = line.rstrip('\n').split('\t')
        I = dict(kv.split('=', 1) for kv in t[7].split(';') if '=' in kv)
        ROWS.setdefault(t[0], []).append((int(t[1]), t[3], t[4], I, t[9:]))
ESM = {}
for r in csv.DictReader(io.TextIOWrapper(gzip.open(f'{FX}/annotation/grz2023_missense_esm_testregions.tsv.gz')), delimiter='\t'):
    ESM[(r['chrom'], int(r['pos']), r['ref'], r['alt'], r['vcf_protein'], r['variant'])] = r


def entries(I):
    """[(gene, effect, impact, sub, esm_here)]: one per consequence. A single-consequence site
    keeps its ESM whatever its TYPE (Data reads such a row as is)."""
    G, T = I.get('GENEMODEL', '').split(','), I.get('TYPE', '').split(',')
    E, U = I.get('EFFECT', '').split(','), I.get('SUB', '').split(',')
    esm_at = next((i for i, x in enumerate(T) if 'missense' in x), -1)
    multi = len(G) > 1
    return [(clean(G[i]) or '—', clean(T[i]) or 'intergenic',
             E[i].strip().upper() if i < len(E) and E[i].strip().upper() in SEV else 'MODIFIER',
             clean(U[i]) if i < len(U) else '', (i == esm_at) if multi else True) for i in range(len(G))]


def for_gene(I, gene):
    best = None
    for e in entries(I):
        if e[0] == gene and (best is None or SEV[e[2]] > SEV[best[2]]): best = e
    return best


def scores(I, e):
    return tuple(num(I.get(k)) if e[4] else None for k in ('ESM1_score', 'ESM2_score', 'ESM3_score'))


def residue(sub):
    m = re.search(r'\d+', sub or '')
    return int(m.group()) if m else None


summary = []
for gene, got in P['genes'].items():
    iv = got['interval']; chrom, lo, hi = iv['chr'], int(iv['start']), int(iv['end'])
    rows = [r for r in ROWS.get(chrom, []) if lo <= r[0] <= hi]
    exp_v, exp_fold, first_rule = [], [], 0
    for pos, ref, alt, I, cells in rows:
        if clean(I.get('GENEMODEL', '').split(',')[0]) == gene: first_rule += 1
        e = for_gene(I, gene)
        if not e: continue
        cls, severe = klass(e[1]); esm = scores(I, e)
        d = [dose(c) for c in cells]
        hom, het = sum(x == 2 for x in d), sum(x == 1 for x in d)
        exp_v.append({'pos': pos, 'ref': ref, 'alt': alt, 'cls': cls, 'esm': esm, 'hom': hom, 'het': het,
                      'severe': severe, 'sub': e[3]})
        if fold_coding(e[1]) and residue(e[3]) is not None:
            exp_fold.append({'pos': pos, 'resi': residue(e[3]), 'esm': esm[0]})
        # (c) a score shown for this gene must be the ESM table's score for this gene's substitution
        if esm[0] is not None and cls == 'missense':
            r = ESM.get((chrom, pos, ref, alt, gene, e[3]))
            if r is None or f'{round(float(r["llr_esm1b"]), 1):.1f}' != f'{esm[0]:.1f}':
                bad += 1; print('ESM misattributed', gene, pos, e[3], esm[0], r and r['llr_esm1b'])
    fn = got['fn']
    gv = {(v['pos'], v['ref'], v['alt']): v for v in fn['v']}
    if fn['n'] != len(exp_v) or set(gv) != {(v['pos'], v['ref'], v['alt']) for v in exp_v}:
        bad += 1; print('SNPFunction sites', gene, fn['n'], len(exp_v))
    by = {}
    for v in exp_v: by[v['cls']] = by.get(v['cls'], 0) + 1
    if {k: n for k, n in fn['byClass'].items() if n} != by:
        bad += 1; print('burden by class', gene, fn['byClass'], by)
    for v in exp_v:
        g = gv.get((v['pos'], v['ref'], v['alt']))
        if not g: continue
        if g['cls'] != v['cls'] or (g['hom'], g['het']) != (v['hom'], v['het']) or \
           (g['esm'], g['esm2'], g['esm3']) != v['esm']:
            bad += 1; print('SNPFunction variant', gene, v['pos'], g, v)
        if v['cls'] == 'missense' and g['sub'] != v['sub']:
            bad += 1; print('substitution', gene, v['pos'], g['sub'], v['sub'])
    dmg = sum(v['severe'] or v['cls'] == 'lof' or (v['cls'] == 'missense' and v['esm'][0] is not None and v['esm'][0] <= -4) for v in exp_v)
    if fn['damaging'] != dmg: bad += 1; print('damaging', gene, fn['damaging'], dmg)
    if fn['nAcc'] != 933: bad += 1; print('panel', gene, fn['nAcc'])
    key = lambda v: (v['pos'], v['resi'], 'none' if v['esm'] is None else f"{v['esm']:.1f}")   # JSON writes -3.0 as -3
    fold = sorted(key(v) for v in got['fold'])
    if fold != sorted(key(v) for v in exp_fold):
        bad += 1; print('SNPFold variants', gene, len(fold), len(exp_fold))
    summary.append(f"{gene}: {len(exp_v)} sites (first-listed rule kept {first_rule}), "
                   f"{by.get('missense', 0)} missense ({sum(v['esm'][0] is not None for v in exp_v if v['cls'] == 'missense')} with ESM), "
                   f"{len(exp_fold)} on the structure, {dmg} damaging")

# (d) SNPGeo, gene search Zm00001eb374230
iv = P['genes']['Zm00001eb374230']['interval']
exp = {}
for pos, ref, alt, I, _ in ROWS[iv['chr']]:
    if int(iv['start']) <= pos <= int(iv['end']):
        e = for_gene(I, 'Zm00001eb374230') or entries(I)[0]
        exp[pos] = (e[0], e[1], e[3], scores(I, e)[0])
geo = {r['pos']: (r['gene'], r['effect'], r['sub'], r['esm1']) for r in P['geo']}
if geo != exp:
    bad += 1; print('SNPGeo rows differ', [p for p in exp if geo.get(p) != exp[p]][:5])
if geo.get(13120567, ())[:3] != ('Zm00001eb374230', 'missense variant', 'A471G'):
    bad += 1; print('SNPGeo 13120567', geo.get(13120567))

# (e) SNPVersity window: +N markers and the gene list
V = P['versity']; chrom, lo, hi = V['window']
rows = [r for r in ROWS[chrom] if lo <= r[0] <= hi]
marks, genes = [], set()
for pos, ref, alt, I, _ in sorted(rows, key=lambda r: r[0]):
    E = entries(I)
    genes |= {e[0] for e in E if e[0] != '—' and not re.search(r'\s|\.\.', e[0])}
    if len(E) < 2: marks.append(''); continue
    seen, n = set(), 0
    for i, e in enumerate(E):
        k = (e[0], e[1], e[3])
        if k in seen: continue
        seen.add(k); n += (i > 0)
    marks.append(f'+{n}' if n else '')
if V['rows'] != len(rows) or V['markers'] != marks:
    bad += 1; print('SNPVersity +N markers', V['rows'], len(rows), sum(bool(m) for m in V['markers']), sum(bool(m) for m in marks))
if sorted(V['genes']) != sorted(genes):
    bad += 1; print('SNPVersity gene list', V['genes'], sorted(genes))

for s in summary: print(s)
print(f"SNPGeo Zm00001eb374230: {len(geo)} sites, 13120567 reads {geo.get(13120567)}; "
      f"SNPVersity {chrom}:{lo}-{hi}: {sum(bool(m) for m in marks)} of {len(rows)} sites marked +N, genes {', '.join(sorted(genes))}; mismatches={bad}")
sys.exit(1 if bad else 0)
