#!/usr/bin/env python3
"""annotate_release_info.py -- add the per-variant annotation SNPVersity shows to a
GRIN-linked release VCF (dataset zmgrin2026_imp) before it is converted to HDF5.

    python3 tools/annotate_release_info.py --vcf in.vcf.gz --snpeff scored_sites.vcf.gz \
        --esm grz2023_missense_esm.tsv.gz --out out.vcf.gz

Writes INFO keys in the names js/data.js parseVcf() reads (the same names the
schnable2023 and MaizeGDB 2026 stores carry):

  TYPE, EFFECT, GENEMODEL, SUB   SnpEff 5.2a consequence lists, copied unchanged from the
                                 MaizeGDB "Schnable scored" VCFs (Box Data_with_dna_aa_scores/
                                 schnable, chr<N>_schnable_scored.vcf.gz; the Grzybowski et al.
                                 2023 sites). Matched on CHROM/POS/REF/ALT.
  MAF                            minor-allele frequency over the called genotypes of THIS
                                 file's samples (the 933 release samples), 4 decimals.
  ESM1_score, ESM2_score,        missense only, from grz2023_missense_esm.tsv(.gz): llr_esm1b,
  ESM3_score                     llr_esm2 (= the Full_ESM_stack store ESM-2 650M layer, i.e.
                                 esm2_store_score; NOT the published esm2_score), llr_esm3;
                                 rounded to 1 decimal like the published ESM tables. For a site
                                 with several missense consequences the value belongs to the
                                 first missense consequence of the TYPE list (the gene/SUB
                                 SNPVersity displays).

  plantcad1_score, plantcad2_score, evo2_score   optional (--dna-scores): the grz2023 Atlas
                                 PlantCAD1/PlantCAD2 (all SNPs) and Evo2 (genic +/-1 kb subset)
                                 tables, joined on CHROM/POS/REF/ALT; indels get none.
  ESMC_score                     llr_esmc from the same missense table, 1 decimal.

  MAXR2                          optional (--maxr2): highest pairwise LD r2 (PLINK 1.9 --r2) of the
                                 site with any variant 400-5,000 bp away, computed from the 933
                                 release v1.4 genotypes with no filtering (tools: maxr2_chr.sh);
                                 written as up to 6 decimals like the MaizeGDB 2026 stores.

Not written, on purpose:
  MQ, CVC, CVP                   the Schnable scored VCFs carry them as "." and no per-site
                                 value exists for the Grzybowski call set; values from the
                                 MaizeGDB 2026 call set describe different reads/samples.
                                 SNPTools shows them as NA.
  DNA_SCORE/AA_SCORE              of the scored VCF are dropped so parseVcf does not fall back
                                 to them (PlantCAD comes only from --dna-scores).

Both VCF inputs must be sorted by position within each chromosome (release files are);
the join streams, so whole-chromosome files work in constant memory apart from the ESM
table (~416k rows).
"""
import argparse, gzip, sys, collections

KEEP_FROM_SNPEFF = ('TYPE', 'EFFECT', 'GENEMODEL', 'SUB')
DNA_KEYS = ('plantcad1_score', 'plantcad2_score', 'evo2_score')
ADDED = KEEP_FROM_SNPEFF + ('MAF', 'ESM1_score', 'ESM2_score', 'ESM3_score', 'ESMC_score') + DNA_KEYS + ('MAXR2',)
HEADER = [
    '##INFO=<ID=TYPE,Number=.,Type=String,Description="SnpEff 5.2a consequence(s) (Sequence Ontology), from the MaizeGDB Schnable scored VCF of the Grzybowski et al. 2023 sites">',
    '##INFO=<ID=EFFECT,Number=.,Type=String,Description="SnpEff putative impact per consequence (HIGH/MODERATE/LOW/MODIFIER), same source">',
    '##INFO=<ID=GENEMODEL,Number=.,Type=String,Description="B73 v5 gene model(s) per consequence, same source">',
    '##INFO=<ID=SUB,Number=.,Type=String,Description="Amino-acid substitution per consequence (missense), same source">',
    '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency over the called genotypes of the release samples in this file">',
    '##INFO=<ID=ESM1_score,Number=1,Type=Float,Description="ESM-1b 650M WT-marginal LLR (grz2023_missense_esm llr_esm1b), 1 decimal; missense only">',
    '##INFO=<ID=ESM2_score,Number=1,Type=Float,Description="ESM-2 650M WT-marginal LLR from the Full_ESM_stack store (esm2_store_score; not the published esm2_score), 1 decimal; missense only">',
    '##INFO=<ID=ESM3_score,Number=1,Type=Float,Description="ESM3-open sequence-track WT-marginal LLR (grz2023_missense_esm llr_esm3), 1 decimal; missense only">',
    '##INFO=<ID=ESMC_score,Number=1,Type=Float,Description="ESM C 600M WT-marginal LLR (grz2023_missense_esm llr_esmc), 1 decimal; missense only">',
]

DNA_HEADER = [
    '##INFO=<ID=plantcad1_score,Number=1,Type=Float,Description="PlantCAD1 (PlantCaduceus) zero-shot score, grz2023 Atlas re-score (--dna-scores), rounded to 1 decimal (--pc-decimals); SNPs only">',
    '##INFO=<ID=plantcad2_score,Number=1,Type=Float,Description="PlantCAD2 zero-shot score, grz2023 Atlas re-score (--dna-scores), rounded to 1 decimal (--pc-decimals); SNPs only">',
    '##INFO=<ID=evo2_score,Number=1,Type=Float,Description="Evo2 7B log-likelihood ratio (256-bp left context), genic +/-1 kb SNP subset only (--dna-scores), rounded to 1 decimal (--pc-decimals)">',
]
MAXR2_HEADER = '##INFO=<ID=MAXR2,Number=1,Type=Float,Description="Highest pairwise LD r2 (PLINK 1.9 --r2) with any variant 400-5,000 bp away, from the 933 release v1.4 genotypes, no MAF/missingness filtering (--maxr2)">'


def fmt_r2(v):
    s = f'{float(v):.6f}'.rstrip('0')
    return s + '0' if s.endswith('.') else s


def opener(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p)


def parse_info(s):
    out = collections.OrderedDict()
    if s in ('.', ''):
        return out
    for kv in s.split(';'):
        k, _, v = kv.partition('=')
        out[k] = v if _ else None
    return out


def fmt_info(d):
    if not d:
        return '.'
    return ';'.join(k if v is None else f'{k}={v}' for k, v in d.items())


def snpeff_stream(path):
    """yield (chrom, pos, ref, alt, {TYPE,EFFECT,GENEMODEL,SUB}) in file order."""
    with opener(path) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            t = line.rstrip('\n').split('\t', 8)
            info = parse_info(t[7])
            yield t[0], int(t[1]), t[3], t[4], {k: info[k] for k in KEEP_FROM_SNPEFF if info.get(k) not in (None, '.', '')}


def dna_stream(path, keys_wanted=DNA_KEYS, opt='--dna-scores'):
    """yield (chrom, pos, ref, alt, {plantcad1_score, plantcad2_score, evo2_score}) from a
    position-sorted TSV with header chr,pos,ref,alt,<any of DNA_KEYS>; empty/nan values skipped."""
    with opener(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for c in ('chr', 'pos', 'ref', 'alt'):
            if c not in ix:
                sys.exit(f'{opt}: column {c} missing')
        keys = [k for k in keys_wanted if k in ix]
        if not keys:
            sys.exit(f'{opt}: none of {keys_wanted} present')
        for line in fh:
            t = line.rstrip('\n').split('\t')
            yield t[ix['chr']], int(t[ix['pos']]), t[ix['ref']], t[ix['alt']], \
                {k: t[ix[k]] for k in keys if t[ix[k]] not in ('', 'nan', 'NA', '.')}


class SortedJoin:
    """Advance a position-sorted record stream alongside a position-sorted VCF."""
    def __init__(self, it):
        self.it = it; self.cur = next(it, None); self.key = None; self.buf = {}

    def get(self, chrom, pos, ref, alt):
        if self.key != (chrom, pos):
            while self.cur is not None and self.cur[0] != chrom:
                self.cur = next(self.it, None)
            while self.cur is not None and self.cur[0] == chrom and self.cur[1] < pos:
                self.cur = next(self.it, None)
            self.key, self.buf = (chrom, pos), {}
            while self.cur is not None and self.cur[0] == chrom and self.cur[1] == pos:
                self.buf[(self.cur[2], self.cur[3])] = self.cur[4]
                self.cur = next(self.it, None)
        return self.buf.get((ref, alt), {})

def load_esm(path):
    esm = collections.defaultdict(list)
    with opener(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for c in ('chrom', 'pos', 'ref', 'alt', 'vcf_protein', 'variant', 'status', 'llr_esm1b', 'llr_esm2', 'llr_esm3', 'llr_esmc'):
            if c not in ix:
                sys.exit(f'--esm: column {c} missing')
        for line in fh:
            t = line.rstrip('\n').split('\t')
            if t[ix['status']] != 'OK':
                continue
            esm[(t[ix['chrom']], int(t[ix['pos']]), t[ix['ref']], t[ix['alt']])].append(
                (t[ix['vcf_protein']], t[ix['variant']], t[ix['llr_esm1b']], t[ix['llr_esm2']], t[ix['llr_esm3']], t[ix['llr_esmc']]))
    return esm


def r1(v):
    try:
        return f'{round(float(v), 1):.1f}'
    except ValueError:
        return None


def pick_esm(rows, ann, stats):
    """ESM row for the first missense consequence of TYPE (the gene/SUB shown)."""
    types = (ann.get('TYPE') or '').split(',')
    genes = (ann.get('GENEMODEL') or '').split(',')
    subs = (ann.get('SUB') or '').split(',')
    for i, ty in enumerate(types):
        if 'missense' in ty:
            g = genes[i] if i < len(genes) else None
            s = subs[i] if i < len(subs) else None
            for r in rows:
                if r[0] == g and r[1] == s:
                    return r
            stats['esm_no_exact_consequence_match'] += 1
            return rows[0]
    stats['esm_without_missense_type'] += 1
    return rows[0]


def maf_of(cells):
    alt = n = 0
    for c in cells:
        g = c.split(':', 1)[0]
        if len(g) < 3 or g[0] == '.' or g[2] == '.':
            continue
        a, b = g[0], g[2]
        n += 2
        alt += (a != '0') + (b != '0')
    if not n:
        return None
    p = alt / n
    return f'{min(p, 1 - p):.4f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--snpeff', required=True, help='Schnable scored VCF (or its sites-only copy), same sites')
    ap.add_argument('--esm', required=True, help='grz2023_missense_esm.tsv(.gz)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--pc-decimals', type=int, default=1, help='decimals for plantcad1/2_score and evo2_score (default 1; full precision stays in the score tables)')
    ap.add_argument('--dna-scores', help='position-sorted TSV chr,pos,ref,alt,plantcad1_score,plantcad2_score,evo2_score (optional)')
    ap.add_argument('--maxr2', help='position-sorted TSV chr,pos,ref,alt,MAXR2 (optional; maxr2_chr.sh output)')
    a = ap.parse_args()
    esm = load_esm(a.esm)
    stats = collections.Counter()
    se = snpeff_stream(a.snpeff)
    dj = SortedJoin(dna_stream(a.dna_scores)) if a.dna_scores else None
    mj = SortedJoin(dna_stream(a.maxr2, ('MAXR2',), '--maxr2')) if a.maxr2 else None
    cur = next(se, None)
    posbuf_key, posbuf = None, {}
    out = gzip.open(a.out, 'wt', compresslevel=6) if a.out.endswith('.gz') else open(a.out, 'w')
    with opener(a.vcf) as fh, out:
        for line in fh:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=') and line[11:].split(',', 1)[0] in ADDED:
                    continue
                out.write(line)
                continue
            if line.startswith('#CHROM'):
                out.write('\n'.join(HEADER) + '\n')
                if dj is not None:
                    out.write('\n'.join(DNA_HEADER) + '\n')
                if mj is not None:
                    out.write(MAXR2_HEADER + '\n')
                out.write('##annotate_release_info=TYPE/EFFECT/GENEMODEL/SUB from ' + a.snpeff.split('/')[-1] +
                          '; ESM from ' + a.esm.split('/')[-1] + '; MAF from this file\'s genotypes\n')
                out.write(line)
                continue
            t = line.rstrip('\n').split('\t')
            chrom, pos, ref, alt = t[0], int(t[1]), t[3], t[4]
            stats['sites'] += 1
            if posbuf_key != (chrom, pos):
                # advance the sorted SnpEff stream to this chromosome/position and buffer every
                # record at this position (split multi-allelics may come in any ALT order)
                while cur is not None and cur[0] != chrom:
                    cur = next(se, None)
                while cur is not None and cur[0] == chrom and cur[1] < pos:
                    cur = next(se, None)
                posbuf_key, posbuf = (chrom, pos), {}
                while cur is not None and cur[0] == chrom and cur[1] == pos:
                    posbuf[(cur[2], cur[3])] = cur[4]
                    cur = next(se, None)
            ann = posbuf.get((ref, alt), {})
            stats['snpeff_matched' if ann else 'snpeff_unmatched'] += 1
            info = parse_info(t[7])
            for k in ADDED + ('DNA_SCORE', 'AA_SCORE'):
                info.pop(k, None)
            for k in KEEP_FROM_SNPEFF:
                if k in ann:
                    info[k] = ann[k]
            if dj is not None:
                dv = dj.get(chrom, pos, ref, alt)
                stats['dna_matched' if dv else 'dna_unmatched'] += 1
                for k in DNA_KEYS:
                    if k in dv:
                        v = dv[k]
                        if k in DNA_KEYS:
                            v = f'{round(float(v), a.pc_decimals):.{a.pc_decimals}f}'
                            if v in ('-0.0', '-0'):
                                v = v[1:]
                        info[k] = v; stats[k] += 1
            if mj is not None:
                mv = mj.get(chrom, pos, ref, alt).get('MAXR2')
                if mv is not None:
                    info['MAXR2'] = fmt_r2(mv); stats['MAXR2'] += 1
                else:
                    stats['maxr2_missing'] += 1
            m = maf_of(t[9:])
            if m is not None:
                info['MAF'] = m
            else:
                stats['maf_no_calls'] += 1
            rows = esm.get((chrom, pos, ref, alt))
            if rows:
                r = pick_esm(rows, ann, stats)
                stats['esm_sites'] += 1
                if len(rows) > 1:
                    stats['esm_sites_multi_consequence'] += 1
                for k, v in (('ESM1_score', r[2]), ('ESM2_score', r[3]), ('ESM3_score', r[4]), ('ESMC_score', r[5])):
                    v = r1(v)
                    if v is not None:
                        info[k] = v
            t[7] = fmt_info(info)
            out.write('\t'.join(t) + '\n')
    print(' '.join(f'{k}={v}' for k, v in sorted(stats.items())), file=sys.stderr)


if __name__ == '__main__':
    main()
