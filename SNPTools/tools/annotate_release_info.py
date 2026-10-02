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

Not written, on purpose:
  MQ, CVC, CVP, MAXR2            the Schnable scored VCFs carry them as "." and no per-site
                                 value exists for the Grzybowski call set; values from the
                                 MaizeGDB 2026 call set describe different reads/samples.
  plantcad1_score/plantcad2_score  pending the PlantCAD Atlas merge. DNA_SCORE/AA_SCORE of the
                                 scored VCF are dropped so parseVcf does not fall back to them.

Both VCF inputs must be sorted by position within each chromosome (release files are);
the join streams, so whole-chromosome files work in constant memory apart from the ESM
table (~416k rows).
"""
import argparse, gzip, sys, collections

KEEP_FROM_SNPEFF = ('TYPE', 'EFFECT', 'GENEMODEL', 'SUB')
ADDED = KEEP_FROM_SNPEFF + ('MAF', 'ESM1_score', 'ESM2_score', 'ESM3_score')
HEADER = [
    '##INFO=<ID=TYPE,Number=.,Type=String,Description="SnpEff 5.2a consequence(s) (Sequence Ontology), from the MaizeGDB Schnable scored VCF of the Grzybowski et al. 2023 sites">',
    '##INFO=<ID=EFFECT,Number=.,Type=String,Description="SnpEff putative impact per consequence (HIGH/MODERATE/LOW/MODIFIER), same source">',
    '##INFO=<ID=GENEMODEL,Number=.,Type=String,Description="B73 v5 gene model(s) per consequence, same source">',
    '##INFO=<ID=SUB,Number=.,Type=String,Description="Amino-acid substitution per consequence (missense), same source">',
    '##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency over the called genotypes of the release samples in this file">',
    '##INFO=<ID=ESM1_score,Number=1,Type=Float,Description="ESM-1b 650M WT-marginal LLR (grz2023_missense_esm llr_esm1b), 1 decimal; missense only">',
    '##INFO=<ID=ESM2_score,Number=1,Type=Float,Description="ESM-2 650M WT-marginal LLR from the Full_ESM_stack store (esm2_store_score; not the published esm2_score), 1 decimal; missense only">',
    '##INFO=<ID=ESM3_score,Number=1,Type=Float,Description="ESM3-open sequence-track WT-marginal LLR (grz2023_missense_esm llr_esm3), 1 decimal; missense only">',
]


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


def load_esm(path):
    esm = collections.defaultdict(list)
    with opener(path) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in enumerate(hdr)}
        for c in ('chrom', 'pos', 'ref', 'alt', 'vcf_protein', 'variant', 'status', 'llr_esm1b', 'llr_esm2', 'llr_esm3'):
            if c not in ix:
                sys.exit(f'--esm: column {c} missing')
        for line in fh:
            t = line.rstrip('\n').split('\t')
            if t[ix['status']] != 'OK':
                continue
            esm[(t[ix['chrom']], int(t[ix['pos']]), t[ix['ref']], t[ix['alt']])].append(
                (t[ix['vcf_protein']], t[ix['variant']], t[ix['llr_esm1b']], t[ix['llr_esm2']], t[ix['llr_esm3']]))
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
    a = ap.parse_args()
    esm = load_esm(a.esm)
    stats = collections.Counter()
    se = snpeff_stream(a.snpeff)
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
                for k, v in (('ESM1_score', r[2]), ('ESM2_score', r[3]), ('ESM3_score', r[4])):
                    v = r1(v)
                    if v is not None:
                        info[k] = v
            t[7] = fmt_info(info)
            out.write('\t'.join(t) + '\n')
    print(' '.join(f'{k}={v}' for k, v in sorted(stats.items())), file=sys.stderr)


if __name__ == '__main__':
    main()
