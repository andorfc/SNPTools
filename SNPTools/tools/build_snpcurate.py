#!/usr/bin/env python3
"""build_snpcurate.py -- SNPCurate's data file from the hand-edited registry of published alleles.

    python3 tools/build_snpcurate.py                         # writes js/snpcurate.data.js
    python3 tools/build_snpcurate.py --table localdev/out/snpcurate_benchmark.tsv
                                                             # also the benchmark table (.tsv + .md)

Inputs (all read-only): data/curate/snpcurate.source.json (the registry; its fields are described in
that file's "fields" block), the variant stores hdf5/version3/zmgrin2026_<chr>_impute.h5 and their
site-QC sidecars (tools/build_site_qc.py), js/zmgrin.catalog.js (the 933 lines: subpopulation,
country) and data/traits/zmgrin2026.traits.json (GRIN trait summaries per line), gff/genes_index.txt
(gene intervals).

For an entry with status "site" it finds exactly one site with that chr, pos, ref and alt (else it
stops) and computes: heterozygous / homozygous-alternate / missing calls over all lines, the Site QC
class, the alternate-allele frequency among called lines, MAF and MAXR2; the consequence and
substitution for site_gene if given, else for gene (that gene's most severe SnpEff entry, as the app
reads it: Data.rowForGene), with the seven scores (ESM scores only on the entry they were computed
for) and the priority as Data.impactPriority computes it (combined = mean of PlantCAD1 and ESM1 when
both exist); the homozygous and heterozygous carriers, homozygous carriers per subpopulation and
carriers per country with each group's total; the trait by genotype when grin_descriptor is set
(numeric: mean and n of the line means for homozygous-reference and homozygous-alternate lines; coded:
per code, homozygous-alternate lines of the lines with that code); and, for a missense entry, its rank
by each score (1 = most negative) among the gene's missense variants that have a carrier.
For any other entry it confirms that no site exists at its pos or codon_positions, and stops if one does.

Marks: gold = causal, validated, a site and nothing to review; outline = any other site; grey = not a site.
Needs h5py + numpy.
"""
import argparse, datetime, json, os, re, sys

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIXED = ('CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'INFO')
SEVERITY = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
SCORES = (('pc1', 'plantcad1_score'), ('pc2', 'plantcad2_score'), ('evo2', 'evo2_score'),
          ('esm1', 'ESM1_score'), ('esm2', 'ESM2_score'), ('esm3', 'ESM3_score'), ('esmc', 'ESMC_score'))
ESM_KEYS = ('esm1', 'esm2', 'esm3', 'esmc')
MARKS = {'gold': 'Validated causal change, genotyped in this release.',
         'outline': 'Published marker, associated change or tagging site, genotyped in this release.',
         'grey': 'Known allele that this release cannot show.'}


class BuildError(SystemExit):
    pass


def site_qc(het, hom):          # the rule of tools/build_site_qc.py, one site at a time
    if het + hom == 0: return 'NO_CARRIER'
    if hom == 0: return 'HET_ONLY'
    if het > hom: return 'HET_EXCESS'
    if het * 4 >= het + hom: return 'HET_ELEVATED'
    return 'PASS'


def impact_class(effect):       # Data.impactClass in js/data.js
    e = re.sub(r'\s+', '_', str(effect or '').lower())
    for pat, klass, label, severe in (
            (r'missense|protein_altering|non[_-]?synonymous', 'missense', 'Missense', False),
            (r'stop_gained|nonsense', 'lof', 'Stop gained', True), (r'frameshift', 'lof', 'Frameshift', True),
            (r'start_lost|initiator', 'lof', 'Start lost', True), (r'stop_lost', 'lof', 'Stop lost', True),
            (r'splice_(acceptor|donor)', 'splice', 'Splice site', True), (r'splice', 'splice', 'Splice region', False),
            (r'inframe_insertion', 'indel', 'In-frame insertion', False), (r'inframe_deletion', 'indel', 'In-frame deletion', False),
            (r'synonymous|stop_retained', 'syn', 'Synonymous', False), (r'intron', 'other', 'Intron', False),
            (r'5_prime_utr|five_prime', 'other', "5' UTR", False), (r'3_prime_utr|three_prime', 'other', "3' UTR", False),
            (r'upstream', 'other', 'Upstream', False), (r'downstream', 'other', 'Downstream', False),
            (r'intergenic', 'other', 'Intergenic', False)):
        if re.search(pat, e):
            return {'klass': klass, 'label': label, 'severe': severe}
    return {'klass': 'other', 'label': clean_tok(effect) or 'Other', 'severe': False}


def impact_priority(cls, combined, level):   # Data.impactPriority
    if cls['severe']: return 'TOP'
    if cls['klass'] in ('missense', 'indel'):
        if combined is not None:
            if combined <= -7: return 'TOP'
            if combined <= -4: return 'HIGH'
            if combined <= -1: return 'MODERATE'
            return 'LOW'
        return level if level in ('HIGH', 'MODERATE') else 'LOW'
    return 'HIGH' if level == 'HIGH' else 'LOW'


def clean_tok(s):               # Data's cleanTok
    t = re.split(r'[;,]+', str(s or ''))[0].replace('_', ' ').strip()
    return '' if t == '.' else t


def parse_info(raw):
    s = raw.decode() if isinstance(raw, bytes) else str(raw)
    return dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in s.split(';') if kv and s != '.')


def num(v):
    if v in (None, '', '.'): return None
    try: return float(str(v).split(',')[0])
    except ValueError: return None


def annotations(info):
    """Data.annotationsOf: one entry per SnpEff consequence (gene, effect, impact, sub, esm_here).
    A single-consequence site reads its most severe EFFECT token, as Data.parseVcf does."""
    G, T = info.get('GENEMODEL', ''), info.get('TYPE', '')
    if ',' not in G:
        best = max((t.strip().upper() for t in re.split(r'[;,]+', info.get('EFFECT', '')) if t.strip().upper() in SEVERITY),
                   key=lambda t: SEVERITY[t], default='MODIFIER')
        return [{'gene': clean_tok(G) or '—', 'effect': clean_tok(T) or 'intergenic', 'impact': best,
                 'sub': clean_tok(info.get('SUB')), 'esm': True}]
    G, T, E, U = (info.get(k, '').split(',') for k in ('GENEMODEL', 'TYPE', 'EFFECT', 'SUB'))
    if len(T) != len(G):
        return annotations({k: v.split(',')[0] for k, v in info.items()})
    esm_at = next((i for i, t in enumerate(T) if 'missense' in t), -1)
    return [{'gene': clean_tok(G[i]) or '—', 'effect': clean_tok(T[i]) or 'intergenic',
             'impact': (E[i].strip().upper() if i < len(E) and E[i].strip().upper() in SEVERITY else 'MODIFIER'),
             'sub': clean_tok(U[i]) if i < len(U) else '', 'esm': i == esm_at} for i in range(len(G))]


def entry_for(info, gene):      # Data.rowForGene: the gene's most severe entry, the first on a tie
    best = None
    for a in annotations(info):
        if a['gene'] == gene and (best is None or SEVERITY[a['impact']] > SEVERITY[best['impact']]):
            best = a
    return best


def scores_of(info, entry):
    out = {k: num(info.get(key)) for k, key in SCORES}
    if not entry['esm']:
        for k in ESM_KEYS: out[k] = None
    return out


def combined_of(sc):
    pc, es = sc['pc1'], sc['esm1']
    if pc is not None and es is not None: return round((pc + es) / 2, 2)
    return pc if pc is not None else es


class Store:
    """One chromosome's store, read-only, with its sidecar (checked against the store)."""
    def __init__(self, root, chrom):
        self.path = os.path.join(root, 'hdf5', 'version3', f'zmgrin2026_{chrom}_impute.h5')
        if not os.path.isfile(self.path): raise BuildError(f'no store for {chrom}: {self.path}')
        self.f = h5py.File(self.path, 'r')
        self.pos = self.f['POS'][:].astype(np.int64)
        self.samples = [k for k in self.f.keys() if k not in FIXED]
        side = self.path[:-3] + '.siteqc.h5'
        if not os.path.isfile(side): raise BuildError(f'no site-QC sidecar beside {self.path}: run make site-qc')
        self.s = h5py.File(side, 'r')
        if int(self.s.attrs['n_sites']) != len(self.pos) or int(self.s.attrs['store_bytes']) != os.path.getsize(self.path):
            raise BuildError(f'{side} does not match its store: run make site-qc')

    def find(self, pos, ref, alt):
        lo, hi = np.searchsorted(self.pos, pos, 'left'), np.searchsorted(self.pos, pos, 'right')
        return [int(i) for i in range(lo, hi) if self.f['REF'][i].decode() == ref and self.f['ALT'][i].decode() == alt]

    def any_at(self, pos):
        return int(np.searchsorted(self.pos, pos, 'right') - np.searchsorted(self.pos, pos, 'left'))

    def genotypes(self, i):
        return {k: int(self.f[k][i]) for k in self.samples}

    def window(self, lo, hi):
        a, b = int(np.searchsorted(self.pos, lo, 'left')), int(np.searchsorted(self.pos, hi, 'right'))
        return a, b, self.s['NHET'][a:b], self.s['NHOM'][a:b], self.f['INFO'][a:b]


def load_catalog(root):
    s = open(os.path.join(root, 'js', 'zmgrin.catalog.js'), encoding='utf-8').read()
    fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', s, re.S).group(1))
    return {a['id']: a for p in fam['projects'] for g in p['groups'] for a in g['accessions']}


def gene_index(root):
    out = {}
    for line in open(os.path.join(root, 'gff', 'genes_index.txt')):
        if line.startswith('#'): continue
        t = [x.strip() for x in line.split('\t')]
        if len(t) >= 4: out[t[0]] = (t[1], int(t[2]), int(t[3]))
    return out


def trait_by_genotype(desc, gt, traits):
    d = traits['dictionary'].get(desc)
    if not d: raise BuildError(f'descriptor {desc} is not in data/traits/zmgrin2026.traits.json')
    out = {'descriptor': desc, 'name': d.get('name') or desc, 'scale': d.get('scale'), 'unit': d.get('unit') or '', 'url': d.get('url')}
    val = lambda sid: (traits['samples'].get(sid) or {}).get('traits', {}).get(desc)
    if d.get('scale') == 'numeric':
        for key, code in (('ref', 0), ('alt', 2)):
            xs = [val(sid)[1] for sid, g in gt.items() if g == code and val(sid) and val(sid)[1] is not None]
            out[key] = {'n': len(xs), 'mean': round(float(np.mean(xs)), 4) if xs else None}
    else:
        rows = []
        for code, label in (d.get('codes') or {}).items():
            lines = [sid for sid in gt if val(sid) and str(val(sid)[1]) == str(code)]
            if lines:
                rows.append({'code': code, 'label': label, 'lines': len(lines), 'alt': sum(gt[s] == 2 for s in lines)})
        out['codes'] = rows
    return out


def build(root):
    src = json.load(open(os.path.join(root, 'data', 'curate', 'snpcurate.source.json'), encoding='utf-8'))
    cat = load_catalog(root)
    traits = json.load(open(os.path.join(root, 'data', 'traits', 'zmgrin2026.traits.json'), encoding='utf-8'))
    genes = gene_index(root)
    stores = {}
    store = lambda c: stores.setdefault(c, Store(root, c))
    sub_tot, cty_tot = {}, {}
    for a in cat.values():
        sub_tot[a.get('subpop') or 'Unknown'] = sub_tot.get(a.get('subpop') or 'Unknown', 0) + 1
        cty_tot[a.get('country') or 'Unknown'] = cty_tot.get(a.get('country') or 'Unknown', 0) + 1
    out = []
    for e in src['entries']:
        r = dict(e)
        if e['status'] != 'site':
            r['mark'] = 'grey'
            st = store(e['chr'])
            for p in ([e['pos']] if e.get('pos') else []) + list(e.get('codon_positions') or []):
                if st.any_at(p):
                    raise BuildError(f"{e['id']}: a site exists at {e['chr']}:{p}, but the entry says {e['status']}")
            out.append(r); continue
        r['mark'] = 'gold' if (e['kind'] == 'causal' and e['evidence'] == 'validated' and not e.get('needs_review')) else 'outline'
        st = store(e['chr'])
        hit = st.find(e['pos'], e['ref'], e['alt'])
        if len(hit) != 1:
            raise BuildError(f"{e['id']}: {len(hit)} sites for {e['chr']}:{e['pos']} {e['ref']}>{e['alt']} (need exactly one)")
        i = hit[0]
        gt = st.genotypes(i)
        if set(gt) != set(cat): raise BuildError(f"{e['id']}: the store's samples differ from the catalogue")
        het = sum(g == 1 for g in gt.values()); hom = sum(g == 2 for g in gt.values()); miss = sum(g == 3 for g in gt.values())
        if (het, hom) != (int(st.s['NHET'][i]), int(st.s['NHOM'][i])):
            raise BuildError(f"{e['id']}: genotype counts differ from the sidecar")
        info = parse_info(st.f['INFO'][i])
        target = e.get('site_gene') or e['gene']
        ent = entry_for(info, target)
        first = annotations(info)[0]
        ann = ent or first
        cls = impact_class(ann['effect'])
        sc = scores_of(info, ann)
        comb = combined_of(sc)
        called = len(gt) - miss
        homs = sorted(s for s, g in gt.items() if g == 2); hets = sorted(s for s, g in gt.items() if g == 1)
        subpop = {}
        for s in homs: subpop[cat[s].get('subpop') or 'Unknown'] = subpop.get(cat[s].get('subpop') or 'Unknown', 0) + 1
        country = {}
        for s in homs + hets: country[cat[s].get('country') or 'Unknown'] = country.get(cat[s].get('country') or 'Unknown', 0) + 1
        site = {
            'nHet': het, 'nHom': hom, 'nMiss': miss, 'qc': site_qc(het, hom),
            'af': round((het + 2 * hom) / (2 * called), 4) if called else None,
            'maf': num(info.get('MAF')), 'maxr2': num(info.get('MAXR2')),
            'gene': target, 'geneEntry': bool(ent), 'effect': ann['effect'], 'consequence': cls['label'], 'consClass': cls['klass'],
            'impact': ann['impact'], 'sub': ann['sub'] or None, 'scores': sc, 'combined': comb,
            'priority': impact_priority(cls, comb, ann['impact']),
            'carriers': {'hom': homs, 'het': hets},
            'subpop': sorted(({'name': k, 'hom': subpop.get(k, 0), 'total': n} for k, n in sub_tot.items()),
                             key=lambda x: (-x['hom'], -x['total'], x['name'])),
            'country': sorted(({'name': k, 'carriers': v, 'total': cty_tot[k]} for k, v in country.items()),
                              key=lambda x: (-x['carriers'], x['name'])),
        }
        if e.get('grin_descriptor'):
            site['trait'] = trait_by_genotype(e['grin_descriptor'], gt, traits)
        if cls['klass'] == 'missense' and ent:
            if target not in genes: raise BuildError(f"{e['id']}: gene {target} is not in gff/genes_index.txt")
            gc, glo, ghi = genes[target]
            a, b, nh, nm, infos = st.window(glo, ghi)
            pool = []
            for k in range(b - a):
                if int(nh[k]) + int(nm[k]) == 0: continue
                I = parse_info(infos[k]); x = entry_for(I, target)
                if x and impact_class(x['effect'])['klass'] == 'missense':
                    pool.append(scores_of(I, x))
            mine = sc
            site['rank'] = {k: ([1 + sum(1 for p in pool if p[k] is not None and p[k] < mine[k]), sum(1 for p in pool if p[k] is not None)]
                                if mine[k] is not None else None) for k, _ in SCORES}
            site['rankPool'] = len(pool)
        exp = e.get('expect')
        if exp and (exp['het'], exp['hom'], exp['site_qc']) != (het, hom, site['qc']):
            raise BuildError(f"{e['id']}: het/hom/class {het}/{hom}/{site['qc']} differ from its expect block {exp}")
        r['site'] = site
        out.append(r)
    return src, out


def write_js(root, src, entries):
    today = datetime.date.today().isoformat()
    obj = {'release': src.get('release'), 'reference': src.get('reference'), 'built': today,
           'source': 'data/curate/snpcurate.source.json', 'marks': MARKS, 'entries': entries}
    p = os.path.join(root, 'js', 'snpcurate.data.js')
    with open(p, 'w', encoding='utf-8') as fh:
        fh.write(f'/* Generated by tools/build_snpcurate.py on {today} from data/curate/snpcurate.source.json, the\n'
                 f'   variant stores and their site-QC sidecars, js/zmgrin.catalog.js and data/traits. Do not edit by\n'
                 f'   hand: edit the source and run `make curate`. {len(entries)} entries. */\n')
        fh.write('window.SNP_CURATE = ' + json.dumps(obj, ensure_ascii=False, separators=(',', ':')) + ';\n')
    return p


def write_table(path, entries):
    cols = ['id', 'gene', 'symbol', 'change', 'reference', 'status', 'consequence', 'ESM2', 'PlantCAD2', 'priority',
            'het/hom', 'site_qc', 'mark', 'literature_match']
    fmt = lambda v: '' if v is None else (f'{v:.1f}' if isinstance(v, float) else str(v))
    rows = []
    for e in entries:
        s = e.get('site') or {}
        rows.append([e['id'], e['gene'], e['symbol'], e['label'], '; '.join(r['cite'] for r in e.get('refs', [])),
                     e['status'], s.get('consequence') or e.get('variant_type') or '', fmt((s.get('scores') or {}).get('esm2')),
                     fmt((s.get('scores') or {}).get('pc2')), s.get('priority') or '',
                     f"{s['nHet']}/{s['nHom']}" if s else '', s.get('qc') or '', e['mark'], e.get('note') or ''])
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows: fh.write('\t'.join(x.replace('\t', ' ') for x in r) + '\n')
    md = re.sub(r'\.tsv$', '', path) + '.md'
    with open(md, 'w', encoding='utf-8') as fh:
        fh.write('| ' + ' | '.join(c.replace('_', ' ') for c in cols) + ' |\n|' + '---|' * len(cols) + '\n')
        for r in rows: fh.write('| ' + ' | '.join(x.replace('|', '/') for x in r) + ' |\n')
    return path, md


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--root', default=ROOT, help='SNPTools web root (default: this tree)')
    ap.add_argument('--table', help='also write the benchmark table to this .tsv (and the same name .md)')
    ap.add_argument('--no-js', action='store_true', help='do not write js/snpcurate.data.js')
    a = ap.parse_args()
    src, entries = build(a.root)
    marks = {m: sum(e['mark'] == m for e in entries) for m in ('gold', 'outline', 'grey')}
    sites = sum(e['status'] == 'site' for e in entries)
    for e in entries:
        s = e.get('site')
        print(f"{e['id']} {e['symbol']:<16} {e['mark']:<7} {e['status']:<13}" +
              (f" {s['nHet']:>3}/{s['nHom']:<3} {s['qc']:<12} {s['consequence']:<18} {s['priority']}" if s else ''))
    print(f'{len(entries)} entries: {sites} sites, {len(entries) - sites} not sites; marks {marks}')
    if not a.no_js:
        print('wrote', os.path.relpath(write_js(a.root, src, entries), a.root))
    if a.table:
        print('wrote', *write_table(a.table, entries))


if __name__ == '__main__':
    main()
