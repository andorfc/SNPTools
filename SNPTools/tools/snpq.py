#!/usr/bin/env python3
"""snpq.py -- the query library behind the SNPMaize MCP server (mcp.php).

One call = one process, like processForm.php -> h5_to_vcf.py: mcp.php validates the JSON-RPC
request, then runs

    python tools/snpq.py            (request JSON on stdin: {"tool": ..., "args": {...}, "site": "https://.../"})

and relays the JSON answer. Every tool reads the same files the pages read -- the variant stores
and their site-QC sidecars (hdf5/grin2026), js/zmgrin.catalog.js (passport), data/traits
(GRIN descriptors), js/zmgrin.lineqc.js, js/snpcurate.data.js (the SNPCurate registry as built),
data/qc/zmgrin2026.scorecal.json (composite calibration), gff/genes_index.txt, data/gwas and
distance/zmgrin2026 -- through the mirrors of the page logic in tools/build_snpcurate.py (site QC
rule, consequence classes, composite z and priority tiers, annotation parsing), so an answer here
agrees with what the same question shows on the page.

Answers are bounded (CAPS below) and every one carries `links`: SNPMaize permalinks (#tool/thing
hash routes) that open the same view on the site. Errors are {"ok": false, "error": {"code",
"message"}}; the codes are invalid_params, not_found, too_large and unavailable.

    python tools/snpq.py --list            the tool names
    python tools/snpq.py <tool> '<json>'   one call from the command line
"""
import bisect, csv, glob, json, os, re, signal, sys, time
from collections import Counter, OrderedDict
T_START = time.time()

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np                                                   # noqa: E402
import h5py                                                          # noqa: E402
from build_snpcurate import (site_qc, impact_class, impact_priority, annotations, entry_for,   # noqa: E402
                             scores_of, composite_of, parse_info, num, clean_tok, SEVERITY,
                             FIXED, GT_WORD, B73, load_catalog, gene_index, trait_by_genotype)

FAMILY = 'zmgrin2026'
CHROMS = ['chr%d' % i for i in range(1, 11)]
DATASETS = OrderedDict([
    ('zmgrin2026_imp', {'suffix': '_impute', 'name': 'MaizeGDB GRIN-linked 2026 · Imputed (Grzybowski 2023 sites)',
                        'note': 'Beagle 5-imputed calls at the Grzybowski et al. 2023 sites; every line has a call at every site.'}),
    ('zmgrin2026_hq', {'suffix': '_HQ', 'name': 'MaizeGDB GRIN-linked 2026 · High Quality',
                       'note': 'The same lines at the sites that pass the High Quality rule (>= 1 carrier and, at sites with >= 10 carriers, LD max R² >= 0.2).'}),
])
CAPS = {'max_span_bp': 250000, 'max_flank_bp': 50000, 'max_rows': 500, 'default_rows': 100,
        'max_lines': 100, 'max_sites': 50, 'max_names': 200, 'max_list': 200, 'max_top': 100}
CLASSES = ('missense', 'lof', 'splice', 'indel', 'syn', 'other')
TIERS = ('TOP', 'HIGH', 'MODERATE', 'LOW')
QC_SETS = {'usable': ('PASS', 'HET_ELEVATED'), 'pass': ('PASS',),
           'all': ('PASS', 'HET_ELEVATED', 'HET_ONLY', 'HET_EXCESS', 'NO_CARRIER')}
PANELS = {'nam': 'inNAM', 'ames282': 'inAmes282', 'widiv': 'inWiDiv'}
PASSPORT = ('id', 'label', 'strain', 'grin', 'genotypeId', 'run', 'namFounder', 'panel', 'inNAM', 'inAmes282', 'inWiDiv',
            'country', 'state', 'originType', 'improvement', 'pedigree', 'subpop', 'era', 'species', 'kernelType',
            'cobColor', 'kernelColor', 'nTraits')
MAIZEGDB_GENE = 'https://www.maizegdb.org/gene_center/gene/'
MAIZEGDB_MCP = 'https://claude.maizegdb.org/mcp'
GRIN_SEARCH = 'https://npgsweb.ars-grin.gov/gringlobal/search?q='


class QueryError(Exception):
    def __init__(self, message, code='invalid_params', **extra):
        super().__init__(message)
        self.code, self.extra = code, extra


# ---------------------------------------------------------------- loaders (each file once per process)
_C = {}


def cached(name, fn):
    if name not in _C: _C[name] = fn()
    return _C[name]


def js_object(path, assign_re):
    s = open(os.path.join(ROOT, path), encoding='utf-8').read()
    m = re.search(assign_re, s)
    if not m: raise QueryError(f'{path} does not carry the expected assignment', 'unavailable')
    return json.JSONDecoder().raw_decode(s, m.end())[0]


def catalog():       return cached('catalog', lambda: load_catalog(ROOT))
def family():        return cached('family', lambda: js_object('js/zmgrin.catalog.js', r'families\["zmgrin2026"\]\s*=\s*'))
def traits():        return cached('traits', lambda: json.load(open(os.path.join(ROOT, 'data', 'traits', 'zmgrin2026.traits.json'), encoding='utf-8')))
def line_qc():       return cached('lineqc', lambda: js_object('js/zmgrin.lineqc.js', r'window\.SNP_LINE_QC\["zmgrin2026"\]\s*=\s*'))
def curate():        return cached('curate', lambda: js_object('js/snpcurate.data.js', r'window\.SNP_CURATE\s*=\s*'))
def scorecal():      return cached('scorecal', lambda: json.load(open(os.path.join(ROOT, 'data', 'qc', 'zmgrin2026.scorecal.json'), encoding='utf-8')))
def genes():         return cached('genes', lambda: gene_index(ROOT))
def siteqc_summary():
    p = os.path.join(ROOT, 'data', 'qc', 'zmgrin2026.siteqc.summary.json')
    return cached('siteqc', lambda: json.load(open(p, encoding='utf-8')) if os.path.isfile(p) else None)
def gwas_manifest(): return cached('gwas', lambda: json.load(open(os.path.join(ROOT, 'data', 'gwas', 'manifest.json'), encoding='utf-8')))


def norm_key(s):   return str('' if s is None else s).replace('﻿', '').strip().strip('"\'').lower()
def loose_key(s):  return re.sub(r'[^a-z0-9]+', '', norm_key(s))


def acc_index():      # buildAccIndex in js/snpversity.js: id, run, founder, label, GRIN accession, line name, genotype id
    def build():
        idx = {}
        def add(k, i):
            if k: idx.setdefault(k, set()).add(i)
        for a in catalog().values():
            keys = [a.get(k) for k in ('id', 'run', 'founder', 'label', 'grin', 'strain', 'genotypeId')]
            if a.get('founder') and a.get('run'):
                keys += [a['founder'] + '_' + a['run'], a['run'] + '_' + a['founder']]
            for v in keys:
                if v: add(norm_key(v), a['id']); add(loose_key(v), a['id'])
        return idx
    return cached('accidx', build)


def resolve_line(raw):
    """One name -> {query, ids, via} (ids: the accession ids it names) or None."""
    line = str(raw or '').replace('﻿', '').strip()
    if not line: return None
    idx = acc_index()
    for t in [line] + [x.strip().strip('"\'') for x in re.split(r'[,\t;|]+', line)]:
        hit = idx.get(norm_key(t)) or idx.get(loose_key(t))
        if hit: return {'query': line, 'ids': sorted(hit), 'via': t}
    return {'query': line, 'ids': [], 'via': None}


def resolve_lines(names, cap=CAPS['max_lines']):
    names = [str(n) for n in as_list(names) if str(n).strip()]
    if not names: raise QueryError('lines: give at least one line name, GRIN accession (PI 550473) or id (ZmG_B73).')
    if len(names) > cap: raise QueryError(f'lines: at most {cap} names per call ({len(names)} given).', 'too_large')
    ids, seen, unresolved = [], set(), []
    for n in names:
        r = resolve_line(n)
        if not r: continue
        if not r['ids']: unresolved.append(n); continue
        for i in r['ids']:
            if i not in seen: seen.add(i); ids.append(i)
    return ids, unresolved


def one_line(name):
    r = resolve_line(name)
    if not r or not r['ids']:
        cands = search_lines(str(name or ''))[:8]
        hint = ('; did you mean ' + ', '.join(f"{a['label']} ({a['id']})" for a in cands) + '?') if cands else ''
        raise QueryError(f'No line answers to "{name}"{hint}', 'not_found')
    if len(r['ids']) > 1:
        raise QueryError(f'"{name}" names {len(r["ids"])} lines: ' + ', '.join(r['ids'][:10]) + '. Give one id.', 'invalid_params')
    return catalog()[r['ids'][0]]


def search_lines(q, pedigree=False):
    """Lines whose id, name, GRIN accession or genotype id contain q (or whose pedigree does, with pedigree=True);
    the lines q names exactly come first."""
    q = norm_key(q)
    if not q: return []
    exact = (resolve_line(q) or {}).get('ids') or []
    out = []
    for a in catalog().values():
        keys = ('id', 'label', 'strain', 'grin', 'run', 'genotypeId') + (('pedigree',) if pedigree else ())
        hay = ' '.join(str(a.get(k) or '') for k in keys).lower()
        if q in hay or a['id'] in exact: out.append(a)
    return sorted(out, key=lambda a: (a['id'] not in exact, a.get('label') or a['id']))


def line_qc_of(sid):
    q = (line_qc().get('lines') or {}).get(sid)
    if not q: return None
    labels = {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'}
    return {'het_share': q[0], 'class': q[1], 'label': labels.get(q[1], q[1])}


def line_brief(a):
    return {'id': a['id'], 'label': a.get('label'), 'grin': a.get('grin'), 'subpop': a.get('subpop'), 'country': a.get('country'),
            'kernel_type': a.get('kernelType'), 'panel': a.get('panel'),
            'panels': [n for n, k in (('NAM', 'inNAM'), ('Ames282', 'inAmes282'), ('WiDiv', 'inWiDiv')) if a.get(k) == 'yes'],
            'nam_founder': a.get('namFounder'), 'line_qc': line_qc_of(a['id'])}


# ---------------------------------------------------------------- arguments
def as_list(v):
    if v is None: return []
    if isinstance(v, (list, tuple)): return list(v)
    return [x.strip() for x in re.split(r'[,\n;]+', str(v)) if x.strip()]


def arg_int(args, key, default=None, lo=None, hi=None):
    v = args.get(key, default)
    if v in (None, ''): return default
    try: v = int(v)
    except (TypeError, ValueError): raise QueryError(f'{key}: an integer is required, not {v!r}')
    if lo is not None and v < lo: raise QueryError(f'{key}: at least {lo}')
    if hi is not None and v > hi: raise QueryError(f'{key}: at most {hi}', 'too_large')
    return v


def arg_float(args, key, default=None):
    v = args.get(key, default)
    if v in (None, ''): return default
    try: return float(v)
    except (TypeError, ValueError): raise QueryError(f'{key}: a number is required, not {v!r}')


def arg_enum(args, key, allowed, default=None):
    v = args.get(key, default)
    if v in (None, ''): return default
    v = str(v).strip()
    for a in allowed:
        if v.lower() == a.lower(): return a
    raise QueryError(f'{key}: one of {", ".join(allowed)}')


def arg_dataset(args):
    return arg_enum(args, 'dataset', list(DATASETS), 'zmgrin2026_imp')


def parse_region(text):
    s = str(text or '').strip().replace(',', '').replace(' ', '')
    m = re.match(r'^(?:chr)?(\d{1,2})[:\-_.](\d+)(?:[\-_.]{1,2}|:)(\d+)$', s, re.I)
    if not m: raise QueryError(f'region: use chr:start-end on B73 v5, e.g. chr4:46640000-46660000 (got {text!r})')
    chrom, start, end = 'chr' + m.group(1), int(m.group(2)), int(m.group(3))
    if chrom not in CHROMS: raise QueryError(f'region: chromosome must be chr1..chr10 (got {chrom})')
    if start < 1 or end < start: raise QueryError('region: start must be >= 1 and end >= start')
    return chrom, start, end


def check_span(chrom, start, end, cap=CAPS['max_span_bp']):
    if end - start + 1 > cap:
        raise QueryError(f'{chrom}:{start}-{end} spans {end - start + 1:,} bp; at most {cap:,} bp per call. '
                         'Ask for a gene, or a smaller window, or page through the region.', 'too_large')


def resolve_gene(text):
    """Gene id (Zm00001eb175150, a transcript/protein id of it) or a registry symbol (tga1) -> (id, chrom, start, end)."""
    s = str(text or '').strip()
    if not s: raise QueryError('gene: a B73 v5 gene id (Zm00001eb175150) or a symbol from the SNPCurate registry (tga1) is required.')
    g = re.sub(r'_[TP]\d{3}$', '', s)
    gi = genes()
    if g in gi: return (g,) + gi[g]
    for k in gi:
        if k.lower() == g.lower(): return (k,) + gi[k]
    for e in curate()['entries']:
        if (e.get('symbol') or '').lower() == s.lower() and e.get('gene') in gi:
            return (e['gene'],) + gi[e['gene']]
    raise QueryError(f'gene {s!r} is not a B73 v5 (Zm00001eb) gene model in this release and not a registry symbol. '
                     'Resolve names to gene ids with the MaizeGDB MCP (maizegdb_resolve).', 'not_found')


def region_or_gene(args, need=True):
    """-> dict(chrom, start, end, gene?, flank, label)."""
    flank = arg_int(args, 'flank', 0, 0, CAPS['max_flank_bp'])
    if args.get('gene'):
        gid, chrom, gs, ge = resolve_gene(args['gene'])
        r = {'chrom': chrom, 'start': max(1, gs - flank), 'end': ge + flank, 'gene': gid, 'gene_start': gs, 'gene_end': ge, 'flank': flank}
    elif args.get('region'):
        chrom, s, e = parse_region(args['region'])
        r = {'chrom': chrom, 'start': s, 'end': e, 'gene': None, 'flank': 0}
    elif need:
        raise QueryError('Give a gene (id or registry symbol) or a region (chr:start-end on B73 v5).')
    else:
        return None
    check_span(r['chrom'], r['start'], r['end'])
    r['label'] = f"{r['chrom']}:{r['start']:,}-{r['end']:,}"
    return r


# ---------------------------------------------------------------- stores
def Store(chrom, dataset='zmgrin2026_imp'):
    """One open store per (chromosome, dataset) per process."""
    return cached(('store', chrom, dataset), lambda: _Store(chrom, dataset))


class _Store:
    """One chromosome of one dataset, read-only, with its site-QC sidecar (checked against the store as the builder does)."""
    def __init__(self, chrom, dataset='zmgrin2026_imp'):
        suf = DATASETS[dataset]['suffix']
        self.dataset, self.chrom = dataset, chrom
        self.path = os.path.join(ROOT, 'hdf5', 'grin2026', f'{FAMILY}_{chrom}{suf}.h5')
        if not os.path.isfile(self.path):
            raise QueryError(f'The {DATASETS[dataset]["name"]} store for {chrom} is not installed on this server.', 'unavailable')
        self.f = h5py.File(self.path, 'r')
        self.pos = self.f['POS'][:].astype(np.int64)
        self.samples = [k for k in self.f.keys() if k not in FIXED]
        side = self.path[:-3] + '.siteqc.h5'
        self.s = None
        if os.path.isfile(side):
            s = h5py.File(side, 'r')
            if int(s.attrs['n_sites']) == len(self.pos) and int(s.attrs['store_bytes']) == os.path.getsize(self.path):
                self.s = s

    def span(self, lo, hi):
        return int(np.searchsorted(self.pos, lo, 'left')), int(np.searchsorted(self.pos, hi, 'right'))

    def qc_arrays(self, a, b):
        if self.s is not None:
            return self.s['NHET'][a:b], self.s['NHOM'][a:b], self.s['NMISS'][a:b]
        # no usable sidecar: count from the genotypes (slow path, small windows only)
        het = np.zeros(b - a, dtype=np.int64); hom = np.zeros(b - a, dtype=np.int64); mis = np.zeros(b - a, dtype=np.int64)
        for k in self.samples:
            g = self.f[k][a:b]
            het += (g == 1); hom += (g == 2); mis += (g >= 3)
        return het, hom, mis

    def column(self, sid, a, b):
        return self.f[sid][a:b] if sid in self.f else None

    def at(self, pos, ref=None, alt=None):
        a, b = self.span(pos, pos)
        rows = [i for i in range(a, b) if (ref is None or self.f['REF'][i].decode() == ref) and (alt is None or self.f['ALT'][i].decode() == alt)]
        return rows

    def record(self, i):
        return {'chr': self.chrom, 'pos': int(self.pos[i]), 'ref': self.f['REF'][i].decode(), 'alt': self.f['ALT'][i].decode(),
                'info': parse_info(self.f['INFO'][i])}


def site_row(st, i, info, het, hom, mis, gene=None):
    """One variant as the pages describe it: consequence for `gene` (or the most severe), scores, composite, tier, QC."""
    cal = scorecal()
    ann = annotations(info)
    e = entry_for(info, gene) if gene else None
    if e is None:
        e = max(ann, key=lambda x: SEVERITY[x['impact']]) if ann else {'gene': '—', 'effect': 'intergenic', 'impact': 'MODIFIER', 'sub': '', 'esm': True}
    cls = impact_class(e['effect'])
    sc = scores_of(info, e)
    comb, n = composite_of(sc, cal)
    tier = impact_priority(cls, comb, e['impact'], cal['composite']['tiers'])
    pct = None
    if comb is not None:
        pct = int(min(100, bisect.bisect_left(cal['composite']['percentiles'], comb)))
    return {'chr': st.chrom, 'pos': int(st.pos[i]), 'ref': st.f['REF'][i].decode(), 'alt': st.f['ALT'][i].decode(),
            'gene': None if e['gene'] == '—' else e['gene'], 'consequence': cls['label'], 'class': cls['klass'],
            'impact': e['impact'], 'sub': e['sub'] or None, 'other_genes': sorted({x['gene'] for x in ann if x['gene'] not in (e['gene'], '—')}),
            'scores': {k: v for k, v in sc.items() if v is not None}, 'composite': comb, 'composite_n': n,
            'composite_percentile': pct, 'tier': tier,
            'maf': num(info.get('MAF')), 'maxr2': num(info.get('MAXR2')),
            'n_het': int(het), 'n_hom': int(hom), 'n_missing': int(mis), 'site_qc': site_qc(int(het), int(hom))}


def scan(st, r, gene=None):
    a, b = st.span(r['start'], r['end'])
    het, hom, mis = st.qc_arrays(a, b)
    infos = st.f['INFO'][a:b]
    for j, i in enumerate(range(a, b)):
        yield i, site_row(st, i, parse_info(infos[j]), het[j], hom[j], mis[j], gene)


# ---------------------------------------------------------------- links
def links_for(site, **routes):
    out = {}
    for k, h in routes.items():
        if h: out[k] = (site or '') + h
    return out


def route(tool, thing=None, query=None):
    from urllib.parse import quote
    return '#' + tool + ('/' + quote(str(thing), safe='') if thing else '') + ('?' + query if query else '')


def region_route(r, dataset):
    return route('snpversity', None, f"ds={dataset}&chr={r['chrom']}&start={r['start']}&end={r['end']}")


# ---------------------------------------------------------------- tools
def t_datasets(args, site):
    fam, sq, cal, cu, lq = family(), siteqc_summary(), scorecal(), curate(), line_qc()
    installed = {d: [c for c in CHROMS if os.path.isfile(os.path.join(ROOT, 'hdf5', 'grin2026', f'{FAMILY}_{c}{m["suffix"]}.h5'))]
                 for d, m in DATASETS.items()}
    marks = Counter(e.get('mark') for e in cu['entries'])
    out = {
        'release': {'family': FAMILY, 'label': fam.get('label'), 'reference': cu.get('reference') or 'B73 v5 (Zm-B73-REFERENCE-NAM-5.0)',
                    'lines': fam.get('count'), 'nam_founders': fam.get('namFounders'),
                    'sites': (sq or {}).get('total', {}).get('sites'), 'coordinates': '1-based, B73 v5; chromosomes chr1..chr10'},
        'datasets': [{'id': d, 'name': m['name'], 'note': m['note'], 'chromosomes_installed': installed[d], 'default': d == 'zmgrin2026_imp'}
                     for d, m in DATASETS.items()],
        'site_qc': {'rule': (sq or {}).get('rule'), 'codes': (sq or {}).get('codes'), 'usable': (sq or {}).get('usable'),
                    'flagged': (sq or {}).get('flagged'), 'counts': (sq or {}).get('total'),
                    'default_filter': 'usable (PASS or HET_ELEVATED); ask qc=all to see every site'},
        'line_qc': {'thresholds': lq.get('thresholds'), 'labels': {'I': 'Inbred', 'E': 'Elevated heterozygosity', 'H': 'Heterozygous sample'},
                    'counts': Counter(v[1] for v in lq.get('lines', {}).values()), 'share': 'heterozygous calls at clean sites / clean sites'},
        'scores': {'per_site': {'plantcad1_score': 'pc1', 'plantcad2_score': 'pc2', 'plantcad2_onepass_score': 'pc2op', 'evo2_score': 'evo2',
                                'ESM1_score': 'esm1', 'ESM2_score': 'esm2', 'ESM3_score': 'esm3', 'ESMC_score': 'esmc'},
                   'direction': 'more negative = more likely damaging; ESM scores are set only at missense sites',
                   'composite': {'of': cal['composite']['of'], 'definition': cal['composite']['definition'],
                                 'reference': f"{cal['n_reference']:,} {cal['reference']}", 'tiers': cal['composite']['tiers'],
                                 'tier_percentiles': cal['composite']['tier_percentiles'],
                                 'rule': 'severe consequences (stop gained, frameshift, start/stop lost, splice site) are TOP regardless of score; '
                                         'missense and in-frame indels take the tier of their composite; other consequences follow SnpEff impact'}},
        'registry': {'entries': len(cu['entries']), 'by_mark': dict(marks), 'marks': cu.get('marks'), 'built': cu.get('built')},
        'traits': {'descriptors': len(traits()['dictionary']), 'keys': sorted(traits()['dictionary']),
                   'lines_with_traits': sum(1 for s in traits()['samples'].values() if s.get('traits'))},
        'gwas': {'studies': len(gwas_manifest()), 'populations': sorted({g.get('population') for g in gwas_manifest()}),
                 'publications': sorted({g.get('publication') for g in gwas_manifest()})},
        'tools': ['snpmaize_datasets', 'snpmaize_lines', 'snpmaize_line', 'snpmaize_variants', 'snpmaize_gene', 'snpmaize_carriers',
                  'snpmaize_curated', 'snpmaize_alleles', 'snpmaize_similarity', 'snpmaize_gwas'],
        'caps': CAPS,
        'notes': ['Genotypes: 0/0 ref, 1/0 het, 1/1 alt; the Imputed set has no missing calls, the High Quality set may.',
                  'Line names resolve by id (ZmG_B73), line name (B73), GRIN accession (PI 550473) or genotype id (ZmG:B73).',
                  'Gene records, sequences, expression, GO, QTL and literature live on the MaizeGDB MCP at ' + MAIZEGDB_MCP
                  + '; use it to turn a gene name into a Zm00001eb id, then query SNPMaize by that id.'],
        'links': links_for(site, snpversity=route('snpversity'), snpcurate=route('snpcurate'), snptrait=route('snptrait')),
    }
    if site: out['links']['site'] = site
    return out


def trait_value(desc, v):
    d = traits()['dictionary'].get(desc) or {}
    if not v: return None
    out = {'descriptor': desc, 'name': d.get('name') or desc, 'category': d.get('category'), 'scale': d.get('scale'), 'unit': d.get('unit') or None, 'n_obs': v[0]}
    if d.get('scale') == 'numeric':
        out.update({'mean': v[1], 'sd': v[2] if len(v) > 2 else None, 'min': v[3] if len(v) > 3 else None, 'max': v[4] if len(v) > 4 else None})
    elif d.get('scale') == 'ordinal_or_nominal_coded':
        out.update({'code': v[1], 'label': (d.get('codes') or {}).get(str(v[1]), str(v[1]))})
    else:
        out.update({'text': v[1]})
    return out


def find_descriptor(text):
    s = str(text or '').strip()
    dic = traits()['dictionary']
    if s in dic: return s
    for k, d in dic.items():
        if s.lower() in (k.lower(), (d.get('name') or '').lower()): return k
    hits = [k for k, d in dic.items() if s.lower() in k.lower() or s.lower() in (d.get('name') or '').lower()]
    if len(hits) == 1: return hits[0]
    if hits: raise QueryError(f'trait {s!r} matches several descriptors: ' + ', '.join(hits[:12]), 'invalid_params')
    raise QueryError(f'trait {s!r} is not a GRIN descriptor of this release; snpmaize_datasets lists the keys.', 'not_found')


def descriptor_summary(desc):
    d = traits()['dictionary'][desc]
    vals = [(sid, s['traits'][desc]) for sid, s in traits()['samples'].items() if s.get('traits', {}).get(desc)]
    out = {'descriptor': desc, 'name': d.get('name'), 'category': d.get('category'), 'scale': d.get('scale'), 'unit': d.get('unit') or None,
           'url': d.get('url'), 'lines_measured': len(vals)}
    if d.get('scale') == 'numeric':
        xs = np.array([v[1] for _, v in vals if v[1] is not None], dtype=float)
        if len(xs):
            q = np.percentile(xs, [0, 25, 50, 75, 100])
            out.update({'mean': round(float(xs.mean()), 3), 'min': float(q[0]), 'q1': float(q[1]), 'median': float(q[2]), 'q3': float(q[3]), 'max': float(q[4])})
        fb = (traits().get('facetBins') or {}).get(desc)
        if fb: out['facet_bins'] = fb.get('labels')
    else:
        c = Counter(str(v[1]) for _, v in vals)
        codes = d.get('codes') or {}
        out['values'] = [{'code': k, 'label': codes.get(k, k), 'lines': n} for k, n in c.most_common()]
    return out


def t_lines(args, site):
    limit = arg_int(args, 'limit', 50, 1, CAPS['max_list'])
    q = str(args.get('query') or '').strip()
    desc = find_descriptor(args['trait']) if args.get('trait') else None
    tval = args.get('trait_value')
    tmin, tmax = arg_float(args, 'trait_min'), arg_float(args, 'trait_max')
    panel = arg_enum(args, 'panel', ['NAM', 'Ames282', 'WiDiv'])
    filters = {k: str(args[k]).strip() for k in ('subpop', 'country', 'kernel_type', 'improvement', 'era', 'species', 'state') if args.get(k)}
    key_of = {'kernel_type': 'kernelType'}
    pedigree = str(args.get('pedigree', '')).lower() in ('1', 'true', 'yes')
    rows = search_lines(q, pedigree) if q else sorted(catalog().values(), key=lambda a: a.get('label') or a['id'])
    samples = traits()['samples']
    out_rows = []
    for a in rows:
        ok = True
        for k, v in filters.items():
            if norm_key(a.get(key_of.get(k, k))) != norm_key(v) and norm_key(v) not in norm_key(a.get(key_of.get(k, k))): ok = False; break
        if not ok: continue
        if panel and a.get(PANELS[panel.lower()]) != 'yes': continue
        tv = None
        if desc:
            tv = (samples.get(a['id']) or {}).get('traits', {}).get(desc)
            if not tv: continue
            d = traits()['dictionary'][desc]
            if d.get('scale') == 'numeric':
                x = tv[1]
                if x is None or (tmin is not None and x < tmin) or (tmax is not None and x > tmax): continue
            elif tval is not None:
                lab = (d.get('codes') or {}).get(str(tv[1]), '')
                tl, ll = norm_key(tval), norm_key(lab)       # code, label, or a label the value abbreviates (Pop / Popcorn, Red / Red / uniform)
                if tl != norm_key(tv[1]) and tl != ll and not ll.startswith(tl + ' ') and not (len(ll) >= 3 and tl.startswith(ll)): continue
        b = line_brief(a)
        if desc: b['trait'] = trait_value(desc, tv)
        out_rows.append(b)
    facets = {k: Counter(r.get(k) for r in out_rows).most_common(12) for k in ('subpop', 'country', 'kernel_type', 'panel')}
    out = {'matched': len(out_rows), 'returned': min(limit, len(out_rows)), 'truncated': len(out_rows) > limit,
           'lines': out_rows[:limit], 'facets': {k: [{'value': v, 'lines': n} for v, n in c] for k, c in facets.items()},
           'links': links_for(site, snptrait=route('snptrait'), snpversity=route('snpversity'))}
    if desc:
        out['descriptor'] = descriptor_summary(desc)
        if desc.endswith('-BLUP') or 'BLUP' in desc: out['descriptor']['note'] = 'BLUP values from the source study, not raw GRIN observations.'
    if not any([q, desc, panel, filters]):
        out['note'] = 'No filter given: the first lines of the panel. Filter by query, subpop, country, panel, kernel_type, improvement, era, species or a trait.'
    return out


def t_line(args, site):
    a = one_line(args.get('line') or args.get('name') or args.get('id'))
    s = traits()['samples'].get(a['id']) or {}
    tr = [trait_value(k, v) for k, v in sorted((s.get('traits') or {}).items())]
    out = {'line': line_brief(a), 'passport': {k: a.get(k) for k in PASSPORT if a.get(k) not in (None, '')},
           'is_nam_founder': bool(a.get('namFounder')), 'traits': tr, 'n_traits': len(tr),
           'links': links_for(site, snptrait=route('snptrait'), snpcompare=route('snpcompare'), snpversity=route('snpversity'))}
    if a.get('grin'): out['links']['grin_search'] = GRIN_SEARCH + a['grin'].replace(' ', '%20')
    reg = [e for e in curate()['entries'] if a['id'] in ((e.get('site') or {}).get('carriers') or {}).get('hom', [])
           or a['id'] in ((e.get('site') or {}).get('carriers') or {}).get('het', [])]
    if reg:
        out['curated_alleles_carried'] = [{'id': e['id'], 'symbol': e.get('symbol'), 'gene': e.get('gene'), 'label': e.get('label'), 'trait': e.get('trait'),
                                           'genotype': 'alt' if a['id'] in e['site']['carriers']['hom'] else 'het', 'favorable': e.get('favorable'),
                                           'link': (site or '') + route('snpcurate', e['id'], 'lines=' + a['id'])} for e in reg]
    return out


def t_variants(args, site):
    r = region_or_gene(args)
    ds = arg_dataset(args)
    limit = arg_int(args, 'limit', CAPS['default_rows'], 1, CAPS['max_rows'])
    classes = [arg_enum({'c': c}, 'c', CLASSES) for c in as_list(args.get('consequence'))]
    tier = arg_enum(args, 'tier', TIERS)
    qc = arg_enum(args, 'qc', list(QC_SETS), 'usable')
    min_maf = arg_float(args, 'min_maf')
    sort = arg_enum(args, 'sort', ['composite', 'position', 'maf'], 'composite')
    gene_only = str(args.get('gene_only', 'true')).lower() not in ('0', 'false', 'no') if r['gene'] else False
    ids, unresolved = ([], [])
    if args.get('lines'): ids, unresolved = resolve_lines(args['lines'])
    st = Store(r['chrom'], ds)
    a, b = st.span(r['start'], r['end'])
    cols = {sid: st.column(sid, a, b) for sid in ids}
    b73 = st.column(B73, a, b)
    total, kept = 0, []
    by_class, by_qc, by_tier = Counter(), Counter(), Counter()
    tier_rank = {t: i for i, t in enumerate(TIERS)}
    for j, (i, row) in enumerate(scan(st, r, r['gene'])):
        total += 1
        by_qc[row['site_qc']] += 1
        if gene_only and row['gene'] != r['gene']: continue
        if row['site_qc'] not in QC_SETS[qc]: continue
        if classes and row['class'] not in classes: continue
        if tier and tier_rank[row['tier']] > tier_rank[tier]: continue
        if min_maf is not None and (row['maf'] is None or row['maf'] < min_maf): continue
        if ids:
            carry = {sid: int(cols[sid][j]) for sid in ids if cols[sid] is not None}
            alt_in = [sid for sid, g in carry.items() if g in (1, 2)]
            if not alt_in: continue
            row['lines'] = {'alt_carriers': alt_in, 'genotypes': {sid: GT_WORD.get(g, 'missing') for sid, g in carry.items()}}
        if b73 is not None: row['b73'] = GT_WORD.get(int(b73[j]), 'missing')
        by_class[row['class']] += 1; by_tier[row['tier']] += 1
        kept.append(row)
    if sort == 'position': kept.sort(key=lambda x: x['pos'])
    elif sort == 'maf': kept.sort(key=lambda x: -(x['maf'] or 0))
    else: kept.sort(key=lambda x: (x['composite'] if x['composite'] is not None else 9e9, -SEVERITY[x['impact']], x['pos']))
    out = {'dataset': ds, 'region': {k: r[k] for k in ('chrom', 'start', 'end', 'label', 'gene', 'flank') if k in r},
           'sites_in_window': total, 'matched': len(kept), 'returned': min(limit, len(kept)), 'truncated': len(kept) > limit,
           'filters': {'qc': qc, 'consequence': classes or None, 'tier': tier, 'min_maf': min_maf, 'gene_only': gene_only, 'lines': ids or None,
                       'unresolved_lines': unresolved or None},
           'by_site_qc_in_window': dict(by_qc), 'by_class': dict(by_class), 'by_tier': dict(by_tier),
           'sort': sort, 'variants': kept[:limit],
           'links': links_for(site, snpversity=region_route(r, ds),
                              snpfunction=route('snpfunction', r['gene']) if r['gene'] else None,
                              snpfold=route('snpfold', r['gene']) if r['gene'] else None,
                              maizegdb_gene=None)}
    if r['gene']: out['links']['maizegdb_gene'] = MAIZEGDB_GENE + r['gene']
    if st.s is None: out['note'] = 'Site QC counted from the genotypes (no sidecar for this store).'
    return out


def t_gene(args, site):
    gid, chrom, gs, ge = resolve_gene(args.get('gene'))
    ds = arg_dataset(args)
    r = {'chrom': chrom, 'start': gs, 'end': ge, 'gene': gid, 'label': f'{chrom}:{gs:,}-{ge:,}'}
    check_span(chrom, gs, ge, CAPS['max_span_bp'] * 4)
    st = Store(chrom, ds)
    by_class, by_qc, by_tier = Counter(), Counter(), Counter()
    top, total, in_gene = [], 0, 0
    for i, row in scan(st, r, gid):
        total += 1
        if row['gene'] != gid: continue
        in_gene += 1
        by_qc[row['site_qc']] += 1
        if row['site_qc'] not in QC_SETS['usable']: continue
        by_class[row['class']] += 1
        if row['class'] in ('missense', 'indel', 'lof', 'splice'): by_tier[row['tier']] += 1
        if row['tier'] in ('TOP', 'HIGH'): top.append(row)
    top.sort(key=lambda x: (0 if x['impact'] == 'HIGH' and x['tier'] == 'TOP' else 1, x['composite'] if x['composite'] is not None else 9e9))
    reg = [{'id': e['id'], 'label': e.get('label'), 'mark': e.get('mark'), 'role': e.get('role'), 'status': e.get('status'), 'trait': e.get('trait'),
            'favorable': e.get('favorable'), 'site': f"{e['chr']}:{e['pos']} {e['ref']}>{e['alt']}" if e.get('pos') else None,
            'link': (site or '') + route('snpcurate', e['id'])} for e in curate()['entries'] if e.get('gene') == gid]
    gw = gwas_hits({'chrom': chrom, 'start': max(1, gs - 10000), 'end': ge + 10000}, None, None, True, 5)
    out = {'gene': gid, 'symbol': next((e.get('symbol') for e in curate()['entries'] if e.get('gene') == gid and e.get('symbol')), None),
           'location': {'chr': chrom, 'start': gs, 'end': ge, 'length_bp': ge - gs + 1}, 'dataset': ds,
           'sites_in_gene': in_gene, 'sites_in_window_other_genes': total - in_gene,
           'by_site_qc': dict(by_qc), 'by_class_usable': dict(by_class), 'by_tier_coding_usable': dict(by_tier),
           'top_variants': top[:10], 'n_top_or_high': len(top),
           'curated_alleles': reg,
           'gwas_within_10kb': {'significant_hits': gw['matched'], 'by_study': gw['by_study'], 'top': gw['hits']},
           'protein_sites_available': os.path.isfile(os.path.join(ROOT, 'data', 'sites', 'by_gene', gid + '.tsv')),
           'links': links_for(site, snpfunction=route('snpfunction', gid), snpfold=route('snpfold', gid), snpgeo=route('snpgeo', gid),
                              paneffect=route('paneffect', gid), snpversity=region_route(r, ds))}
    out['links']['maizegdb_gene'] = MAIZEGDB_GENE + gid
    return out


def registry_site(chrom, pos, ref=None, alt=None):
    for e in curate()['entries']:
        if e.get('chr') == chrom and e.get('pos') == pos and (ref is None or e.get('ref') == ref) and (alt is None or e.get('alt') == alt):
            return e
    return None


def locate_site(args, ds):
    """site args -> (Store, index, registry entry or None)."""
    if args.get('curated'):
        e = next((x for x in curate()['entries'] if x['id'].lower() == str(args['curated']).lower()), None)
        if not e: raise QueryError(f'registry id {args["curated"]!r} not found (ids look like SC0001).', 'not_found')
        if not e.get('pos') or e.get('status') == 'absent' or e.get('mark') == 'white':
            raise QueryError(f'{e["id"]} ({e.get("label")}) is not a site of this marker set; there are no calls to show.', 'not_found')
        chrom, pos, ref, alt = e['chr'], int(e['pos']), e['ref'], e['alt']
    else:
        s = str(args.get('site') or '').strip().replace(',', '')
        m = re.match(r'^(?:chr)?(\d{1,2})[:\s_-](\d+)(?:[:\s_/>-]+([ACGTN]+)[:\s_/>-]+([ACGTN]+))?$', s, re.I)
        if not m: raise QueryError('site: chr:pos or chr:pos:ref:alt on B73 v5, e.g. chr4:46648374:C:G; or curated=SC0001.')
        chrom, pos = 'chr' + m.group(1), int(m.group(2))
        ref = (args.get('ref') or m.group(3) or '').upper() or None
        alt = (args.get('alt') or m.group(4) or '').upper() or None
        if chrom not in CHROMS: raise QueryError('site: chromosome must be chr1..chr10')
    st = Store(chrom, ds)
    rows = st.at(pos, ref, alt)
    if not rows:
        near = st.span(pos - 50, pos + 50)
        hint = ', '.join(f'{int(st.pos[i])}:{st.f["REF"][i].decode()}>{st.f["ALT"][i].decode()}' for i in range(near[0], min(near[1], near[0] + 6)))
        raise QueryError(f'No site at {chrom}:{pos}' + (f' {ref}>{alt}' if ref and alt else '') + f' in {ds}.' + (f' Nearby: {hint}' if hint else ''), 'not_found')
    if len(rows) > 1:
        alts = ', '.join(f'{st.f["REF"][i].decode()}>{st.f["ALT"][i].decode()}' for i in rows)
        raise QueryError(f'{len(rows)} records at {chrom}:{pos} ({alts}); give ref and alt.', 'invalid_params')
    i = rows[0]
    return st, i, registry_site(chrom, pos, st.f['REF'][i].decode(), st.f['ALT'][i].decode())


def t_carriers(args, site):
    ds = arg_dataset(args)
    limit = arg_int(args, 'limit', 50, 0, CAPS['max_names'])
    st, i, reg = locate_site(args, ds)
    info = parse_info(st.f['INFO'][i])
    gt = {sid: int(st.f[sid][i]) for sid in st.samples}
    het, hom, mis = sum(g == 1 for g in gt.values()), sum(g == 2 for g in gt.values()), sum(g >= 3 for g in gt.values())
    row = site_row(st, i, info, het, hom, mis, reg.get('gene') if reg else None)
    cat = catalog()
    def grp(key):
        c = {}
        for sid, g in gt.items():
            k = (cat.get(sid) or {}).get(key) or 'unknown'
            d = c.setdefault(k, {'lines': 0, 'hom': 0, 'het': 0, 'missing': 0})
            d['lines'] += 1
            if g == 2: d['hom'] += 1
            elif g == 1: d['het'] += 1
            elif g >= 3: d['missing'] += 1
        rows = [{key: k, **v, 'alt_allele_freq': round((2 * v['hom'] + v['het']) / (2 * max(1, v['lines'] - v['missing'])), 4)}
                for k, v in c.items()]
        return sorted(rows, key=lambda x: (-(x['hom'] + x['het']), -x['lines']))
    def names(code):
        ids = sorted(sid for sid, g in gt.items() if g == code)
        return {'n': len(ids), 'lines': [line_brief(cat[s]) if s in cat else {'id': s} for s in ids[:limit]], 'truncated': len(ids) > limit}
    out = {'dataset': ds, 'site': row,
           'genotypes': {'ref_hom': sum(g == 0 for g in gt.values()), 'het': het, 'alt_hom': hom, 'missing': mis, 'lines': len(gt),
                         'alt_allele_freq': round((2 * hom + het) / (2 * max(1, len(gt) - mis)), 4)},
           'b73': GT_WORD.get(gt.get(B73, 3), 'missing'),
           'carriers': {'alt_hom': names(2), 'het': names(1)},
           'by_subpop': grp('subpop'), 'by_country': grp('country')[:25],
           'links': links_for(site, snpgeo=route('snpgeo', row['gene']) if row['gene'] else None,
                              snpversity=region_route({'chrom': st.chrom, 'start': max(1, row['pos'] - 2000), 'end': row['pos'] + 2000}, ds))}
    carriers_qc = Counter((line_qc_of(sid) or {}).get('class') for sid, g in gt.items() if g in (1, 2))
    out['carrier_line_qc'] = dict(carriers_qc)
    if reg:
        out['curated'] = {'id': reg['id'], 'symbol': reg.get('symbol'), 'label': reg.get('label'), 'mark': reg.get('mark'), 'role': reg.get('role'),
                          'trait': reg.get('trait'), 'favorable': reg.get('favorable'), 'favorable_note': reg.get('favorable_note')}
        out['links']['snpcurate'] = (site or '') + route('snpcurate', reg['id'])
    if args.get('trait'):
        desc = find_descriptor(args['trait'])
        out['trait_by_genotype'] = trait_by_genotype(desc, gt, traits())
        out['trait_by_genotype']['note'] = 'Means of the GRIN descriptor by genotype class at this one site: a description of the panel, not an association test.'
    return out


def t_curated(args, site):
    cu = curate()
    if args.get('id'):
        e = next((x for x in cu['entries'] if x['id'].lower() == str(args['id']).lower()), None)
        if not e: raise QueryError(f'registry id {args["id"]!r} not found (ids look like SC0001).', 'not_found')
        out = dict(e)
        out['links'] = links_for(site, snpcurate=route('snpcurate', e['id']),
                                 snpfunction=route('snpfunction', e['gene']) if e.get('gene') else None,
                                 snpgeo=route('snpgeo', e['gene']) if e.get('gene') else None)
        if e.get('gene'): out['links']['maizegdb_gene'] = MAIZEGDB_GENE + e['gene']
        out['marks'] = cu.get('marks')
        return out
    limit = arg_int(args, 'limit', 50, 1, CAPS['max_list'])
    mark = arg_enum(args, 'mark', ['gold', 'silver', 'white'])
    role = str(args.get('role') or '').strip().lower() or None
    gene = str(args.get('gene') or '').strip()
    q = norm_key(args.get('query') or args.get('trait') or '')
    rows = []
    for e in cu['entries']:
        if mark and e.get('mark') != mark: continue
        if role and (e.get('role') or '').lower() != role: continue
        if gene and gene.lower() not in ((e.get('gene') or '').lower(), (e.get('symbol') or '').lower()): continue
        if q:
            hay = ' '.join(str(e.get(k) or '') for k in ('symbol', 'gene', 'label', 'trait', 'note', 'published_as', 'role')).lower()
            if q not in hay: continue
        s = e.get('site') or {}
        rows.append({'id': e['id'], 'symbol': e.get('symbol'), 'gene': e.get('gene'), 'label': e.get('label'), 'mark': e.get('mark'),
                     'role': e.get('role'), 'variant_type': e.get('variant_type'), 'status': e.get('status'), 'trait': e.get('trait'),
                     'site': f"{e['chr']}:{e['pos']} {e['ref']}>{e['alt']}" if e.get('pos') else None,
                     'favorable': e.get('favorable'), 'b73': s.get('b73'), 'carriers_hom': (s.get('nHom')), 'carriers_het': s.get('nHet'),
                     'site_qc': s.get('qc'), 'tier': s.get('priority'), 'composite': s.get('combined'),
                     'link': (site or '') + route('snpcurate', e['id'])})
    return {'matched': len(rows), 'returned': min(limit, len(rows)), 'truncated': len(rows) > limit, 'entries': rows[:limit],
            'marks': cu.get('marks'), 'genes': sorted({r['gene'] for r in rows if r['gene']}),
            'links': links_for(site, snpcurate=route('snpcurate'))}


def t_alleles(args, site):
    ds = arg_dataset(args)
    if str(args.get('founders', '')).lower() in ('1', 'true', 'yes'):
        ids, unresolved = ['ZmG_' + f for f in family().get('namFounders', []) if 'ZmG_' + f in catalog()], []
    else:
        ids, unresolved = resolve_lines(args.get('lines'))
    sites = as_list(args.get('sites'))
    if not sites: raise QueryError('sites: registry ids (SC0001) and/or chr:pos[:ref:alt] sites, at most %d.' % CAPS['max_sites'])
    if len(sites) > CAPS['max_sites']: raise QueryError(f'sites: at most {CAPS["max_sites"]} per call.', 'too_large')
    cols, carried = [], {sid: 0 for sid in ids}
    for s in sites:
        a = {'curated': s} if re.match(r'^SC\d+$', s, re.I) else {'site': s}
        try:
            st, i, reg = locate_site(a, ds)
        except QueryError as ex:
            cols.append({'site': s, 'error': ex.args[0]}); continue
        gt = {sid: int(st.f[sid][i]) if sid in st.f else 3 for sid in ids}
        fav = (reg or {}).get('favorable')
        def carry(g):
            if fav not in ('ref', 'alt'): return None
            if g >= 3: return 'missing'
            want = 2 if fav == 'alt' else 0
            return 'yes' if g == want else ('het' if g == 1 else 'no')
        col = {'site': f"{st.chrom}:{int(st.pos[i])}:{st.f['REF'][i].decode()}:{st.f['ALT'][i].decode()}",
               'curated': ({'id': reg['id'], 'symbol': reg.get('symbol'), 'label': reg.get('label'), 'trait': reg.get('trait'), 'favorable': fav,
                            'favorable_note': reg.get('favorable_note'), 'mark': reg.get('mark')} if reg else None),
               'b73': GT_WORD.get(int(st.f[B73][i]) if B73 in st.f else 3, 'missing'),
               'genotypes': {sid: GT_WORD.get(g, 'missing') for sid, g in gt.items()},
               'favorable_carry': {sid: carry(g) for sid, g in gt.items()} if fav in ('ref', 'alt') else None}
        for sid, g in gt.items():
            if col['favorable_carry'] and col['favorable_carry'][sid] == 'yes': carried[sid] += 1
        if reg: col['link'] = (site or '') + route('snpcurate', reg['id'], 'lines=' + ','.join(ids))
        cols.append(col)
    cat = catalog()
    return {'dataset': ds, 'lines': [line_brief(cat[s]) for s in ids if s in cat], 'unresolved_lines': unresolved or None,
            'sites': cols, 'favorable_sites_per_line': carried,
            'legend': {'genotypes': 'ref = homozygous reference (B73 v5 allele), het, alt = homozygous alternate, missing',
                       'favorable_carry': 'yes = homozygous for the favorable allele of the registry entry, het, no, missing; null when the registry does not name one'},
            'links': links_for(site, snpcurate=route('snpcurate'))}


def similarity_row(focal, which):
    d = os.path.join(ROOT, 'distance', FAMILY)
    ids = [l.strip() for l in open(os.path.join(d, 'ids.txt')) if l.strip()]
    if focal not in ids: raise QueryError(f'{focal} is not in the genome-wide IBS set.', 'not_found')
    k = ids.index(focal)
    csvp = os.path.join(d, 'similarity.csv' if which == 'all' else 'similarity_snp.csv')
    if not os.path.isfile(csvp): raise QueryError('The genome-wide similarity matrix is not installed on this server.', 'unavailable')
    with open(csvp) as fh:
        for n, line in enumerate(fh):
            if n == k:
                vals = [float(x) for x in line.strip().split(',')]
                break
    return ids, vals


def t_similarity(args, site):
    mode = arg_enum(args, 'mode', ['genome', 'region'], 'region' if (args.get('region') or args.get('gene')) else 'genome')
    top = arg_int(args, 'top', 20, 1, CAPS['max_top'])
    which = arg_enum(args, 'sites', ['all', 'snp'], 'all')
    cat = catalog()
    if mode == 'genome':
        focal = one_line(args.get('line') or args.get('focal'))
        ids, vals = similarity_row(focal['id'], which)
        subset = None
        if args.get('lines'):
            subset, _ = resolve_lines(args['lines'])
        pairs = [(ids[j], v) for j, v in enumerate(vals) if ids[j] != focal['id'] and (subset is None or ids[j] in subset)]
        pairs.sort(key=lambda x: -x[1])
        fmt = lambda p: {'line': line_brief(cat[p[0]]) if p[0] in cat else {'id': p[0]}, 'ibs_similarity': round(p[1], 6)}
        ibs = json.load(open(os.path.join(ROOT, 'distance', FAMILY, 'ibs_summary.json'))) if os.path.isfile(os.path.join(ROOT, 'distance', FAMILY, 'ibs_summary.json')) else {}
        return {'mode': 'genome', 'focal': line_brief(focal), 'sites': which,
                'sites_compared': ibs.get('sites_all' if which == 'all' else 'sites_snp'),
                'most_similar': [fmt(p) for p in pairs[:top]], 'least_similar': [fmt(p) for p in pairs[-min(5, len(pairs)):]][::-1],
                'mean_similarity': round(float(np.mean([p[1] for p in pairs])), 6) if pairs else None,
                'definition': 'genome-wide identity-by-state: 1 - mean allele distance / 2 over sites called in both lines',
                'links': links_for(site, snpcompare=route('snpcompare'), snptree=route('snptree'))}
    r = region_or_gene(args)
    ds = arg_dataset(args)
    st = Store(r['chrom'], ds)
    a, b = st.span(r['start'], r['end'])
    if b <= a: raise QueryError(f'No sites in {r["label"]} ({ds}).', 'not_found')
    focal = one_line(args.get('line') or args.get('focal')) if (args.get('line') or args.get('focal')) else None
    if args.get('lines'):
        ids, unresolved = resolve_lines(args['lines'])
    else:
        ids, unresolved = list(st.samples), []
    if focal and focal['id'] not in ids: ids = [focal['id']] + ids
    if focal is None and len(ids) < 2: raise QueryError('Give a focal line, or at least two lines to compare in the region.')
    if focal is None and len(ids) > CAPS['max_lines']: raise QueryError(f'At most {CAPS["max_lines"]} lines for an all-pairs region comparison.', 'too_large')
    M = np.stack([st.f[s][a:b] if s in st.f else np.full(b - a, 3, dtype=np.int8) for s in ids]).astype(np.int16)
    called = M <= 2
    def sim(i, j):
        m = called[i] & called[j]
        n = int(m.sum())
        return (round(1.0 - float(np.abs(M[i][m] - M[j][m]).mean()) / 2.0, 6) if n else None), n
    out = {'mode': 'region', 'dataset': ds, 'region': {k: r[k] for k in ('chrom', 'start', 'end', 'label', 'gene', 'flank') if k in r},
           'sites_in_region': b - a, 'unresolved_lines': unresolved or None,
           'definition': 'identity-by-state in the region: 1 - mean |g1 - g2| / 2 over sites called in both lines (g: 0 ref, 1 het, 2 alt)',
           'links': links_for(site, snpcompare=route('snpcompare'), snpversity=region_route(r, ds))}
    if focal:
        k = ids.index(focal['id'])
        rows = []
        for j, s in enumerate(ids):
            if j == k: continue
            v, n = sim(k, j)
            rows.append({'line': line_brief(cat[s]) if s in cat else {'id': s}, 'ibs_similarity': v, 'sites_co_called': n})
        rows.sort(key=lambda x: -(x['ibs_similarity'] if x['ibs_similarity'] is not None else -1))
        out.update({'focal': line_brief(focal), 'lines_compared': len(rows), 'most_similar': rows[:top], 'least_similar': rows[-min(5, len(rows)):][::-1]})
    else:
        pairs = []
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                v, n = sim(i, j)
                pairs.append({'a': ids[i], 'b': ids[j], 'ibs_similarity': v, 'sites_co_called': n})
        pairs.sort(key=lambda x: -(x['ibs_similarity'] if x['ibs_similarity'] is not None else -1))
        out.update({'lines': [line_brief(cat[s]) for s in ids if s in cat], 'pairs': pairs[:max(top, 1) * 5] if len(pairs) > top * 5 else pairs,
                    'n_pairs': len(pairs)})
    return out


def gwas_hits(r, study, trait, significant, limit):
    man = gwas_manifest()
    studies = [g for g in man if (not study or study.lower() in (g['id'].lower(), g.get('trait', '').lower()) or study.lower() in g['id'].lower())
               and (not trait or trait.lower() in (g.get('trait', '') + ' ' + g.get('traitCategory', '')).lower())]
    if not studies: raise QueryError('No GWAS study matches; call snpmaize_gwas with no arguments to list them.', 'not_found')
    hits, by_study = [], Counter()
    pre = (r['chrom'][3:] + ',') if r else None
    for g in studies:
        p = os.path.join(ROOT, 'data', 'gwas', g['file'])
        if not os.path.isfile(p): continue
        thr = float(g.get('threshP') or 0)
        with open(p) as fh:
            head = fh.readline().strip().split(',')
            for line in fh:
                if pre and not line.startswith(pre): continue
                t = line.rstrip('\n').split(',')
                if len(t) < 9: continue
                bp = int(t[2])
                if r and (bp < r['start'] or bp > r['end']): continue
                try: pv = float(t[8])
                except ValueError: continue
                sig = pv <= thr if thr else None
                if significant and not sig: continue
                by_study[g['id']] += 1
                hits.append({'study': g['id'], 'trait': g.get('trait'), 'population': g.get('population'), 'chr': 'chr' + t[0], 'pos': bp, 'snp': t[1],
                             'a1': t[3], 'a2': t[4], 'freq_a1': float(t[5]) if t[5] else None, 'beta': float(t[6]) if t[6] else None,
                             'se': float(t[7]) if t[7] else None, 'p': pv, 'significant': sig})
    hits.sort(key=lambda h: h['p'])
    return {'matched': len(hits), 'by_study': dict(by_study), 'hits': hits[:limit], 'truncated': len(hits) > limit,
            'studies': [{'id': g['id'], 'trait': g.get('trait'), 'category': g.get('traitCategory'), 'measure': g.get('measure'), 'population': g.get('population'),
                         'publication': g.get('publication'), 'threshold_p': g.get('threshP'), 'method': g.get('gwasMethod')} for g in studies]}


def t_gwas(args, site):
    limit = arg_int(args, 'limit', CAPS['default_rows'], 1, CAPS['max_rows'])
    study = str(args.get('study') or '').strip() or None
    trait = str(args.get('trait') or '').strip() or None
    significant = str(args.get('significant_only', 'true')).lower() not in ('0', 'false', 'no')
    r = region_or_gene(dict(args, flank=args.get('flank', 10000 if args.get('gene') else 0)), need=False)
    if not r and not study and not trait:
        man = gwas_manifest()
        return {'studies': [{'id': g['id'], 'trait': g.get('trait'), 'category': g.get('traitCategory'), 'measure': g.get('measure'),
                             'population': g.get('population'), 'publication': g.get('publication'), 'publication_url': g.get('publicationUrl'),
                             'threshold_p': g.get('threshP'), 'threshold_method': g.get('thresholdMethod'), 'markers': g.get('totalMarkers'),
                             'method': g.get('gwasMethod'), 'link': (site or '') + route('snpgwas', g['id'])} for g in man],
                'note': 'Give a gene or region (and optionally study or trait) for hits; study alone lists its top hits genome-wide.',
                'links': links_for(site, snpgwas=route('snpgwas'))}
    if r:
        check_span(r['chrom'], r['start'], r['end'], CAPS['max_span_bp'] * 8)
    out = gwas_hits(r, study, trait, significant, limit)
    out.update({'region': ({k: r[k] for k in ('chrom', 'start', 'end', 'label', 'gene', 'flank') if k in r} if r else None),
                'significant_only': significant, 'returned': min(limit, out['matched']),
                'note': 'Positions are B73 v5. Significance uses each study\'s own SimpleM threshold; the GWAS panels (NAM, Ames) are not the SNPMaize lines.',
                'links': links_for(site, snpgwas=route('snpgwas', out['studies'][0]['id']) if len(out['studies']) == 1 else route('snpgwas'))})
    if r and r.get('gene'): out['links']['snpfunction'] = (site or '') + route('snpfunction', r['gene'])
    return out


TOOLS = {'snpmaize_datasets': t_datasets, 'snpmaize_lines': t_lines, 'snpmaize_line': t_line, 'snpmaize_variants': t_variants,
         'snpmaize_gene': t_gene, 'snpmaize_carriers': t_carriers, 'snpmaize_curated': t_curated, 'snpmaize_alleles': t_alleles,
         'snpmaize_similarity': t_similarity, 'snpmaize_gwas': t_gwas}


def run(tool, args, site=''):
    t0 = time.time()
    if tool not in TOOLS: raise QueryError(f'unknown tool {tool!r}', 'invalid_params')
    if not isinstance(args, dict): raise QueryError('arguments must be an object')
    out = TOOLS[tool](args, site)
    out['meta'] = {'release': family().get('label'), 'reference': 'B73 v5', 'seconds': round(time.time() - t0, 3), 'tool': tool}
    return out


def jsonable(o):
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray): return o.tolist()
    if isinstance(o, (set, tuple)): return list(o)
    if isinstance(o, bytes): return o.decode()
    raise TypeError(f'not JSON serialisable: {type(o).__name__}')


def main(argv):
    if argv[1:2] == ['--list']:
        print('\n'.join(TOOLS)); return 0
    if argv[1:]:
        req = {'tool': argv[1], 'args': json.loads(argv[2]) if len(argv) > 2 else {}, 'site': os.environ.get('SNPMAIZE_SITE', '')}
    else:
        req = json.loads(sys.stdin.read() or '{}')
    budget = float(os.environ.get('SNPQ_TIMEOUT') or 0)      # wall-time budget set by mcp.php, from process start
    def expire(*_):
        raise QueryError(f'The query took longer than {budget:g} s and was stopped. Ask for a smaller region or fewer lines.', 'timeout')
    try:
        if budget > 0:
            left = budget - (time.time() - T_START)
            if left <= 0: expire()
            signal.signal(signal.SIGALRM, expire)
            signal.setitimer(signal.ITIMER_REAL, left)
        out = run(req.get('tool'), req.get('args') or {}, req.get('site') or '')
        out['ok'] = True
    except QueryError as e:
        out = {'ok': False, 'error': {'code': e.code, 'message': str(e), **e.extra}}
    finally:
        if budget > 0: signal.setitimer(signal.ITIMER_REAL, 0)
    json.dump(out, sys.stdout, default=jsonable, ensure_ascii=False)
    sys.stdout.write('\n')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv))
