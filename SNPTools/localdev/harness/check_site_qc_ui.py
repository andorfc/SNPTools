#!/usr/bin/env python3
"""check_site_qc_ui.py <results.json> <site_root> -- Site QC as the pages show it.

run_scenarios.js (siteqcui) drives the pages with the full stores, once as is and once with
SNPTOOLS_SITEQC=0 for the PHP process. Recomputed here from the VCFs the pages' queries wrote:
(a) SNPVersity, eight gene intervals x 26 NAM lines: each row's class (from its NHET / NHOM, the
    panel-wide counts), the rows left under All sites / Hide flagged and no-carrier / Passing
    only, the "<n> flagged, <m> with no carrier ... of <N> sites." line, the Site QC column after
    MAF and the text of its cells.
(b) SNPFunction, the same genes, from the 933-line VCF its query wrote: every site with an
    entry for the gene (that gene's most severe SnpEff entry), carriers and class; "<v> of <n>
    sites vary"; the allele list (damaging alleles at usable sites), the "Flagged calls" and "No
    carrier in this release" groups, knockout lines (homozygous for a usable loss-of-function
    allele) and those only in flagged calls, the burden's usable-site count, and the banner
    (more than half of at least 5 protein-changing sites with a carrier are flagged). The page's
    DOM must show the same lists, labels and pills as its data.
(c) SNPImpact on su1: the rows under each choice and the hidden-count note. SNPFold on
    Zm00001eb406050: rows with "Usable alleles only" on and off. SNPGeo: a pill on every row.
(d) SNPTree, SNPMatrix, SNPCompare on chr6:91,593,082-91,793,082: sites used under both settings.
(e) The JS rule (Data.siteQc) agrees with the rule on all 3,721 (het, hom) pairs from 0 to 60.
(f) With SNPTOOLS_SITEQC=0: no Site QC column, select, note, pill or site control in SNPVersity,
    SNPImpact, SNPFold, SNPGeo, SNPTree, SNPMatrix or SNPCompare, and SNPFunction gives the same
    lists from its own counts.
Standard library only."""
import gzip, json, os, re, sys

res, root = sys.argv[1:3]
R = json.load(open(res))
U = R.get('siteqcui')
NAMES = ['NO_CARRIER', 'HET_ONLY', 'HET_EXCESS', 'HET_ELEVATED', 'PASS']
LABEL = {'PASS': 'pass', 'HET_ELEVATED': 'Het elevated', 'HET_EXCESS': 'Het excess', 'HET_ONLY': 'Het only', 'NO_CARRIER': 'No carrier'}
SEV = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond:
        fails.append(msg)


def rule(het, hom):
    if het + hom == 0: return 'NO_CARRIER'
    if hom == 0: return 'HET_ONLY'
    if het > hom: return 'HET_EXCESS'
    if 4 * het >= het + hom: return 'HET_ELEVATED'
    return 'PASS'


usable = lambda q: q in ('PASS', 'HET_ELEVATED')
flagged = lambda q: q in ('HET_ONLY', 'HET_EXCESS')
keep = {'all': lambda q: True, 'usable': lambda q: q is None or usable(q), 'pass': lambda q: q is None or q == 'PASS'}


def dose(cell):
    a = cell.split(':')[0].replace('|', '/').split('/')
    if len(a) < 2 or '.' in a or '' in a: return None
    return (a[0] != '0') + (a[1] != '0')


def read_vcf(v):
    p = os.path.join(root, v[2:] if v.startswith('./') else v)
    samples, rows = [], []
    with gzip.open(p, 'rt') as fh:
        for line in fh:
            if line.startswith('##'): continue
            t = line.rstrip('\n').split('\t')
            if line.startswith('#'): samples = t[9:]; continue
            I = dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in t[7].split(';') if kv)
            rows.append((int(t[1]), t[3], t[4], I, t[9:]))
    return samples, rows


def info_qc(I):
    return rule(int(I['NHET']), int(I['NHOM'])) if 'NHET' in I and 'NHOM' in I else None


def clean(tok):
    t = (tok or '').split(';')[0].replace('_', ' ').strip()
    return '' if t == '.' else t


def klass(effect):         # Data.impactClass
    e = re.sub(r'\s+', '_', (effect or '').lower())
    if re.search(r'missense|protein_altering|non[_-]?synonymous', e): return 'missense', False
    if re.search(r'stop_gained|nonsense|frameshift|start_lost|initiator|stop_lost', e): return 'lof', True
    if re.search(r'splice_(acceptor|donor)', e): return 'splice', True
    if 'splice' in e: return 'splice', False
    if re.search(r'inframe_(insertion|deletion)', e): return 'indel', False
    if re.search(r'synonymous|stop_retained', e): return 'syn', False
    return 'other', False


def entries(I):            # Data.annotationsOf: (gene, effect, impact, esm_here)
    G, T = I.get('GENEMODEL', '').split(','), I.get('TYPE', '').split(',')
    E = I.get('EFFECT', '').split(',')
    esm_at = next((i for i, x in enumerate(T) if 'missense' in x), -1)
    multi = len(G) > 1 and len(T) == len(G)
    if len(G) > 1 and len(T) != len(G): G, T, E = G[:1], T[:1], E[:1]
    return [(clean(G[i]) or '—', clean(T[i]) or 'intergenic',
             E[i].strip().upper() if i < len(E) and E[i].strip().upper() in SEV else 'MODIFIER',
             (i == esm_at) if multi else True) for i in range(len(G))]


def for_gene(I, gene):     # Data.rowForGene: the gene's most severe entry, the first on a tie
    best = None
    for e in entries(I):
        if e[0] == gene and (best is None or SEV[e[2]] > SEV[best[2]]): best = e
    return best


def num(v):
    return None if v in (None, '', '.') else float(v.split(',')[0])


def gene_function(vcf, gene):
    samples, rows = read_vcf(vcf)
    V = []
    for pos, ref, alt, I, cells in rows:
        e = for_gene(I, gene)
        if not e: continue
        cls, severe = klass(e[1])
        d = [dose(c) for c in cells]
        het, hom = d.count(1), d.count(2)
        pc, es = num(I.get('plantcad1_score')), (num(I.get('ESM1_score')) if e[3] else None)
        comb = round((pc + es) / 2, 2) if pc is not None and es is not None else (pc if pc is not None else es)
        V.append(dict(key=(pos, ref, alt), cls=cls, severe=severe, het=het, hom=hom, qc=rule(het, hom), combined=comb,
                      homIds={samples[i] for i, x in enumerate(d) if x == 2}))
    dmg = [v for v in V if v['cls'] == 'lof' or v['severe'] or (v['cls'] == 'missense' and v['combined'] is not None and v['combined'] <= -4)]
    ko = set().union(*[v['homIds'] for v in dmg if usable(v['qc']) and v['cls'] == 'lof']) if dmg else set()
    kof = set().union(*[v['homIds'] for v in dmg if flagged(v['qc']) and v['cls'] == 'lof']) - ko if dmg else set()
    pcv = [v for v in V if (v['cls'] in ('missense', 'lof', 'indel') or (v['cls'] == 'splice' and v['severe'])) and v['qc'] != 'NO_CARRIER']
    nf = sum(flagged(v['qc']) for v in pcv)
    return dict(n=len(V), nLines=len(samples), variable=sum(v['qc'] != 'NO_CARRIER' for v in V), used=sum(usable(v['qc']) for v in V),
                list={v['key'] for v in dmg if usable(v['qc'])}, flagged={v['key'] for v in dmg if flagged(v['qc'])},
                nocarrier={v['key'] for v in dmg if v['qc'] == 'NO_CARRIER'}, ko=len(ko), kof=len(kof),
                banner=bool(pcv) and nf / len(pcv) > 0.5 and len(pcv) >= 5, share=(nf, len(pcv)),
                qc={v['key']: v['qc'] for v in V}, het={v['key']: (v['het'], v['hom']) for v in V})


if not U:
    print('site QC in the pages: skipped (the full chr2/chr4/chr6/chr8/chr10 stores are not installed)')
    sys.exit(0)

# (e) the JS rule
diff = [p for p in U['jsRule'] if p[2] != rule(p[0], p[1])]
check(len(U['jsRule']) == 3721 and not diff, f"JS Data.siteQc agrees with the rule on all {len(U['jsRule'])} pairs {diff[:3]}")
check(U['defaults'] == {'versity': 'all', 'impact': 'usable', 'func': 'usable', 'fold': 'usable', 'geo': 'all', 'distance': 'usable'},
      f"Data.SITE_QC_DEFAULTS {U['defaults']}")

FUNC = {}
for mode in ('on', 'off'):
    O = U[mode]; on = mode == 'on'
    tag = '' if on else ' [SNPTOOLS_SITEQC=0]'
    # (a) SNPVersity
    for sym, v in O['versity'].items():
        _, rows = read_vcf(v['vcf'])
        qc = [info_qc(I) for _, _, _, I, _ in rows]
        cnt = {k: qc.count(k) for k in NAMES}
        filt = {c: sum(keep[c](q) for q in qc) for c in keep}
        hdr = v['headers']
        if on:
            ok = (v['qcRows'] == qc and all(q is not None for q in qc) and v['filter'] == filt and v['select'] and v['selected'] == 'all'
                  and hdr[hdr.index('MAF') + 1] == 'Site QC'
                  and v['note'] == f"{cnt['HET_ONLY'] + cnt['HET_EXCESS']:,} flagged, {cnt['NO_CARRIER']:,} with no carrier in this release, of {len(rows):,} sites."
                  and v['cellText'] == [LABEL[q] for q in qc[:400]])
        else:
            ok = (all(q is None for q in v['qcRows']) and all(q is None for q in qc) and 'Site QC' not in hdr and not v['select']
                  and v['note'] is None and v['filter'] == {c: len(rows) for c in keep} and not v['cellText'])
        check(ok, f"SNPVersity {sym}{tag}: {len(rows)} sites {cnt if on else ''} filters {v['filter']}")
    # (b) SNPFunction
    for sym, f in O['func'].items():
        d = f['data']
        exp = gene_function(f['vcf'], d['gene'])
        key = lambda lst: {(x[0], x[1], x[2]) for x in lst}
        ok = (d['nVariants'] == exp['n'] and d['nVariable'] == exp['variable'] and key(d['damaging']) == exp['list']
              and key(d['damagingFlagged']) == exp['flagged'] and key(d['damagingNoCarrier']) == exp['nocarrier']
              and d['koLines'] == exp['ko'] and d['koLinesFlagged'] == exp['kof'] and d['burden']['nUsed'] == exp['used']
              and (d['share']['flagged'], d['share']['withCarrier']) == exp['share']
              and all(x[6] == exp['qc'][(x[0], x[1], x[2])] and (x[4], x[5]) == exp['het'][(x[0], x[1], x[2])]
                      for x in d['damaging'] + d['damagingFlagged'] + d['damagingNoCarrier']))
        check(ok, f"SNPFunction {sym}{tag} data: {exp['variable']} of {exp['n']} vary, list {len(exp['list'])}, flagged {len(exp['flagged'])}, "
                  f"no carrier {len(exp['nocarrier'])}, KO {exp['ko']} (+{exp['kof']} flagged only), banner {exp['banner']} {exp['share']}")
        dom_list = [(r['variant'], r['qc']) for r in f['list']]
        want = lambda lst: [(x[3], LABEL[x[6]]) for x in lst]
        ok = (dom_list == want(d['damaging']) and [(r['variant'], r['qc']) for r in f['flagged']] == want(d['damagingFlagged'])
              and [(r['variant'], r['qc']) for r in f['nocarrier']] == want(d['damagingNoCarrier'])
              and f['variable'] == f"{exp['variable']:,} of {exp['n']:,} sites vary among the {exp['nLines']:,} lines"
              and f['burden'] == f"over {exp['used']:,} usable sites of {exp['n']:,}" and f['banner'] == exp['banner']
              and f['koLines'] == str(exp['ko']) and f['koFlagged'] == (f"{exp['kof']:,} more lines only in flagged calls" if exp['kof'] else None))
        check(ok, f"SNPFunction {sym}{tag} page: list {[r['variant'] for r in f['list']]}, flagged {[r['variant'] for r in f['flagged']]}, "
                  f"no carrier {[r['variant'] for r in f['nocarrier']]}, KO {f['koLines']} {f['koFlagged'] or ''}, banner {f['banner']}")
        FUNC.setdefault(sym, {})[mode] = {k: d[k] for k in ('nVariable', 'damaging', 'damagingFlagged', 'damagingNoCarrier', 'koLines', 'koLinesFlagged')}
    # (c) SNPImpact, SNPFold, SNPGeo
    im = O['impact']
    _, rows = read_vcf(im['vcf'])
    qc = [info_qc(I) for _, _, _, I, _ in rows]
    if on:
        hid = lambda c: {g: sum(not keep[c](q) and q in ks for q in qc) for g, ks in
                         (('flagged', ('HET_ONLY', 'HET_EXCESS')), ('nocarrier', ('NO_CARRIER',)), ('caution', ('HET_ELEVATED',)))}
        def note(c):
            h = hid(c)
            if not any(h.values()): return None
            parts = [f"{h['flagged']:,} flagged", f"{h['nocarrier']:,} no-carrier"] + ([f"{h['caution']:,} het-elevated"] if h['caution'] else [])
            lst = ', '.join(parts[:-1]) + ' and ' + parts[-1] if len(parts) > 2 else ' and '.join(parts)
            return f"{lst} variants hidden. Show all."
        ok = all(im[c]['rows'] == sum(keep[c](q) for q in qc) and im[c]['note'] == note(c) for c in ('usable', 'pass', 'all')) \
            and 'Site QC' in im['header'] and im['select'] == 5
        check(ok, f"SNPImpact su1: rows usable {im['usable']['rows']}, pass {im['pass']['rows']}, all {im['all']['rows']}; note \"{im['usable']['note']}\"")
    else:
        check(all(im[c]['rows'] == len(rows) and im[c]['note'] is None for c in ('usable', 'pass', 'all')) and 'Site QC' not in im['header'] and im['select'] == 4,
              f"SNPImpact su1{tag}: {len(rows)} rows under every choice, no Site QC select, column or note")
    fo = O['fold']
    fq = [x[3] for x in fo['data']['fold']]
    n_use = sum(usable(q) for q in fq)
    if on:
        check(fo['checkbox'] and fo['checked'] and fo['rows'] == n_use and fo['rowsAll'] == len(fq) and 'Site QC' in fo['header']
              and fo['count'] == (f"{n_use} of {len(fq)} coding variants shown" if n_use < len(fq) else f"{len(fq)} coding variants"),
              f"SNPFold Zm00001eb406050: {fo['rows']} of {fo['rowsAll']} with Usable alleles only ({fo['count']})")
    else:
        # SNPFold reads SNPFunction's full-panel counts, so it keeps the class without the sidecar
        check(fo['checkbox'] and fo['rows'] == n_use and fo['rowsAll'] == len(fq),
              f"SNPFold Zm00001eb406050{tag}: classes from SNPFunction's own counts, {fo['rows']} of {fo['rowsAll']} shown")
    ge = O['geo']
    check((ge['pills'] == ge['rows'] == ge['withQc'] and ge['label']) if on else (ge['pills'] == 0 and ge['withQc'] == 0 and not ge['label']),
          f"SNPGeo su1{tag}: {ge['pills']} pills on {ge['rows']} rows, label pill {ge['label']}")
    # (d) distance tools
    di = O['distance']
    _, rows = read_vcf(di['vcf'])
    qc = [info_qc(I) for _, _, _, I, _ in rows]
    n_all, n_use = len(rows), (sum(usable(q) for q in qc) if on else len(rows))
    t = lambda a, b: f"{a:,} of {b:,} sites used"
    if on:
        ok = (di['treeDefault'] == 'usable' and di['treeControl'] and di['mtxControl'] and di['cmpControl']
              and (di['tree_usable']['sites'], di['tree_all']['sites']) == (n_use, n_all) and di['tree_usable']['text'] == t(n_use, n_all)
              and di['mtx_usable'] == t(n_use, n_all) and di['mtx_all'] == t(n_all, n_all) and (di['cmp_usable'], di['cmp_all']) == (n_use, n_all))
    else:
        ok = (not di['treeControl'] and not di['mtxControl'] and not di['cmpControl']
              and di['tree_usable']['sites'] == di['tree_all']['sites'] == n_all and di['cmp_usable'] == di['cmp_all'] == n_all)
    check(ok, f"SNPTree / SNPMatrix / SNPCompare chr6:91,593,082-91,793,082{tag}: {n_use:,} of {n_all:,} sites used (usable only), {n_all:,} (all)")
    check(not O.get('consoleErrors'), f"no script errors{tag} {(O.get('consoleErrors') or [])[:2]}")

# (f) SNPFunction classifies from its own counts with or without the sidecar
for sym, m in FUNC.items():
    check(m.get('on') == m.get('off'), f"SNPFunction {sym}: same lists with SNPTOOLS_SITEQC=0")

print(f"site QC in the pages: {len(fails)} failed; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
