#!/usr/bin/env python3
"""check_snpcurate.py <results.json> <site_root> -- every number SNPCurate shows, recomputed.

From data/curate/snpcurate.source.json, the variant stores (read through h5_to_vcf.py with
PYTHON_PATH: each site for all 933 lines, each missense entry's gene interval for one line, whose INFO
carries the panel-wide NHET / NHOM), js/zmgrin.catalog.js, data/traits/zmgrin2026.traits.json and
gff/genes_index.txt, independently of tools/build_snpcurate.py:
(a) each site entry: het / hom / missing calls, Site QC class (and the source's expect block), the
    allele frequency among called lines, MAF, MAXR2, the consequence and substitution of site_gene
    (else gene) read as the app reads it (that gene's most severe SnpEff entry; ESM scores only on the
    entry they were computed for), the seven scores, the combined score and priority
    (Data.impactPriority), the carriers, homozygous carriers per subpopulation and carriers per
    country with their totals, the trait by genotype, and for a missense entry its rank by each
    score among the gene's missense variants that have a carrier;
(b) each other entry: no site at its pos or codon positions;
(c) the marks (gold / outline / grey) and the counts the page states;
(d) the page as drawn (run_scenarios.js, curate): the data file, the table rows, the filters and
    sorts, every record's facts, scores and ranks, subpopulation, country and trait rows, the reason
    shown for an allele the release cannot show, the references as DOI links, and the three buttons
    of SC0001 (carriers to SNPVersity, gene to SNPFunction, site to SNPGeo).
Standard library only; the store reads run h5_to_vcf.py."""
import json, math, os, re, subprocess, sys, tempfile

res, root = sys.argv[1:3]
root = os.path.abspath(root)
K = json.load(open(res)).get('curate')
PY = os.environ.get('PYTHON_PATH') or 'python3'
SEV = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
QC_LABEL = {'PASS': 'pass', 'HET_ELEVATED': 'Het elevated', 'HET_EXCESS': 'Het excess', 'HET_ONLY': 'Het only', 'NO_CARRIER': 'No carrier'}
SCORES = [('pc1', 'plantcad1_score', 'PlantCAD1'), ('pc2', 'plantcad2_score', 'PlantCAD2'), ('evo2', 'evo2_score', 'Evo2'),
          ('esm1', 'ESM1_score', 'ESM1'), ('esm2', 'ESM2_score', 'ESM2'), ('esm3', 'ESM3_score', 'ESM3'), ('esmc', 'ESMC_score', 'ESM-C')]
STATUS = {'site': 'Site in this release', 'not_a_site': 'Not a site', 'not_locatable': 'Not locatable', 'structural': 'Structural variant',
          'unreliable': 'Unreliable calls', 'not_annotated': 'Gene not annotated'}
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond: fails.append(msg)


if not K:
    print('SNPCurate: skipped (the full stores are not installed)'); sys.exit(0)

src = json.load(open(f'{root}/data/curate/snpcurate.source.json', encoding='utf-8'))
cat_js = open(f'{root}/js/zmgrin.catalog.js', encoding='utf-8').read()
fam = json.loads(re.search(r'families\["zmgrin2026"\]\s*=\s*(\{.*?\});\n', cat_js, re.S).group(1))
CAT = {a['id']: a for p in fam['projects'] for g in p['groups'] for a in g['accessions']}
IDS = list(CAT)
TR = json.load(open(f'{root}/data/traits/zmgrin2026.traits.json', encoding='utf-8'))
GENES = {}
for line in open(f'{root}/gff/genes_index.txt'):
    t = [x.strip() for x in line.split('\t')]
    if not line.startswith('#') and len(t) >= 4: GENES[t[0]] = (t[1], int(t[2]), int(t[3]))
tmp = tempfile.mkdtemp(prefix='snpt_curate_')
ids_all = os.path.join(tmp, 'all.json'); json.dump(IDS, open(ids_all, 'w'))
ids_one = os.path.join(tmp, 'one.json'); json.dump(['ZmG_B73'], open(ids_one, 'w'))


def vcf(chrom, lo, hi, ids):
    """Rows of h5_to_vcf.py for the interval: (pos, ref, alt, info dict, {sample: dose or None})."""
    out = os.path.join(tmp, 'q.vcf')
    if os.path.exists(out): os.remove(out)
    p = subprocess.run([PY, os.path.join(root, 'h5_to_vcf.py'), f'{root}/hdf5/grin2026/zmgrin2026_{chrom}_impute.h5', out, str(lo), str(hi), ids],
                       capture_output=True, text=True, cwd=root)
    if not os.path.exists(out): return []
    rows, samples = [], []
    for line in open(out):
        if line.startswith('##'): continue
        t = line.rstrip('\n').split('\t')
        if line.startswith('#'): samples = t[9:]; continue
        I = dict(kv.split('=', 1) if '=' in kv else (kv, '') for kv in t[7].split(';') if kv)
        rows.append((int(t[1]), t[3], t[4], I, {s: dose(c) for s, c in zip(samples, t[9:])}))
    return rows


def dose(cell):
    a = cell.split(':')[0].replace('|', '/').split('/')
    if len(a) < 2 or '.' in a or '' in a: return None
    return (a[0] != '0') + (a[1] != '0')


def rule(het, hom):
    if het + hom == 0: return 'NO_CARRIER'
    if hom == 0: return 'HET_ONLY'
    if het > hom: return 'HET_EXCESS'
    if 4 * het >= het + hom: return 'HET_ELEVATED'
    return 'PASS'


def clean(t):
    t = re.split(r'[;,]+', t or '')[0].replace('_', ' ').strip()
    return '' if t == '.' else t


def klass(effect):
    e = re.sub(r'\s+', '_', (effect or '').lower())
    if re.search(r'missense|protein_altering|non[_-]?synonymous', e): return 'missense', False, 'Missense'
    if re.search(r'stop_gained|nonsense', e): return 'lof', True, 'Stop gained'
    if 'frameshift' in e: return 'lof', True, 'Frameshift'
    if re.search(r'start_lost|initiator', e): return 'lof', True, 'Start lost'
    if 'stop_lost' in e: return 'lof', True, 'Stop lost'
    if re.search(r'splice_(acceptor|donor)', e): return 'splice', True, 'Splice site'
    if 'splice' in e: return 'splice', False, 'Splice region'
    if 'inframe_insertion' in e: return 'indel', False, 'In-frame insertion'
    if 'inframe_deletion' in e: return 'indel', False, 'In-frame deletion'
    if re.search(r'synonymous|stop_retained', e): return 'syn', False, 'Synonymous'
    for k, lab in (('intron', 'Intron'), ('5_prime_utr', "5' UTR"), ('3_prime_utr', "3' UTR"), ('upstream', 'Upstream'), ('downstream', 'Downstream'), ('intergenic', 'Intergenic')):
        if k in e: return 'other', False, lab
    return 'other', False, clean(effect) or 'Other'


def entry(I, gene):
    G = I.get('GENEMODEL', '')
    if ',' not in G:
        imp = max((x.strip().upper() for x in re.split(r'[;,]+', I.get('EFFECT', '')) if x.strip().upper() in SEV), key=SEV.get, default='MODIFIER')
        e = (clean(G) or '—', clean(I.get('TYPE')) or 'intergenic', imp, clean(I.get('SUB')), True)
        return e if e[0] == gene else None, e
    G, T, E, U = (I.get(k, '').split(',') for k in ('GENEMODEL', 'TYPE', 'EFFECT', 'SUB'))
    esm_at = next((i for i, t in enumerate(T) if 'missense' in t), -1)
    es = [(clean(G[i]) or '—', clean(T[i]) or 'intergenic', E[i].strip().upper() if i < len(E) and E[i].strip().upper() in SEV else 'MODIFIER',
           clean(U[i]) if i < len(U) else '', i == esm_at) for i in range(len(G))]
    best = None
    for e in es:
        if e[0] == gene and (best is None or SEV[e[2]] > SEV[best[2]]): best = e
    return best, es[0]


def num(v):
    try: return None if v in (None, '', '.') else float(v.split(',')[0])
    except ValueError: return None


def jsround(x): return int(math.floor(x + 0.5))
def pct(a, b): return f" ({jsround(100 * a / b)}%)" if b else ''
def sc1(v): return '—' if v is None else ('+' if v > 0 else '') + f'{v:.1f}'
def fmtm(v): return '—' if v is None else (f'{jsround(v):,}' if abs(v) >= 100 else f'{v:.2f}')


def priority(cls, severe, comb, level):
    if severe: return 'TOP'
    if cls in ('missense', 'indel'):
        if comb is not None: return 'TOP' if comb <= -7 else 'HIGH' if comb <= -4 else 'MODERATE' if comb <= -1 else 'LOW'
        return level if level in ('HIGH', 'MODERATE') else 'LOW'
    return 'HIGH' if level == 'HIGH' else 'LOW'


sub_tot, cty_tot = {}, {}
for a in CAT.values():
    sub_tot[a.get('subpop') or 'Unknown'] = sub_tot.get(a.get('subpop') or 'Unknown', 0) + 1
    cty_tot[a.get('country') or 'Unknown'] = cty_tot.get(a.get('country') or 'Unknown', 0) + 1
D = {e['id']: e for e in K['data']['entries']}
EXP = {}
for e in src['entries']:
    x = {'mark': 'grey'}
    if e['status'] != 'site':
        hits = [p for p in ([e['pos']] if e.get('pos') else []) + list(e.get('codon_positions') or []) if vcf(e['chr'], p, p, ids_one)]
        check(not hits, f"{e['id']} {e['symbol']} ({e['status']}): no site at {e.get('pos') or e.get('codon_positions') or 'its interval (nothing to test)'}")
        EXP[e['id']] = x; continue
    x['mark'] = 'gold' if e['kind'] == 'causal' and e['evidence'] == 'validated' and not e.get('needs_review') else 'outline'
    rows = [r for r in vcf(e['chr'], e['pos'], e['pos'], ids_all) if (r[0], r[1], r[2]) == (e['pos'], e['ref'], e['alt'])]
    if len(rows) != 1:
        check(False, f"{e['id']}: {len(rows)} sites for {e['chr']}:{e['pos']} {e['ref']}>{e['alt']}"); continue
    pos, ref, alt, I, G = rows[0]
    het = sum(v == 1 for v in G.values()); hom = sum(v == 2 for v in G.values()); miss = sum(v is None for v in G.values())
    target = e.get('site_gene') or e['gene']
    ent, first = entry(I, target)
    a = ent or first
    cls, severe, label = klass(a[1])
    sc = {k: (num(I.get(key)) if (a[4] or not k.startswith('esm')) else None) for k, key, _ in SCORES}
    comb = round((sc['pc1'] + sc['esm1']) / 2, 2) if sc['pc1'] is not None and sc['esm1'] is not None else (sc['pc1'] if sc['pc1'] is not None else sc['esm1'])
    homs = sorted(s for s, v in G.items() if v == 2); hets = sorted(s for s, v in G.items() if v == 1)
    sub, cty = {}, {}
    for s in homs: sub[CAT[s].get('subpop') or 'Unknown'] = sub.get(CAT[s].get('subpop') or 'Unknown', 0) + 1
    for s in homs + hets: cty[CAT[s].get('country') or 'Unknown'] = cty.get(CAT[s].get('country') or 'Unknown', 0) + 1
    x.update(nHet=het, nHom=hom, nMiss=miss, qc=rule(het, hom), af=round((het + 2 * hom) / (2 * (len(G) - miss)), 4),
             maf=num(I.get('MAF')), maxr2=num(I.get('MAXR2')), consequence=label, consClass=cls, sub=a[3] or None, impact=a[2],
             scores=sc, combined=comb, priority=priority(cls, severe, comb, a[2]), hom=homs, het=hets,
             subpop=sorted(((k, sub.get(k, 0), n) for k, n in sub_tot.items()), key=lambda t: (-t[1], -t[2], t[0])),
             country=sorted(((k, v, cty_tot[k]) for k, v in cty.items()), key=lambda t: (-t[1], t[0])))
    if e.get('expect'):
        check((het, hom, x['qc']) == (e['expect']['het'], e['expect']['hom'], e['expect']['site_qc']), f"{e['id']} {e['symbol']}: {het} het, {hom} hom, {x['qc']} = the source's expect block")
    if e.get('grin_descriptor'):
        d = TR['dictionary'][e['grin_descriptor']]
        val = lambda s: (TR['samples'].get(s) or {}).get('traits', {}).get(e['grin_descriptor'])
        if d['scale'] == 'numeric':
            t = {}
            for key, code in (('ref', 0), ('alt', 2)):
                xs = [val(s)[1] for s, v in G.items() if v == code and val(s) and val(s)[1] is not None]
                t[key] = (len(xs), sum(xs) / len(xs) if xs else None)
            x['trait'] = ('numeric', d.get('unit') or '', t)
        else:
            x['trait'] = ('coded', None, [(code, lab, sum(1 for s in G if val(s) and str(val(s)[1]) == str(code)),
                                         sum(1 for s, v in G.items() if v == 2 and val(s) and str(val(s)[1]) == str(code)))
                                        for code, lab in d['codes'].items() if any(val(s) and str(val(s)[1]) == str(code) for s in G)])
    if cls == 'missense' and ent:
        gc, glo, ghi = GENES[target]
        pool = []
        for p2, r2, a2, I2, _ in vcf(gc, glo, ghi, ids_one):
            if int(I2.get('NHET', 0)) + int(I2.get('NHOM', 0)) == 0: continue
            e2, _f = entry(I2, target)
            if e2 and klass(e2[1])[0] == 'missense':
                pool.append({k: (num(I2.get(key)) if (e2[4] or not k.startswith('esm')) else None) for k, key, _ in SCORES})
        x['rank'] = {k: ([1 + sum(1 for p in pool if p[k] is not None and p[k] < sc[k]), sum(1 for p in pool if p[k] is not None)] if sc[k] is not None else None)
                     for k, _, _ in SCORES}
        x['pool'] = len(pool)
    EXP[e['id']] = x
    s = D[e['id']]['site']
    ok = ((s['nHet'], s['nHom'], s['nMiss'], s['qc'], s['af'], s['maf'], s['maxr2'], s['consequence'], s['consClass'], s['sub'], s['impact'], s['combined'], s['priority'])
          == (het, hom, miss, x['qc'], x['af'], x['maf'], x['maxr2'], label, cls, x['sub'], a[2], comb, x['priority'])
          and s['scores'] == sc and s['carriers'] == {'hom': homs, 'het': hets}
          and [(t['name'], t['hom'], t['total']) for t in s['subpop']] == x['subpop'] and [(t['name'], t['carriers'], t['total']) for t in s['country']] == x['country']
          and (s.get('rank') == x.get('rank')) and (s.get('rankPool') == x.get('pool')))
    if 'trait' in x:
        kind, unit, t = x['trait']
        if kind == 'numeric':
            ok = ok and all(s['trait'][k]['n'] == t[k][0] and (t[k][1] is None or abs(s['trait'][k]['mean'] - t[k][1]) < 1e-3) for k in ('ref', 'alt'))
        else:
            ok = ok and [(c['code'], c['label'], c['lines'], c['alt']) for c in s['trait']['codes']] == t
    check(ok, f"{e['id']} {e['symbol']:<15} data: {het}/{hom} {x['qc']}, {label}{' ' + x['sub'] if x['sub'] else ''}, priority {x['priority']} ({comb}), "
              f"{len(homs) + len(hets)} carriers{', rank pool ' + str(x['pool']) if 'pool' in x else ''}")

# (c) marks and counts
marks = {m: sorted(i for i, x in EXP.items() if x['mark'] == m) for m in ('gold', 'outline', 'grey')}
check(all(D[i]['mark'] == x['mark'] for i, x in EXP.items()) and len(D) == len(src['entries']) == 21,
      f"marks: gold {marks['gold']}, outline {len(marks['outline'])}, grey {marks['grey']}")
check([l['badge'] for l in K['legend']] == ['cur-badge cur-gold', 'cur-badge cur-outline', 'cur-badge cur-grey']
      and len({l['glyph'] for l in K['legend']}) == 3 and all(f"({len(marks[m])})" in l['text'] for m, l in zip(('gold', 'outline', 'grey'), K['legend'])),
      f"legend: three marks, three glyphs, counts {[len(marks[m]) for m in ('gold', 'outline', 'grey')]}")

# (d) the page as drawn
check(K['nav'] == ['SNPCuratenew'] and not K.get('consoleErrors'), f"nav item {K['nav']}; no script errors {K.get('consoleErrors')}")
# the table opens sorted by mark (gold, outline, grey), then gene name (case-insensitive), then id
MARK_ORDER = {'gold': 0, 'outline': 1, 'grey': 2}
order = [e['id'] for e in sorted(src['entries'], key=lambda e: (MARK_ORDER[EXP[e['id']]['mark']], e['symbol'].lower(), e['id']))]
tbl_ok = [r['id'] for r in K['table']] == order
for r in K['table']:
    e = next(x for x in src['entries'] if x['id'] == r['id']); x = EXP[r['id']]
    cells = r['cells']   # gene, change, trait, status, carriers, site qc, priority, reference
    want = [f"{e['symbol']}{e['gene']}", e['label'], e['trait'], STATUS[e['status']],
            f"{x['nHet']:,} / {x['nHom']:,}" if 'nHet' in x else '—', QC_LABEL.get(x.get('qc'), ''), x.get('priority', '—')]
    if cells[:7] != want or r['badge'] != f"cur-badge cur-{x['mark']}": tbl_ok = False; print('   row', r['id'], cells[:7], want)
check(tbl_ok, f"table: {len(K['table'])} rows by mark then gene ({', '.join(order[:4])}, ...) with mark, gene, change, trait, status, carriers, Site QC, priority")
F = K['filters']
carr = lambda i: EXP[i]['nHet'] + EXP[i]['nHom'] if 'nHet' in EXP[i] else None
src_order = [e['id'] for e in src['entries']]     # a header sort starts from the source order; ties by id
want_c = sorted([i for i in src_order if carr(i) is not None], key=lambda i: (-carr(i), src_order.index(i))) + [i for i in src_order if carr(i) is None]
want_g = sorted(order, key=lambda i: (next(e['symbol'] for e in src['entries'] if e['id'] == i).lower(), i))
check(F['gold'] == [i for i in order if i in marks['gold']] and F['grey'] == [i for i in order if i in marks['grey']] and F['site'] == [i for i in order if 'nHet' in EXP[i]]
      and F['kernel'] == [i for i in order for e in src['entries'] if e['id'] == i and 'kernel' in (e['trait'] + ' ' + e['symbol'] + ' ' + e['gene']).lower()] and len(F['kernel']) > 1
      and F['byCarriers'] == want_c and [i.lower() for i in F['byGene']] == [i.lower() for i in want_g],
      f"filters (gold {len(F['gold'])}, grey {len(F['grey'])}, sites {len(F['site'])}, 'kernel' {len(F['kernel'])}) and sorts (carriers, gene)")
for e in src['entries']:
    r, x = K['records'][e['id']], EXP[e['id']]
    refs_ok = len(r['refs']) == len(e['refs']) and r['dois'] == [f"https://doi.org/{ref['doi']}" for ref in e['refs'] if ref.get('doi')]
    review_ok = (r['review'] == f"To confirm: {e['needs_review']}") if e.get('needs_review') else r['review'] is None
    if 'nHet' not in x:
        ok = r['why'] is not None and r['why'].startswith(STATUS[e['status']] + '.') and refs_ok and review_ok and not r['buttons'] \
             and (not e.get('tagged_by') or any(e['tagged_by'] in p for p in r['note']))
        check(ok, f"{e['id']} {e['symbol']:<15} record: \"{(r['why'] or '')[:70]}\"{' tagged by ' + e['tagged_by'] if e.get('tagged_by') else ''}")
        continue
    facts = r['facts']
    exp_qc = f"Site QC{QC_LABEL[x['qc']]} {x['nHet']:,} het · {x['nHom']:,} hom · {x['nMiss']:,} missing"
    exp_pr = f"Priority{x['priority']} combined {'—' if x['combined'] is None else ('+' if x['combined'] > 0 else '') + format(x['combined'], '.2f')}"
    exp_af = f"Allele frequency{x['af'] * 100:.1f}% among called lines · MAF {('—' if x['maf'] is None else format(x['maf'], 'g'))}"
    scores = [[lab, sc1(x['scores'][k])] + ([f"rank {x['rank'][k][0]} of {x['rank'][k][1]}"] if x.get('rank') and x['rank'].get(k) else []) for k, _, lab in SCORES]
    subs = [f"{n}{h:,} of {t:,}{pct(h, t)}" for n, h, t in x['subpop'] if h]
    ctys = [f"{n}{c:,} of {t:,}{pct(c, t)}" for n, c, t in x['country']]
    trait_ok = True
    if 'trait' in x:
        kind, unit, t = x['trait']
        u = f' {unit}' if unit else ''
        if kind == 'numeric':
            want_t = [f"Homozygous for the alternate allele{fmtm(t['alt'][1])}{u} (n = {t['alt'][0]:,})",
                      f"Homozygous for the reference allele{fmtm(t['ref'][1])}{u} (n = {t['ref'][0]:,})"]
        else:
            want_t = [f"{lab} ({code}){alt:,} of {n:,}{pct(alt, n)}" for code, lab, n, alt in t]
        trait_ok = r['trait'] == want_t and r['traitHead'].startswith(TR['dictionary'][e['grin_descriptor']]['name'])
    ok = (facts[2] == exp_qc and facts[4] == exp_pr and facts[3].startswith(exp_af) and r['scores'] == scores and r['subpop'] == (subs or ['No homozygous carrier'])
          and r['country'] == ctys and trait_ok and refs_ok and review_ok and r['badge'] == f"cur-badge cur-{x['mark']}"
          and r['buttons'] == ['Carriers in SNPVersity', 'Gene in SNPFunction', 'On the map'])
    if not ok:
        for k, g, w in (('qc', facts[2], exp_qc), ('priority', facts[4], exp_pr), ('af', facts[3], exp_af), ('scores', r['scores'], scores),
                        ('subpop', r['subpop'], subs), ('country', r['country'][:3], ctys[:3]), ('trait', r['trait'], 'trait ok' if trait_ok else '')):
            if g != w and not (k == 'af' and str(g).startswith(w)): print('  ', e['id'], k, g, '| want', w)
    check(ok, f"{e['id']} {e['symbol']:<15} record: {exp_qc[7:]}, {exp_pr[8:]}, {len(subs)} subpopulations, {len(ctys)} countries"
              + (f", trait {e['grin_descriptor']}" if 'trait' in x else ''))

# SC0001's three buttons
x = EXP['SC0001']; g = GENES['Zm00001eb175150']
check(K['toVersity'] == {'tool': 'snpversity', 'selected': sorted(x['hom'] + x['het']), 'chr': 'chr4', 'start': min(g[1], 46648374), 'end': max(g[2], 46648374)},
      f"SC0001 Carriers in SNPVersity: {len(K['toVersity']['selected'])} lines on {K['toVersity']['chr']}:{K['toVersity']['start']}-{K['toVersity']['end']}")
check(K['toFunction'] == {'tool': 'snpfunction', 'gene': 'Zm00001eb175150'}, f"SC0001 Gene in SNPFunction: {K['toFunction']}")
check(K['toMap']['tool'] == 'snpgeo' and K['toMap']['rows'] == [[46648374, 'C', 'G']] and K['toMap']['accs'] == len(IDS),
      f"SC0001 On the map: SNPGeo with {K['toMap']['rows']} for {K['toMap']['accs']} lines")
# (e) the marks in the other tools (run_scenarios.js, curateMarks), recomputed from the source
CM = json.load(open(res)).get('curateMarks')
if CM:
    W, WO = CM['with'], CM['without']
    GLY = {'gold': '\u2605', 'outline': '\u2606', 'grey': '\u25cb'}
    SRC = src['entries']
    def in_iv(e, chrom, lo, hi):
        return e['chr'] == chrom and ((e.get('pos') is not None and lo <= e['pos'] <= hi) or any(lo <= p <= hi for p in e.get('codon_positions') or [])
                                      or (e.get('start') is not None and e['start'] <= hi and e['end'] >= lo))
    v = W['versity']
    exp_line = [e for e in SRC if in_iv(e, 'chr4', 46647932, 46652896)]
    exp_badge = [{'pos': f"{e['pos']:,}", 'cls': f"cur-badge cur-{EXP[e['id']]['mark']}", 'id': e['id'], 'glyph': GLY[EXP[e['id']]['mark']]} for e in exp_line if e['status'] == 'site']
    check(v['badges'] == exp_badge and v['line'] == 'Curated alleles in this interval: ' + ' · '.join(f"{GLY[EXP[e['id']]['mark']]} {e['symbol']} {e['label']}" for e in exp_line)
          and v['opened'] == {'tool': 'snpcurate', 'id': exp_badge[0]['id']} and v['lines'] == 33,
          f"SNPVersity tga1 x {v['lines']} lines: {[(b['pos'], b['cls'].split('-')[-1], b['id']) for b in v['badges']]}; the badge opens {v['opened']['id']}")
    for k, g in (('su1', 'Zm00001eb174590'), ('Bx13', 'Zm00001eb116010'), ('DGAT1-2', 'Zm00001eb277490')):
        want = [e for e in SRC if g in (e['gene'], e.get('site_gene'))]
        got = W['func'][k]['block']
        ok = [b['id'] for b in got] == [e['id'] for e in want] and all(b['cls'] == f"cur-badge cur-{EXP[e['id']]['mark']}" for b, e in zip(got, want))
        for b, e in zip(got, want):
            x = EXP[e['id']]
            if 'nHet' in x: ok = ok and f"{e['chr']}:{e['pos']:,} · {x['nHet']:,} het / {x['nHom']:,} hom {QC_LABEL[x['qc']]} {x['priority']}" in b['text']
            else: ok = ok and b['text'].endswith('Not a site in this release.' if e['status'] == 'not_a_site' else '.')
        if k == 'DGAT1-2':   # the Phe469 insertion is listed although the allele list does not rank it
            ok = ok and [115184121, 'T', 'TGAA'] not in W['func'][k]['list']
        check(ok, f"SNPFunction {k} curated block: {[(b['id'], b['cls'].split('-')[-1]) for b in got]}")
    im = W['impactTga1']
    check(im['filter'] == 'usable' and [(c['variant'], c['id'], c['cls']) for c in im['curated']] == [('p.N6K', 'SC0001', 'cur-badge cur-gold')],
          f"SNPImpact chr4:46,647,932-46,652,896, Hide flagged and no-carrier: {[(c['variant'], c['id']) for c in im['curated']]}")
    check(W['geo']['badges'] == ['SC0001'] and W['geo']['label'], f"SNPGeo tga1: badge {W['geo']['badges']} in the table and beside the label")
    t = [c for c in W['impactSu1']['curated'] if c['id'] == 'SC9901']
    check(W['testEntry'] and W['impactSu1']['filter'] == 'usable' and len(t) == 1 and t[0]['qc'] == 'Het only' and t[0]['cls'] == 'cur-badge cur-outline'
          and t[0]['variant'] == 'p.G627W' and 'SC9901' in W['funcSu1Test'],
          f"test-only su1 G627W (HET_ONLY): kept by SNPImpact's default filter with {t[0]['qc'] if t else None} and an outline badge; listed by SNPFunction")
    check(W['foldTest']['ring'] == ['SC9902'] and W['foldTest']['table'] == ['SC9902'],
          f"test-only {W['foldTest']['variant']} of Zm00001eb406050: SNPFold ring {W['foldTest']['ring']} and table badge {W['foldTest']['table']}")
    src_txt = open(f'{root}/data/curate/snpcurate.source.json', encoding='utf-8').read() + open(f'{root}/js/snpcurate.data.js', encoding='utf-8').read()
    check('SC99' not in src_txt, 'the test-only entries stayed in the loaded copy (not in the source or the data file)')
    check(not WO['versity']['badges'] and WO['versity']['line'] is None and all(not WO['func'][k]['block'] for k in WO['func'])
          and not WO['impactTga1']['curated'] and not WO['geo']['badges'] and not WO['geo']['label'] and WO['curatePage']['empty']
          and not WO['consoleErrors'] and not W['consoleErrors'],
          f"without js/snpcurate.data.js: no badge, no block, SNPCurate says '{WO['curatePage']['text']}', no script error {W['consoleErrors'] + WO['consoleErrors']}")
for f in os.listdir(tmp): os.remove(os.path.join(tmp, f))
os.rmdir(tmp)
print(f"SNPCurate: {len(src['entries'])} entries, {sum('nHet' in x for x in EXP.values())} sites; {len(fails)} failed; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
