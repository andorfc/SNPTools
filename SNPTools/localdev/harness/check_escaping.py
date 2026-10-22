#!/usr/bin/env python3
"""check_escaping.py -- data strings reach the pages as text, never as markup (the escaping scenario of
run_scenarios.js): VCF alleles and INFO fields in SNPVersity and SNPImpact, the shared tooltip, a typed
PanEffect id. Static checks: the ontology view's outLink keeps non-http hrefs out, and no tool file
interpolates a data field into innerHTML without an escape (a list of known-risky patterns).

    python3 check_escaping.py <results.json> <root>

Exit 1 on a failure."""
import json, re, sys

res, root = sys.argv[1:3]
X = json.load(open(res)).get('escaping')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not X:
    print('FAIL no escaping scenario in results.json'); sys.exit(1)
v = X['versity']
check(not v['injected'] and v['altShown'] and v['effectShown'] and v['geneShown'] and v['subShown'] and v['effectOption'] == 1,
      f"SNPVersity: markers in ALT, TYPE, GENEMODEL and SUB are shown as text, none becomes an element ({v})")
t = X['tooltip']
check('T<img id="xss1">' in (t['altTip'] or '') and not t['altImg'] and (t['moreTip'] is None or ('<img' in t['moreTip'] and not t['moreImg'])),
      f"the shared tooltip shows its attribute as text: {t['altTip'][:50]!r}, image element {t['altImg']}; +1 genes tip {str(t['moreTip'])[:60]!r}")
i = X['impact']
check(not i['injected'] and not i['injectedDetail'] and i['geneShown'] and i['consShown'] and i['detailText'],
      f"SNPImpact: gene, consequence and variant text shown as text in the table and the detail ({i})")
p = X['paneffect']
check(not p['injected'] and (p['shown'] or not p['hosts']),
      f"PanEffect: a typed id with markup never becomes an element; shown as text where the shell has a message host ({p})")
src = open(f'{root}/js/snpfunction-ontology.js', encoding='utf-8').read()
check("/^(https?:\\/\\/|\\/)/i.test(h)" in src, "ontology outLink: only http(s) and site-relative hrefs become links")
core = open(f'{root}/js/core.js', encoding='utf-8').read()
check('tt.textContent=el.dataset.tt' in core and 'tt.innerHTML=el.dataset.tt;' not in core, 'core.js tooltip uses textContent')
# known-risky interpolations that must stay escaped
risky = [('js/snpversity.js', r'\$\{r\.effect\}'), ('js/snpversity.js', r'\$\{r\.sub\}'), ('js/snpimpact.js', r'>\$\{r\.gene\}<'),
         ('js/snpimpact.js', r'\$\{r\.consequence\}</span>'), ('js/snpfold.js', r'>\$\{v\.variant\}'), ('js/snpfold.js', r'\$\{v\.consequence\}</span>'),
         ('js/snpgwas.js', r'<div class="snp">\' \+ f\[1\]'), ('js/snpgwas.js', r'<td class="snpid">\' \+ f\[1\]'),
         ('js/pe/dom.js', r"' \+ main_id \+ '"), ('js/pe/main.js', r"' \+ main_id \+ '")]
for f, pat in risky:
    s = open(f'{root}/{f}', encoding='utf-8').read()
    check(not re.search(pat, s), f"{f}: no unescaped {pat}")
check(not X['errors'], f"no console errors ({len(X['errors'])})" + (': ' + '; '.join(e[:100] for e in X['errors'][:3]) if X['errors'] else ''))
print(f'check_escaping: {bad} failure(s)')
sys.exit(1 if bad else 0)
