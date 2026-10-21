#!/usr/bin/env python3
"""check_feedback.py <results.json> <site_root> -- Contact us and GWAS Explorer's "Add your results".

run_scenarios.js (R.feedback) drives the page; the Jira collector script is external, so the
harness serves it empty and calls its triggerFunction with a stand-in dialog:
(a) the masthead has Help & FAQ and Contact us (class feedback-form, aria-label), each with an icon:
    index.html hides the labels at 640 px and below, so a button without one would be blank there;
(b) the collector config (ATL_JQ_PAGE_PROPS['883299e6']) has a triggerFunction; once called, a
    click on Contact us opens the dialog once and does not follow href="#";
(c) Help & FAQ is the only button marked while the Help page is open, and nothing stays marked
    after going to a tool;
(d) GWAS Explorer: "Add your results to this tool!" beside "Export results…" opens the dialog,
    whose table lists Chr, bp and p as Required and A1, A2, Freq, b and se as Requested; Close
    closes it; its contact link (class feedback-form) closes it and opens the collector once.
Standard library only."""
import json, re, sys

res, root = sys.argv[1:3]
F = json.load(open(res)).get('feedback')
fails = []


def check(cond, msg):
    print(('ok   ' if cond else 'FAIL ') + msg)
    if not cond: fails.append(msg)


if not F:
    print('feedback: FAIL no results'); sys.exit(1)
h = F['header']
check([x['href'] for x in h] == ['#help', '#'] and 'feedback-form' in h[1]['cls'] and h[1]['aria'] == 'Contact us'
      and [x['label'] for x in h] == ['Help & FAQ', 'Contact us'], f"(a) masthead buttons {[(x['href'], x['label']) for x in h]}")
html = open(f'{root}/index.html', encoding='utf-8').read()
hides = re.search(r'@media \(max-width:640px\)\{[^}]*\}[^@]*?\.mg-help \.mg-lbl\{display:none\}', html, re.S)
check(bool(hides) and all(x['svg'] for x in h), "(a) every masthead button has an icon (labels are hidden at 640 px and below)")
check(F['props'] and F['props']['trigger'] == 'function', f"(b) collector config {F['props']}")
hc = F['headerClick']
check(hc['shown'] == 1 and hc['prevented'] and hc['hash'] == '#snpversity', f"(b) Contact us opens the dialog once, the page's own route untouched {hc}")   # href="#" must not navigate
check(F['helpOn'] == [True, False], f"(c) on the Help page only Help & FAQ is marked {F['helpOn']}")
check(F['afterGo'] == [False, False], f"(c) nothing marked after going to a tool {F['afterGo']}")
g = F['gwas']
want = [['Chr', 'Required'], ['bp', 'Required'], ['p', 'Required'], ['A1', 'Requested'], ['A2', 'Requested'],
        ['Freq', 'Requested'], ['b', 'Requested'], ['se', 'Requested']]
check(g['btn'] == 'Add your results to this tool!' and g['exportBtn'] and g['closedAtFirst'] and g['opened'] and g['doneCloses'],
      "(d) GWAS Explorer: the button opens the dialog, Close closes it")
check(g['rows'] == want, f"(d) submission columns {[r[0] for r in g['rows']]}")
c = g['contact']
check(c['closes'] and c['shown'] == 1 and c['prevented'] and 'feedback-form' in c['cls'], f"(d) its contact link closes it and opens the collector {c}")
check(not F.get('consoleErrors'), f"no script errors {F.get('consoleErrors')}")
print(f"feedback: {len(fails)} failed; mismatches={len(fails)}")
sys.exit(1 if fails else 0)
