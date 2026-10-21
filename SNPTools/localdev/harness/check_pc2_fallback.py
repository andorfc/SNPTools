#!/usr/bin/env python3
"""check_pc2_fallback.py -- the PlantCAD2 columns fall back to the 512-bp score when a store has no
one-pass score (Data.parseVcf, Data.compositeScore, scoreColorValue): the pc2fallback scenario of
run_scenarios.js.

    python3 check_pc2_fallback.py <results.json> <root>

Exit 1 on a failure."""
import json, sys

res, root = sys.argv[1:3]
F = json.load(open(res)).get('pc2fallback')
bad = 0
def check(ok, msg):
    global bad
    print(('ok   ' if ok else 'FAIL ') + msg)
    if not ok: bad += 1

if not F:
    print('FAIL no pc2fallback scenario in results.json'); sys.exit(1)
a = F['without']
check(a['pc2op'] == a['pc2'] == -2.4 and a['is512'] is True and a['flag'] is True,
      f"a store without plantcad2_onepass_score: the PlantCAD2 column holds the 512-bp score ({a['pc2op']}), marked pc2opIs512, Data.pc2Fallback() true")
check(a['compN'] == 3 and a['compN4'] is None,
      f"the composite leaves the stand-in out: {a['compN']} of 4 scores for a missense site with PlantCAD1, ESM1 and ESM2")
check(a['color'] == -2.4 and a['note'].startswith('PlantCAD2 shows the masked 512-bp score'),
      f"the colour value is the 512-bp score itself ({a['color']}) and the note names the fallback")
b = F['with']
check(b['pc2op'] == -4.6 and b['pc2'] == -2.4 and not b['is512'] and b['flag'] is False and b['compN'] == 4 and b['note'] == '' and b['color'] != -4.6,
      f"a store with the one-pass score: PlantCAD2 = {b['pc2op']} (512-bp {b['pc2']}), 4 of 4 scores, no note, colour on the one-pass map")
check(F['pageNotes'] == {'versity': True, 'impact': True},
      f"the note shows in SNPVersity and SNPImpact when the fallback is on: {F['pageNotes']}")
check(F['pageNotesOff'] == {'versity': False, 'impact': False}, f"and not when it is off: {F['pageNotesOff']}")
check(not F['errors'], f"no console errors ({len(F['errors'])})")
print(f'check_pc2_fallback: {bad} failure(s)')
sys.exit(1 if bad else 0)
