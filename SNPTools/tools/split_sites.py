#!/usr/bin/env python3
"""split_sites.py -- one InterProScan sites file per gene for SNPFold, from data/results.sites.tsv.

    python3 tools/split_sites.py [--root DIR] [--src data/results.sites.tsv] [--out data/sites/by_gene]

The source is the header-less InterProScan "sites" TSV for every isoform of every gene (12 columns;
protein ids like Zm00001eb062040_P001). SNPFold used to fetch all ~55 MB of it on its first visit to show
one gene. This writes <out>/<gene>.tsv (every isoform's rows of that gene, in source order) and
<out>/MANIFEST.json ({format, built, source_bytes, genes, rows}); snpfold.js fetches a gene's file when
the manifest is there and falls back to the whole table when it is not. The output folder is a data
drop (not in git): rebuild it after replacing the source. Writes to a temporary folder, then renames.
"""
import argparse, datetime, json, os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ID_RE = re.compile(r'^(.*)_([A-Za-z]\d+)$')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--root', default=ROOT)
    ap.add_argument('--src', default='data/results.sites.tsv')
    ap.add_argument('--out', default='data/sites/by_gene')
    a = ap.parse_args()
    src = os.path.join(a.root, a.src); out = os.path.join(a.root, a.out)
    if not os.path.isfile(src): raise SystemExit(f'no source: {src}')
    tmp = out + '.tmp'
    if os.path.isdir(tmp): shutil.rmtree(tmp)
    os.makedirs(tmp)
    per = {}; rows = 0; skipped = 0
    with open(src, encoding='utf-8') as fh:
        for line in fh:
            if not line.strip(): continue
            f = line.rstrip('\n').split('\t')
            m = ID_RE.match(f[0]) if len(f) >= 12 else None
            if not m: skipped += 1; continue
            per.setdefault(m.group(1), []).append(line if line.endswith('\n') else line + '\n')
            rows += 1
    for gene, lines in per.items():
        with open(os.path.join(tmp, gene + '.tsv'), 'w', encoding='utf-8') as fh: fh.writelines(lines)
    man = {'format': 'snptools-sites-by-gene-1', 'built': datetime.date.today().isoformat(), 'source': a.src,
           'source_bytes': os.path.getsize(src), 'genes': len(per), 'rows': rows, 'skipped_lines': skipped}
    with open(os.path.join(tmp, 'MANIFEST.json'), 'w', encoding='utf-8') as fh: json.dump(man, fh, indent=1); fh.write('\n')
    if os.path.isdir(out): shutil.rmtree(out)
    os.rename(tmp, out)
    print(f'{rows:,} rows of {len(per):,} genes -> {os.path.relpath(out, a.root)}/ ({skipped} lines skipped)')


if __name__ == '__main__':
    main()
